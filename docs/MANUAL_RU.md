# Руководство MBS-fast

MBS-fast считает рассеяние света несферическими ледяными частицами: сначала трассируются геометрические лучи, затем для выходных пучков считается дифракция в приближении Physical Optics. Код поддерживает CPU-расчеты через OpenMP/MPI и CUDA-ускорение дифракции на GPU.

## Содержание

| Раздел | Что описано |
|---|---|
| [Сборка](#сборка) | CPU, GPU, варианты точности, EPYC Zen, AVX2/AVX-512 |
| [Быстрый старт](#быстрый-старт) | Минимальные команды запуска |
| [Физическая модель](#физическая-модель) | Трассировка, диады/Jones, дифракция, Mueller |
| [Частицы и размер](#частицы-и-размер) | `-p`, `--pf`, `--rs`, `--k_eq`, показатель преломления |
| [Ориентационные сетки](#ориентационные-сетки) | `--oldauto`, `--random`, Sobol, Euler, `--pole` |
| [Сетка рассеяния](#сетка-рассеяния) | `--grid`, `--tgrid`, `--nphi`, веса углов |
| [GPU и multi-GPU](#gpu-и-multi-gpu) | CUDA backend, FFT, атомики, сканы по размерам |
| [Все флаги](#все-флаги) | Все флаги, которые парсит текущий `src/main.cpp` |
| [Переменные окружения](#переменные-окружения) | Основные production/debug переключатели |
| [Выходные файлы](#выходные-файлы) | Формат `.dat` и диагностические файлы |

## Сборка

### Требования

| Компонент | Для чего нужен | Комментарий |
|---|---|---|
| GCC >= 9 или Clang >= 14 | CPU и host-часть CUDA | На EPYC обычно используется GCC |
| OpenMP | Многопоточность CPU | Для GCC линкуется `libgomp` |
| MPI | CPU split build | `cpu/Makefile` использует `mpicxx`, если он есть |
| CUDA toolkit | GPU split build | Нужны `nvcc`, `libcudart`, `libcufft` |
| NVIDIA driver | GPU запуск | Драйвер должен поддерживать собранный `sm_XX` |

Рекомендуемый путь сборки - split build:

```bash
make -C cpu -j          # CPU MPI/OpenMP
make -C gpu -j          # GPU FP64 без fast-math по умолчанию
```

Объектные файлы CPU лежат в `cpu/build/`, CUDA - в `gpu/build/`. Общая физика и CLI остаются в `src/`, поэтому CPU и GPU версии не расходятся в разные кодовые базы.

### CPU

```bash
# MPI + OpenMP CPU binary
make -C cpu -j

# Один MPI rank, 64 OpenMP потока
OMP_NUM_THREADS=64 cpu/bin/mbs_po_mpi --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 64 --close -o out_cpu

# Четыре MPI rank, по 16 потоков на rank
mpirun -np 4 cpu/bin/mbs_po_mpi --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 16 --close -o out_cpu_mpi
```

Debug/help сборка:

```bash
make -C cpu debug
cpu/bin/mbs_po_mpi_debug --help-debug
```

### GPU

```bash
# Контрольный GPU-бинарный файл по умолчанию: FP64 без fast-math
PATH=/usr/local/cuda/bin:$PATH make -C gpu -j

# Явные варианты
make -C gpu fp32      -j   # gpu/bin/mbs_po_gpu_float
make -C gpu fp64      -j   # gpu/bin/mbs_po_gpu_double
make -C gpu fp32_fast -j   # gpu/bin/mbs_po_gpu_float_fast
make -C gpu fp64_fast -j   # gpu/bin/mbs_po_gpu_double_fast
```

`gpu/Makefile` работает и без `nvcc` в `PATH`: компилятор берётся из
`$CUDA_PATH/bin/nvcc` (по умолчанию `CUDA_PATH=/usr/local/cuda`). Допустимая
версия host GCC читается из заголовков установленной CUDA, после чего
автоматически выбирается самый новый совместимый компилятор, например
`g++-13` для CUDA 12.6. При необходимости выбор задаётся явно через
`CUDA_HOST_CXX=/path/to/g++`. Если системный GCC новее поддерживаемого CUDA,
достаточно установить совместимую версию рядом; менять системный `gcc` не надо.

| Target | Binary | CUDA точность | Fast math | Когда использовать |
|---|---|---:|---:|---|
| `make -C gpu` или `make -C gpu fp64` | `gpu/bin/mbs_po_gpu_double` | FP64 | нет | Контрольный и основной режим |
| `make -C gpu fp32` | `gpu/bin/mbs_po_gpu_float` | FP32 | нет | FP32 после проверки ошибки |
| `make -C gpu fp32_fast` | `gpu/bin/mbs_po_gpu_float_fast` | FP32 | да | Экспериментальный fast math; сначала измерить и проверить |
| `make -C gpu fp64_fast` | `gpu/bin/mbs_po_gpu_double_fast` | FP64 | да | Ускоренный FP64 после проверки |
| `make -C gpu double_debug` | `gpu/bin/mbs_po_gpu_double_debug` | FP64 | да | `--help-debug` и диагностика |

FP32 применяется к хранению геометрии дифракционных пучков, сумм Джонса и
результата Мюллера. Трассировка видимости и топологии, оптические пути,
поглощение, чувствительные к сокращению моменты многоугольников и вычисление
фазы остаются FP64. Команда `--version` выводит фактический профиль хранения,
фазы, математических функций и архитектуры GPU.

На RTX 3080 Ti выполнен контрольный тест невыпуклой поглощающей частицы:
`k_eq=20`, 64 ориентации в редуцированном SO(3), `N_phi=720`, `N_theta=360`.
Медианы времени по трём запускам:

| Профиль | Время | Ускорение относительно FP64 | Глобальная взвешенная L2-ошибка Мюллера | Взвешенная L2-ошибка M11 |
|---|---:|---:|---:|---:|
| Точный FP64 | 4,92 с | 1,00 раза | эталон | эталон |
| Точный смешанный FP32 | 3,45 с | 1,43 раза | `1,00e-6` | `9,99e-7` |
| FP32 с fast math | 3,63 с | 1,36 раза | `1,02e-6` | `1,01e-6` |

Поэтому для потребительских Ampere разумный начальный профиль скорости —
точный смешанный FP32. Fast math автоматически не выбирается: в этом тесте он
оказался медленнее и сильнее зависит от узкой интерференционной структуры.
FP64 остаётся режимом по умолчанию и обязательным эталоном при смене частицы,
показателя преломления, размера или сетки.

В split GPU build CUDA backend включен по умолчанию. Флаг `--gpu` можно писать, но он не обязателен. Флаг `--cpu` принудительно запускает CPU backend внутри GPU-capable бинарника.

Архитектура GPU определяется через `nvidia-smi`. При необходимости задавайте вручную:

```bash
make -C gpu fp32 -j GPU_ARCH=86   # Ampere: точный смешанный FP32
make -C gpu fp32 -j GPU_ARCH=89   # Ada: начать с точного смешанного FP32
```

### EPYC Zen5, AVX-512 и AVX2

Makefile берет CPU-флаги из `scripts/detect_arch_flags.sh` через `ARCH_FLAGS`. Для расчета на той же машине обычно лучше:

```bash
make -C cpu clean
make -C cpu -j ARCH_FLAGS="-march=native -mtune=native"

make -C gpu clean
make -C gpu double_fast -j ARCH_FLAGS="-march=native -mtune=native"
```

| Цель | Флаги GCC/Clang | Комментарий |
|---|---|---|
| EPYC Zen2 | `ARCH_FLAGS="-march=znver2 -mtune=znver2"` | AVX2/FMA |
| EPYC Zen3 | `ARCH_FLAGS="-march=znver3 -mtune=znver3"` | AVX2/FMA |
| EPYC Zen4 | `ARCH_FLAGS="-march=znver4 -mtune=znver4"` | AVX-512 при поддержке компилятора |
| EPYC Zen5 | `ARCH_FLAGS="-march=znver5 -mtune=znver5"` | Нужен GCC/Clang, который знает `znver5` |
| Portable AVX2 | `ARCH_FLAGS="-O3 -mavx2 -mfma"` | Запускается на AVX2 машинах |
| Явный AVX-512 | `ARCH_FLAGS="-O3 -mavx512f -mavx512dq -mavx512cd -mavx512bw -mavx512vl"` | Только если `lscpu` показывает эти флаги |

Флага `AVX12` у GCC/Clang нет. Если имелся в виду Zen5 с AVX-512, используйте `-march=znver5`. Если компилятор старый, обновите его или задайте явные `-mavx512*` флаги после проверки `lscpu`.

## Быстрый старт

```bash
# GPU double, полный диапазон 0..180 градусов
gpu/bin/mbs_po_gpu_double_fast --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 64 --close -o hex_gpu_double

# То же, но CPU backend из GPU бинарника
gpu/bin/mbs_po_gpu_double_fast --po --cpu -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 12 --oldauto 2 \
    --grid 0 180 600 180 --threads 64 --close -o hex_cpu_from_gpu_bin

# Частица из файла, масштабирование через k_eq
gpu/bin/mbs_po_gpu_float_fast --po --pf particle.dat --k_eq 58.81 \
    --ri 1.6 0.002 -w 1.064 -n 14 --oldauto 2 --pole \
    --grid 0 180 600 180 --threads 16 --close -o particle_keq

# Быстрая проверка двух точек: 0 и 180 градусов
gpu/bin/mbs_po_gpu_double_fast --po -p 1 125.9 78.09 \
    --ri 1.3116 0 -w 0.532 -n 1 --oldauto 2 --pole \
    --grid 0 180 600 1 --threads 64 --close -o check_0_180
```

## Физическая модель

### Общий pipeline

| Этап | Где в коде | Результат |
|---|---|---|
| Геометрия частицы | `src/particle`, `src/main.cpp` | Грани, нормали, площадь, объем, симметрия |
| Трассировка лучей | `src/scattering`, `src/tracer`, `src/Splitting.cpp` | Выходные пучки: направление, полигон, площадь, оптический путь, Jones |
| PO дифракция | `HandlerPO::ApplyDiffraction*`, `Handler::DiffractIncline*`, CUDA kernels | Jones амплитуда пучка в дальнем поле |
| Когерентная сумма | `HandlerPO::AddToMueller`, CUDA fused Mueller kernels | Сумма Jones матриц по пучкам для одной ориентации |
| Jones -> Mueller | `src/math/Mueller.cpp` | Реальная матрица Mueller 4x4 |
| Усреднение ориентаций | `HandlerPOTotal`, `TracerPOTotal` | Финальная матрица для случайно ориентированных частиц |

### Jones матрицы лучей

Каждый пучок несет комплексную матрицу Jones `beam.J`:

```text
J = [ Jvv  Jvh ]
    [ Jhv  Jhh ]
```

При отражениях/преломлениях `Splitting.cpp` домножает `beam.J` на коэффициенты Fresnel в локальном вертикальном/горизонтальном базисе луча. В PO режиме это амплитуда, поэтому в обычном когерентном режиме сначала суммируются Jones матрицы пучков, и только затем считается Mueller.

### Диады и поворот в базис рассеяния

Перед добавлением пучка в направление наблюдения локальный Jones базис надо перевести в базис детектора. Это делают `HandlerPO::RotateJones` и `RotateJonesFast`. Они строят 2x2 матрицу проекций через скалярные произведения между:

| Вектор | Смысл |
|---|---|
| `beam.Direction()` | Направление пучка после трассировки |
| `direction` | Требуемое направление рассеяния |
| `vf` | Forward reference direction |
| `info` polarization vectors | Предвычисленный поляризационный базис пучка |

В `ApplyDiffractionFast` вклад пучка имеет вид:

```text
J_beam(theta, phi) = F_edge(theta, phi) * R_out(theta, phi) * F_n(theta, phi)
```

| Множитель | Код | Смысл |
|---|---|---|
| `F_edge` | `DiffractInclineFast`, `DiffractIncline`, `DiffractInclineAbs` | Скалярный edge-интеграл Kirchhoff для полигона пучка |
| `R_out` | `RotateJones` или `KarczewskiJones` | Поворот/проекция Jones базиса в базис рассеяния |
| `F_n` | `ComputeFnJones` | Fresnel/Jones поправка для нормали пучка и направления |

### Флаг Karczewski для поляризации

`--karczewski` заменяет стандартную проекцию базиса `RotateJones` на
`HandlerPO::KarczewskiJones`. Он меняет только выходной поляризационный
множитель `R_out`; трассировка, площади пучков, optical paths, Fresnel
коэффициенты, дифракционная амплитуда и сетка рассеяния остаются теми же.

Это не универсальный флаг "лучше", а инструмент для проверки поляризационной
конвенции. В коде ветка помечена как experimental: она использует Karczewski
матрицу в aperture-frame, но не полностью повторяет всю coordinate pipeline
GOAD. Включать ее имеет смысл, когда сравниваются поляризационно-чувствительные
элементы Mueller или reference calculation использует Karczewski convention.

Поперечные нормы в этой ветви вычисляются устойчиво через `hypot`, все базисы
проверяются перед нормировкой. В вырожденной полюсной конфигурации программа
автоматически применяет основной `RotateJones`; `NaN` из экспериментальной
ветви в итоговую матрицу не пропускается. Это защищает численный расчёт, но не
делает неполную реализацию конвенции GOAD физическим эталоном.

Ожидаемый эффект:

| Величина | Что должно происходить |
|---|---|
| `M11` | Не должен меняться. `M11 = 0.5 * ||J||_F^2`, а Karczewski/default rotations сохраняют одну и ту же Frobenius norm. |
| `M12`, `M22` | Обычно менее чувствительны, чем нижний поляризационный блок, но могут сдвигаться при другой convention базиса. |
| `M33`, `M34`, `M44` | Могут заметно измениться, потому что меняется структура Jones-матрицы в базисе рассеяния. |
| Интегральные величины, завязанные на интенсивность | Нельзя считать улучшенными только из-за `--karczewski`; нужно смотреть сами Mueller elements. |

Практическое правило: если `M11` совпадает, а расходятся только `M33/M34/M44`,
запусти тот же case с `--karczewski` и без него. Если `--karczewski` приближает
результат к reference, проблема скорее в поляризационной convention. Если
меняется `M11`, причина не в этом флаге, а в геометрии, diffraction, cutoff или
normalization.

### Дифракционный интеграл Kirchhoff

Для каждого выходного пучка освещенная апертура представлена полигоном. PO дальнее поле считается как интеграл по апертуре:

```text
F(q) = integral_A exp(i k q.r) dA
```

где `A` - полигон пучка, `k = 2*pi/lambda`, `q` - разность между направлением пучка и направлением наблюдения. В коде интеграл считается через границы полигона, а не через sampling площади. Для ребра от `a` до `b` вклад выражается через разность комплексных экспонент и фазовый знаменатель; для малых знаменателей используются устойчивые sinc-пределы.

| Функция | Назначение |
|---|---|
| `Handler::DiffractIncline` | CPU scalar edge integral без поглощения |
| `Handler::DiffractInclineFast` | Оптимизированный CPU edge integral с предвычисленными ребрами |
| `Handler::DiffractInclineAbs` | Вариант с поглощением |
| `GpuDiffraction.cu` | CUDA реализация тех же вычислений |

CPU и GPU ветки должны использовать одну геометрию, фазовую конвенцию и сетку theta. Отличия обычно дают точность (`float`/`double`), `--use_fast_math`, FFT-интерполяция, порядок редукции/атомики и специальные случаи полюсов/endpoint.

### Когерентная и некогерентная сумма

По умолчанию:

```text
J_total(theta, phi) = sum_beams J_beam(theta, phi)
M(theta, phi)       = Mueller(J_total(theta, phi))
```

С флагом `--incoh`:

```text
M(theta, phi) = sum_beams Mueller(J_beam(theta, phi))
```

Преобразование Jones -> Mueller реализовано в `src/math/Mueller.cpp`: элементы 4x4 строятся из билинейных комбинаций `|S_i|^2`, `Re(S_i conj(S_j))`, `Im(S_i conj(S_j))`.

Когерентная сумма выполняется только между лучевыми путями одной частицы при одной фиксированной ориентации. Для ансамбля случайно и независимо ориентированных частиц программа сначала строит `Mueller(J_total)` для каждой ориентации, затем усредняет матрицы Mueller с ориентационными весами. Складывать Jones-матрицы разных ориентаций нельзя без координат частиц и относительных фаз падающего поля. Поэтому устаревший `--coherent-orientations` отключен: его прежняя реализация фактически давала обычное некогерентное усреднение ориентаций.

Это также означает, что модель одной частицы не воспроизводит когерентное обратное рассеяние среды, возникающее из интерференции взаимно обратных многократных путей между разными частицами. Для такой задачи нужен многочастичный решатель с положениями частиц и межчастичным распространением поля.

Отсечения `--beam-cutoff*` и `--trace-cutoff*` удаляют амплитуды до когерентной суммы и потому могут исказить перекрестные члены, даже если удаляемые интенсивности малы. Обратное направление следует проверять повторным расчетом с `--cutoff-profile safe` или `--cutoff-profile off`.

### Усреднение ориентаций и `--pole`

Ориентационные режимы задают веса по beta/gamma. На точных beta-полюсах все gamma эквивалентны, поэтому `--pole` считает одну gamma и умножает вес. В текущей версии при `--pole` используются beta endpoints, чтобы точка beta=0 или beta=pi действительно присутствовала в сетке, а не заменялась midpoint.

### Веса scattering grid

Для uniform theta grid выводится `Nth + 1` строк. Колонка `2pi*dcos` - вес кольца телесного угла:

```text
dOmega(theta_j) = 2*pi * (cos(theta_left) - cos(theta_right))
```

Для endpoint строк используется половинная ячейка, обрезанная границей диапазона. Для полной сетки `0..180` сумма всех `2pi*dcos` должна быть `4*pi`.

## Частицы и размер

| Флаг | Аргументы | Описание |
|---|---:|---|
| `-p` | `TYPE L D [extra]` | Built-in частица. Нужно указать ровно один источник: `-p` или `--pf`. |
| `--pf` | `FILE` | Частица из файла. |
| `--rs` | `SIZE` | Масштабировать file particle до `Dmax = SIZE`. |
| `--k_eq` | `X` | Масштабировать так, чтобы `k_eq = 2*pi*r_eq/lambda`. |
| `--ri` | `Re Im` | Комплексный показатель преломления. |
| `-w` | `LAMBDA` | Длина волны в микрометрах. |
| `-n` | `N` | Максимальная глубина внутренних отражений/преломлений. |

Координаты файловой частицы и вся внутренняя геометрия трассировки хранятся в
двойной точности. Загрузчик сдвигает координаты относительно не зависящего от
порядка центра ограничивающего параллелепипеда и удаляет только относительный к
масштабу шум текстовой сериализации. Глобальный перенос меняет лишь общую
оптическую фазу; каноническая форма сохраняет матрицу Мюллера, малые рёбра и
порядок граней при больших абсолютных координатах.

Типы `-p`:

| Type | Форма | Параметры |
|---:|---|---|
| 1 | Hexagonal column/plate | `L D` |
| 2 | Bullet | `L D` |
| 3 | Bullet rosette | `L D [cap]` |
| 4 | Droxtal | `L D extra` |
| 10 | Concave hexagonal | `L D concavity` |
| 12 | Hexagonal aggregate | `L D count` |
| 999 | Built-in aggregate | `extra` |

Масштабирование через equivalent size:

```text
r_eq = (3 V / (4*pi))^(1/3)
k_eq = 2*pi*r_eq/lambda
scale = (k_eq_target * lambda / (2*pi)) / r_eq_original
```

## Ориентационные сетки

| Режим | Аргументы | Когда применять | Комментарий |
|---|---:|---|---|
| `--oldauto` | `DIV` | Основной production режим | Шаг сетки связан с diffraction-limited angular scale; типично `2`, `4`, `8`. |
| `--random` | `Nb Ng` | Ручная beta/gamma сетка | Использует symmetry-reduced domain. |
| `--fixed` | `BETA GAMMA` | Отладка одной ориентации | Углы в градусах. |
| `--orientfile` | `FILE` | Пользовательские ориентации | Одна пара `beta_deg gamma_deg` в градусах на строку. |
| `--sobol` | `N` | Quasi-random average | Хорош для сходимости сканов. |
| `--sobol_seed` | `N S` | Sobol/Owen с seed | Повторяемые проверки сходимости. |
| `--sobol_ring` | `Nb Ng` | Sobol beta + uniform gamma rings | Гибридная сетка. |
| `--so3_quat` | `N` | Изотропное усреднение SO(3) | Hammersley в области симметрии; лабораторный alpha интегрируется по phi рассеяния. |
| `--so3_full_quat` | `N` | Независимая проверка SO(3) | Прямые кватернионы Shoemake по полной группе без сокращения по симметрии. |
| `--hammersley` | `N` | Hammersley orientations | Experimental/debug. |
| `--lattice` | `N` | Rank-1 lattice | Experimental/debug. |
| `--lattice_z` | `N Z` | Rank-1 lattice с явным generator | Experimental/debug. |
| `--euler_quad` | `Nb Ng` | Gauss по cos(beta), periodic gamma | Высокий порядок квадратуры. |
| `--euler_adapt` | `Nb NgMax` | Adaptive gamma count | Меньше работы около beta poles. |
| `--montecarlo` | `N` | Псевдослучайные ориентации | Базовая Monte Carlo проверка. |
| `--adaptive` | `EPS` | Автоподбор числа Sobol orientations | Удваивает число ориентаций до цели. |
| `--auto` | `EPS` | Auto theta, phi, orientations | Удобный режим. |
| `--autofull` | `EPS` | Auto `n`, theta, phi, orientations | Более полный и дорогой поиск. |
| `--oldautofull` | `EPS` | Autofull + oldauto final grid | Когда нужна регулярная финальная сетка. |

Модификаторы:

| Флаг | Аргументы | Описание |
|---|---:|---|
| `--ring_points` | `N` | Точек на diffraction ring для oldauto estimates. |
| `--mirror_gamma` | none | Для проверенной зеркально-симметричной частицы: половина диапазона gamma и восстановление второго вклада с учетом четности Стокса. Поддерживается `--so3_quat`. |
| `--so3_mirror_audit` | none | При четном `--so3_quat N` явно трассировать вложенные отраженные пары для проверки относительно `N/2 --mirror_gamma`. Только контрольный режим. |
| `--sym` | `Sb Sg` | Override symmetry: beta range `pi/Sb`, gamma range `2*pi/Sg`. |
| `--b` | `B1 B2` | Диапазон beta для `--random`, градусы. |
| `--g` | `G1 G2` | Диапазон gamma для `--random`, градусы. |
| `--maxorient` | `N` | Верхняя граница adaptive orientations. |
| `--chunk` | `N` | Chunk size по ориентациям/gamma. |
| `--coh_orient` | none | Отключенный устаревший режим; когерентность между разными ориентациями не определена без относительных фаз. |
| `--pole` | none | Одна gamma на точных beta poles. |
| `--owen_avg` | `K` | Усреднение `K` Owen seeds в `--autofull`. |
| `--owen_seeds` | `S...` | Явный список Owen seeds. |

### SO(3) с учетом симметрии

Для изотропного ансамбля мера Хаара в принятом соглашении углов имеет вид

```text
dR = sin(beta) d(beta) d(gamma) d(alpha).
```

Режим `--so3-quaternion N` использует фундаментальную область частицы
`0 <= beta <= beta_sym`, `0 <= gamma < gamma_sym`. Для `i=0,...,N-1` узлы
строятся как

```text
u_i = (i + 1/2) / N
v_i = radical_inverse_base_2(i)
beta_i  = acos(1 - (1 - cos(beta_sym)) u_i)
gamma_i = gamma_sym v_i
w_i = 1/N,  sum_i w_i = 1.
```

Для парной проверки зеркальной реконструкции к четному полному числу узлов
добавляется `--so3-mirror-audit`. При `K=N/2` используются те же формулы с
`i=0,...,K-1`, `u_i=(i+1/2)/K`, `gamma_i=(gamma_sym/2)v_i`, после чего явно
трассируются `(beta_i,gamma_i)` и `(beta_i,gamma_sym-gamma_i)` с весом `1/N`.
Эта сетка точно вложена в `--so3-quaternion K --mirror-gamma`. Без контрольного
флага обычная полнодоменная последовательность Хаммерсли не меняется.

Следовательно, равномерно распределен `cos(beta)`, а не сам beta. Флаг
`--sym Sb Sg` задает `beta_sym=pi/Sb`, `gamma_sym=2*pi/Sg`; его допустимо
использовать только при реальной симметрии геометрии.

С флагом `--mirror-gamma` вычисляется только диапазон
`0 <= gamma < gamma_sym/2`. Вторая половина восстанавливается до квадратуры по
alpha:

```text
M_full(theta,phi) = [M_half(theta,phi)
                    + P M_half(theta,-phi) P] / 2,
P = diag(1,1,-1,-1).
```

Это преобразование допустимо только при наличии у всей оптической частицы
соответствующей плоскости отражения, включая видимость граней и назначенные
материалы, а также при численной инвариантности выбранного трассировщика ФО
относительно такого отражения. По одному файлу частицы программа не может
доказать это свойство. Для каждого нового класса частиц нужен контрольный
расчет с полным диапазоном gamma. `N` остается числом реально трассируемых
ориентаций beta/gamma. Для строгой парной проверки надо сравнивать
`--so3-quaternion 4096 --so3-mirror-audit` и `--so3-quaternion 2048
--mirror-gamma`. Тогда базовые узлы совпадают, а разность измеряет численную
ошибку зеркальной реконструкции, а не шум двух разных выборок Хаммерсли.
Обычную полнодоменную сетку дополнительно сравнивают с независимым плотным
результатом; точку обратного рассеяния проверяют отдельно от общей нормы L2.

Угол alpha отдельно не трассируется. Поворот вокруг падающего луча эквивалентен
изменению азимута рассеяния phi, поэтому для каждой строки theta вычисляется

```text
M_avg(theta) = (1/Nphi) sum_phi M(theta,phi) L(-phi),
```

где `L` поворачивает базис Стокса и содержит `cos(2 phi)`, `sin(2 phi)` в блоке
Q/U. Матрица умножается справа, поскольку меняется базис падающей поляризации.
В направлениях theta=0 и theta=pi азимутальный базис вырожден: там усреднение
выполняется без `L`, после чего накладывается точная полюсная форма матрицы.
Это тот же проверенный код, который использует обычный variable-phi режим.

Необходимо `Nphi >= 2`. Рекомендуемый запуск:

```bash
mbs_po --method po ... --so3-quaternion 4096 \
    --scattering-diffraction-sampling 2 --latitude-phi-grid
```

`--so3-full-quaternion N` оставлен для независимой проверки. В нем `N` - число
полных трехмерных поворотов; в быстром режиме `N` - число реально трассируемых
ориентаций beta/gamma, а alpha считается квадратурой по phi.

## Сетка рассеяния

| Флаг | Аргументы | Описание |
|---|---:|---|
| `--grid` | `T1 T2 Nphi Nth` | Uniform theta range от `T1` до `T2`; выводит `Nth + 1` theta строк. |
| `--grid` | `R Nphi Nth` | Backscatter cone радиуса `R` градусов. |
| `--tgrid` | `FILE` | Неравномерная theta сетка, градусы, одна строка на theta. |
| `--auto_tgrid` | `EPS` | Adaptive theta grid через bisection. |
| `--auto_phi` | none | Автоматический выбор `Nphi`. |
| `--nphi` | `N` | Override phi count, максимальный приоритет. |
| `--diffraction-limit-grid` | `FACTOR` | Вычислить равномерные шаги theta и экваториальный phi из дифракционной оценки. |
| `--latitude-phi-grid` | none | Уменьшать число phi-точек пропорционально `sin(theta)` вдали от экватора. |
| `--filter` | `DEG` | Ограничить output backscatter cone. |
| `--point` | none | Legacy backscatter point mode. |

Приоритет:

| Величина | Приоритет |
|---|---|
| Theta | `--tgrid` > `--grid` > `--auto_tgrid` > `--auto` > default |
| Phi | `--nphi` > `--grid` > `--auto_phi` > `--auto` > default |

### Сетка по дифракционному пределу

Для PO-расчета с `--diffraction-grid DIV` флаг
`--diffraction-limit-grid FACTOR` строит сетку рассеяния из максимального
максимальной длины ребра граней текущей частицы `lmax` и длины волны
`lambda`. Это не диаметр частицы и не максимальное расстояние между
произвольными вершинами:

```text
xi = FACTOR * 0.69 * lambda / lmax
Ntheta = ceil(pi / xi)
Nphi,eq = round_up_6(max(12, ceil(2*pi / xi)))
```

Здесь `xi` измеряется в радианах, а результат содержит `Ntheta + 1` строк
theta от 0 до pi. `FACTOR=1` буквально использует оценку дифракционной
полосы; `FACTOR=0.5` уменьшает оба шага примерно вдвое, а значение больше
единицы делает сетку грубее. Это априорная оценка разрешения, а не
апостериорная гарантия точности: для контрольного расчета следует сравнить
как минимум `FACTOR=1` и `FACTOR=0.5` по всем нужным элементам Mueller.

С `--latitude-phi-grid` число азимутальных узлов зависит от строки:

```text
Nphi(theta) = min(Nphi,eq,
                  round_up_6(max(12, ceil(2*pi*sin(theta) / xi))))
```

Экваториальный шаг остается тем же, но около полюсов не вычисляются
избыточные азимуты, соответствующие почти одним и тем же направлениям.
Число узлов округляется до кратного шести; полюсная строка связывается с
соседней рабочей группой для устойчивой обработки базиса. Результат хранится
на полной прямоугольной выходной сетке `Nphi,eq`, поэтому формат файлов не
меняется.

Переменная сетка phi поддерживается для `--method po --diffraction-grid DIV` и
для `--so3-quaternion N`, по одному размеру частицы на процесс. Для SO(3)
предпочтителен `--scattering-diffraction-sampling Q`: общий
`--diffraction-sampling Q` сам задает регулярный ориентационный режим.
`--diffraction-limit-grid` конфликтует с явным
`--nphi`; `--latitude-phi-grid` пока не используется в общем serial
multi-size trace. Для каждого размера запускайте отдельный процесс, чтобы
`lmax` и сетка были пересчитаны корректно.

```bash
CUDA_VISIBLE_DEVICES=0,1,2,3 \
MBS_GPU_GROUPS=1 MBS_GPU_MULTI_MAX=4 \
gpu/bin/mbs_po_gpu_double --po -p 1 100 70 --ri 1.3116 0 -w 0.532 \
    -n 8 --diffraction-grid 2 --diffraction-limit-grid 0.5 \
    --latitude-phi-grid --threads 32 --close -o results/column
```

## GPU и multi-GPU

В CPU-расчете PO для фиксированной ориентации `--threads N`
распределяет азимутальные строки сетки направлений между
потоками OpenMP. Каждый поток записывает свои ячейки Jones/Mueller, а
вклады пучков в каждую ячейку суммируются в прежнем порядке. Поэтому
численный результат не зависит от числа потоков. Малые сетки остаются
последовательными, чтобы не тратить время на запуск OpenMP.

### За счет чего GPU быстрее

Самая дорогая часть PO - многократный расчет дифракционного edge-интеграла для большого числа направлений, пучков и ориентаций. CPU готовит пучки, а GPU параллелит дифракцию и, в быстрых режимах, сразу накопление Mueller.

| Операция | CPU path | GPU path | Гранулярность |
|---|---|---|---|
| Трассировка через грани | OpenMP по ориентациям/gamma blocks | `--gpu-trace-prefilter` переносит в CUDA упорядочивание граней и грубый отбор кандидатов; точные пересечения и разделение пучков остаются на CPU | Orientation/chunk |
| Packing пучков | CPU | CPU готовит `GpuBeam` buffers и копирует на device | Beam records |
| Edge diffraction | CPU loops | CUDA kernels | Beam x theta x phi x orientation |
| Когерентная Jones сумма | CPU complex arrays | Atomic или no-atomic device kernels | Grid cell/orientation |
| Jones -> Mueller | CPU после Jones sum | `mueller_batch_kernel` или fused kernel | Grid cell |
| Multi-`k_eq` | Повтор по размерам | Optional fused multi-`k_eq` CUDA pass | Size x beam x grid |
| FFT phi | Direct phi grid | cuFFT low-phi pass + zero-padding | Theta/orientation batches |

Масштаб работы:

```text
Nwork ~= Norient * Nbeams_per_orientation * Ntheta * Nphi
```

Эти вклады почти независимы до момента когерентной суммы Jones, поэтому хорошо ложатся на CUDA.

### Atomic и no-atomic режимы

| Режим | Kernels | Как копится результат | Когда используется |
|---|---|---|---|
| Atomic Jones | `diffraction_kernel` | Thread считает вклад beam/grid и делает `atomicAdd` в complex Jones buffer. Потом `mueller_batch_kernel` делает Mueller. | General fallback, включая часть diagnostic/no-shadow режимов. |
| No-atomic grid | `diffraction_grid_kernel` | Thread владеет output cell и сам проходит по пучкам, без конкуренции writers. | Full-only output, когда подходит layout. |
| Fused Mueller | `diffraction_grid_mueller_kernel`, `*_full_kernel`, `*_full8_kernel`, `*_mixed8_kernel` | Warp делит пучки одного направления между 32 потоками, сводит Jones через shuffle и сразу добавляет Mueller. | Быстрый production путь. |
| Staged orientation Mueller | `diffraction_grid_mueller_orient_kernel` + `reduce_mueller_orient_kernel` | Сначала per-orientation Mueller, затем редукция. | Большие batches, меньше contention. |
| Multi-`k_eq` fused | `diffraction_grid_mueller_multik_kernel` | Несколько близких `k_eq` из одного packed beam batch. | Shared-batch scans. |

В atomic path атомики идут по real/imag компонентам 2x2 Jones:

```text
J00.re, J00.im, J01.re, J01.im,
J10.re, J10.im, J11.re, J11.im
```

Для `_noshadow` тот же набор atomics пишется во второй Jones buffer. Это корректно для когерентной PO суммы, но может быть узким местом, если много пучков пишут в одну detector cell.

Управление:

| Переменная | Действие |
|---|---|
| `MBS_GPU_NO_ATOMICS=1` | Форсировать no-atomic путь, если поддерживается. |
| `MBS_GPU_NO_ATOMICS=0` | Форсировать atomic путь для сравнения/debug. |
| `MBS_GPU_FUSED_MUELLER=1` | Предпочитать fused diffraction-to-Mueller kernels. |
| `MBS_GPU_STAGE_MUELLER=1` | Включить staged per-orientation Mueller reduction. |
| `MBS_GPU_COMPACT_BEAMS=0` | Отключить компактную упаковку пучков с не более чем 8 вершинами и использовать общий формат для диагностики. |
| `MBS_GPU_WARP_BEAMS=0` | Отключить стандартное warp-разбиение пучков и вернуть последовательный обход одним CUDA-потоком для диагностики. |
| `MBS_GPU_WARP_GRID_3D=0` | Отключить стандартную специализацию warp-пути для трёхмерной сетки, сохранив общий warp-режим. |
| `MBS_GPU_TIMING=1` | Печатать breakdown count/pack/copy/kernels/d2h/add. |
| `MBS_GPU_BLOCK=N` | Переопределить размер блока CUDA: 64 (по умолчанию), 128 или 256. |
| `MBS_ORIENTATION_TIMING=1` | Печатать суммарное CPU-время поворота, трассировки и подготовки пучков. |

Для отдельного сравнения FP64-поворота вектора кватернионом и готовой
матрицей на GPU служит диагностическая цель:

```bash
make USE_CUDA=1 gpu_quaternion_probe
CUDA_VISIBLE_DEVICES=0 ./bin/gpu_quaternion_rotation_probe 65536 64 50
```

Пробник не выполняет расчет рассеяния и не заменяет профиль полной задачи.
Кватернионы задают равномерные ориентации на SO(3), но сами по себе не
ускоряют CPU-трассировку или CUDA-дифракцию.

Пакетная CUDA-трассировка автоматически включается для метода физической оптики
с CUDA-бэкендом. Явный `--gpu-trace-prefilter` включает тот же режим, а
`--no-gpu-trace-prefilter` отключает его для контрольного сравнения. В
проверенном автоматическом профиле используются полное топологическое
упорядочивание, кэш точной сортировки, консервативные флаги пропуска больших
списков, 64 потока большого ядра и отдельный неблокирующий поток CUDA для
каждого рабочего потока OpenMP. Точные пересечения многоугольников и разделение
пучков остаются на CPU.

Несколько OpenMP-потоков поддерживаются автоматически. Общий бюджет около
512 пучков делится между ними: при `--threads 8` один поток передаёт не более
64 пучков. Ограничение необходимо из-за скрытой памяти стека CUDA-ядер.
Нельзя назначать большой пакет только по показанию свободной видеопамяти:
восемь одновременных пакетов по 2048 пучков способны исчерпать 12 ГиБ. Если
начальный пакет всё же исчерпывает стек CUDA, он рекурсивно делится, а найденный
уменьшенный предел запоминается для всех следующих пакетов и ориентаций.

При включённом `--trace-limit-retries` рабочие копии также совместно запоминают
наибольший успешно проверенный предел числа пучков. Следующие ориентации не
повторяют заведомо недостаточный первый проход. Число использованных уровней
по-прежнему отсчитывается от исходного `--trace-max-beams`, поэтому верхняя
граница повторов не увеличивается. На тесте из 32 ориентаций это уменьшило
медианное время с 11,12 до 8,00 с и число полных повторов с 30 до 7. Отдельные
неблокирующие потоки уменьшили время 256 ориентаций с 61,39 до 42,26 с. В обоих
A/B-тестах итоговые таблицы и интегральные характеристики совпали побитово.

Независимые процессы на одной физической видеокарте последовательно выполняют
только тяжёлую часть пакета трассировки через блокировку, привязанную к PCI-
адресу карты. Рабочие потоки OpenMP внутри одного процесса по-прежнему
перекрываются в отдельных потоках CUDA. Этот автоматический режим не даёт
суммарному резерву стека CUDA-ядер исчерпать память карты и не блокирует ядра
дифракции.

Экспертные переопределения требуют `--allow-experimental-environment`:

| Переменная | По умолчанию | Действие |
|---|---:|---|
| `MBS_GPU_TRACE_OPENMP=0/1` | `1` | Отключить или включить CUDA-трассировку из нескольких OpenMP-потоков. |
| `MBS_GPU_TRACE_BATCH_BEAMS=N` | начальный `ceil(512/threads)`, минимум 64 | Размер пакета одного потока; после OOM уменьшается автоматически, переопределять только после измерений. |
| `MBS_GPU_TRACE_FULL_SORT=0/1` | `1` | Полное топологическое упорядочивание граней на CUDA. |
| `MBS_GPU_TRACE_SORT_CACHE=0/1` | `1` | Повторное использование точной сортировки для одинаковых ключей. |
| `MBS_GPU_TRACE_LARGE_SKIP_FLAGS=0/1` | `1` | Консервативные флаги пропуска для списков из 25--256 граней. |
| `MBS_GPU_TRACE_LARGE_THREADS=N` | `64` | Размер блока большого ядра: 64, 128 или 256. |
| `MBS_GPU_TRACE_CACHE_FACETS=0/1` | `1` | Сохранять упакованную геометрию граней в CUDA workspace потока. |
| `MBS_GPU_TRACE_MIN_CANDIDATES=N` | `1024` | Минимальное общее число кандидатов для отправки пакета на GPU. |
| `MBS_GPU_TRACE_NONBLOCKING_STREAM=0/1` | `1` | Отдельный неблокирующий поток CUDA; `0` оставлен для диагностики последовательного режима. |
| `MBS_GPU_TRACE_PROCESS_LOCK=0/1` | `1` | Последовательно выполнять тяжёлые пакеты трассировки независимых процессов на одной физической GPU; `0` только для контролируемой диагностики. |
| `MBS_GPU_TRACE_PREFILTER_FIRST=0/1` | `0` | Диагностический грубый отбор перед сортировкой. |
| `MBS_GPU_TRACE_MARGIN=X` | автоматически | Переопределить консервативный запас проекционных границ. |
| `MBS_GPU_TRACE_TIMING=1` | `0` | Печатать времена H2D, малых и больших ядер и D2H. |
| `MBS_GPU_TRACE_VERIFY_SORT=1` | `0` | Сравнивать CUDA-сортировку с точным CPU-результатом. |

### Multi-size и multi-GPU

Переменная `MBS_GPU_GROUPS=1` включает отдельный точный планировщик для
переменной latitude-phi сетки. Независимые группы строк theta распределяются
динамически по видимым GPU. Трассировка ориентаций и подготовленные пучки
остаются общими в host RAM; на каждой карте создается один CUDA workspace.
Для каждого orientation chunk упакованные пучки, веса и смещения копируются
на данную карту один раз и повторно используются следующими theta-группами.

Число карт задается `CUDA_VISIBLE_DEVICES`; `MBS_GPU_MULTI=N` или
`MBS_GPU_MULTI_MAX=N` ограничивает число workers. Значение
`MBS_GPU_MULTI=0` оставляет одну карту. Этот путь требует прямой CUDA-
дифракции: при FFT-интерполяции или theta-zone beam filtering программа
сохраняет прежний путь. Ускорение не обязано быть линейным, поскольку
CPU-трассировка, упаковка, PCIe-копирование, редукция и неодинаковая стоимость
theta-групп остаются в полном времени.

| Флаг | Аргументы | Описание |
|---|---:|---|
| `--multigrid` | `Dmin Dmax N` | Log-spaced scan по `Dmax`. |
| `--multikeq` | `Kmin Kmax N` | Log-spaced scan по `k_eq`. |
| `--multikeq_list` | `FILE` | Точный список `k_eq`, один на строку. |
| `--multigrid_parallel` | `N` | Запуск scan points как child processes; `0` значит auto по списку GPU или числу CPU cores. |
| `--multigrid_threads` | `N` | OpenMP threads на child process. |
| `--gpu_devices` | `LIST` | Список CUDA devices, например `0,1,2,3,4`. |
| `--multikeq_shared_batches` | none | Batch близких `k_eq` на GPU child с reuse tracing от максимального `k_eq`. |
| `--multikeq_batch_ratio` | `R` | Максимальное `kmax/kmin` внутри shared batch; default `1.05`. |

Parallel scheduler работает через отдельные процессы. Parent убирает scan-флаг
из команды, создает подпапку в `-o` для каждого размера или batch, затем
запускает child с `--rs SIZE`, `--k_eq K` или автоматически созданным
`--multikeq_list FILE`. Лог каждого child пишется в
`<output>/<label>.run.log`.

С `--gpu` каждому child назначается одна CUDA-карта через
`CUDA_VISIBLE_DEVICES` и `MBS_GPU_SLOT`. Если `--gpu_devices` не задан, берутся
CUDA devices, видимые процессу. `--multigrid_parallel 0` обычно означает: для
GPU использовать число видимых/перечисленных карт, для CPU ограничиться
доступными threads. На общей машине лучше задавать `--multigrid_parallel N`
явно.

Exact multi-GPU: один процесс на активный GPU slot, каждый размер трассируется
независимо:

```bash
gpu/bin/mbs_po_gpu_float_fast --po --pf particle.dat \
    --multikeq_list keq.txt --ri 1.6 0.002 -w 1.064 -n 12 \
    --oldauto 2 --grid 0 180 600 180 --fft --close \
    --multigrid_parallel 0 --multigrid_threads 16 --gpu_devices 0,1,2,3,4 \
    -o scan_exact
```

Shared-batch `k_eq`:

```bash
gpu/bin/mbs_po_gpu_float_fast --po --pf particle.dat \
    --multikeq_list keq.txt --ri 1.6 0.002 -w 1.064 -n 12 \
    --oldauto 2 --grid 0 180 600 180 --fft --close \
    --multigrid_parallel 0 --multigrid_threads 16 --gpu_devices 0,1,2,3,4 \
    --multikeq_shared_batches --multikeq_batch_ratio 1.05 \
    -o scan_shared
```

`--multikeq_shared_batches` применим только к scan по `k_eq`. Близкие значения
группируются по `--multikeq_batch_ratio`; child трассирует максимальный `k_eq`
в batch и переиспользует подготовленные пучки для меньших `k_eq`. Это быстрее,
но не полностью эквивалентно независимым oldauto-сеткам. Для проверки точности
уменьшай ratio, а если каждому размеру нужна своя oldauto reference grid,
используй exact mode без shared batches.

Fused multi-`k_eq`:

```bash
MBS_GPU_MULTI_K_FULL=1 gpu/bin/mbs_po_gpu_float_fast ...
```

Для больших direct-сеток `theta x phi` главным ограничением может быть host RAM,
даже если GPU memory свободна. В oldauto log появляется строка:

```text
Oldauto/random memory: gamma chunk=... (GPU host-RAM guard, MemAvailable=... MB, VmRSS=... MB, grid-transient=... MB, budget=... MB)
```

Если `grid-transient` сравним с budget, один только меньший `--chunk` проблему
не решит: нужно уменьшать `Nphi/Nth`, включать `--fft`, увеличивать host-memory
budget или делить задачу на меньшие независимые сетки.

Для oldauto/random checkpoint сетка Mueller без тени сохраняется только если
запрошен no-shadow output. В обычных full-output расчетах код не держит и не
пишет второй полный массив `theta x phi x 4 x 4`, что заметно снижает RSS на
больших direct-сетках.

Основные настройки scheduler и memory guard:

| Переменная | Смысл |
|---|---|
| `MBS_PARALLEL_MEM_FRACTION` | Доля текущей `MemAvailable`, делится между GPU child processes; default `0.70`. |
| `MBS_HOST_MEM_BUDGET_MB` | Жесткий host-memory budget для child. Scheduler задает автоматически, если переменная уже не выставлена. |
| `MBS_HOST_MEM_FRACTION` | Доля текущей доступной ОЗУ для блоков подготовленных ориентаций. |
| `MBS_HOST_MEM_RESERVE_MB` | Резерв host RAM, который oldauto/random старается оставить свободным; default `4096`. |
| `MBS_OLDAUTO_GRID_MEM_SAFETY` | Запас для transient массивов direct-grid; default `1.25`. |
| `MBS_OLDAUTO_BYTES_PER_GAMMA_MB` | Оценка памяти на один gamma для prepared beams. |
| `MBS_OLDAUTO_GAMMA_CHUNK` | Default oldauto gamma chunk до автоматического memory clamp. |
| `MBS_PARALLEL_SHARED_KEQ=1` | Environment equivalent для `--multikeq_shared_batches`. |
| `MBS_PARALLEL_KEQ_BATCH_RATIO=R` | Environment default для `--multikeq_batch_ratio`. |

## Все флаги

### Метод, частица, физические параметры

| Флаг | Аргументы | Описание |
|---|---:|---|
| `--po` | none | Physical Optics. |
| `--go` | none | Geometrical Optics. |
| `-p` | `TYPE L D [extra]` | Built-in particle. |
| `--pf` | `FILE` | Particle from file. |
| `--rs` | `SIZE` | Resize file particle by Dmax. |
| `--k_eq` | `X` | Resize by equivalent size parameter. |
| `--ri` | `Re Im` | Refractive index. |
| `-w` | `LAMBDA` | Wavelength, micrometers. |
| `-n` | `N` | Max internal reflections/refractions. |
| `--abs` | none | Enable absorption accounting. |
| `--abs_points` | `N` или `all` | Absorption samples. |

### Ориентации

| Флаг | Аргументы | Описание |
|---|---:|---|
| `--fixed` | `BETA GAMMA` | Одна ориентация в градусах. |
| `--random` | `Nb Ng` | Regular beta/gamma grid. |
| `--sobol` | `N` | Sobol orientations. |
| `--so3_quat` | `N` | SO(3) с сокращением по симметрии; alpha интегрируется по phi. |
| `--so3_mirror_audit` | none | Явные вложенные отраженные пары для проверки `--mirror_gamma`; только четное `N`. |
| `--so3_full_quat` | `N` | Прямая кватернионная проверка полной SO(3). |
| `--sobol_seed` | `N S` | Sobol with Owen seed. |
| `--sobol_ring` | `Nb Ng` | Sobol beta x gamma rings. |
| `--hammersley` | `N` | Hammersley set. |
| `--lattice` | `N` | Rank-1 lattice. |
| `--lattice_z` | `N Z` | Rank-1 lattice with generator. |
| `--euler_quad` | `Nb Ng` | Euler/Gauss quadrature. |
| `--euler_adapt` | `Nb NgMax` | Adaptive gamma count. |
| `--montecarlo` | `N` | Monte Carlo orientations. |
| `--adaptive` | `EPS` | Adaptive orientation convergence. |
| `--auto` | `EPS` | Auto theta, phi, orientation count. |
| `--autofull` | `EPS` | Auto `n`, theta, phi, orientation count. |
| `--oldautofull` | `EPS` | Autofull + oldauto final grid. |
| `--owen_avg` | `K` | Average final Owen seeds. |
| `--owen_seeds` | `S...` | Explicit Owen seeds. |
| `--oldauto` | `DIV` | Physics-based regular grid. |
| `--ring_points` | `N` | Points per diffraction ring estimate. |
| `--mirror_gamma` | none | Mirror half gamma domain. |
| `--orientfile` | `FILE` | Прочитать `beta gamma` в градусах из файла. |
| `--b` | `B1 B2` | Beta range for `--random`. |
| `--g` | `G1 G2` | Gamma range for `--random`. |
| `--maxorient` | `N` | Max adaptive orientations. |
| `--chunk` | `N` | Orientation/gamma chunk size. |
| `--coh_orient` | none | Отключенный устаревший режим; не использовать. |
| `--pole` | none | Exact beta-pole gamma shortcut. |
| `--sym` | `Sb Sg` | Override symmetry. |

### Сетка, ускорение, cutoffs

| Флаг | Аргументы | Описание |
|---|---:|---|
| `--grid` | `T1 T2 Nphi Nth` | Uniform theta range. |
| `--grid` | `R Nphi Nth` | Backscatter cone. |
| `--tgrid` | `FILE` | Non-uniform theta grid. |
| `--auto_tgrid` | `EPS` | Adaptive theta grid. |
| `--auto_phi` | none | Auto phi count. |
| `--nphi` | `N` | Override phi count. |
| `--diffraction-limit-grid` | `FACTOR` | Шаги theta и экваториального phi: `FACTOR*0.69*lambda/lmax`. |
| `--latitude-phi-grid` | none | Число phi-точек по строкам пропорционально `sin(theta)`. |
| `--threads` | `N` | OpenMP worker threads. |
| `--gpu` | none | CUDA backend. |
| `--cpu` | none | Force CPU backend in GPU build. |
| `--fft` | none | cuFFT phi interpolation. |
| `--beam_cutoff` | `EPS` | Set beam `J` and area cutoffs. |
| `--beam_cutoff_j` | `EPS` | Skip by relative `|J|^2`. |
| `--beam_cutoff_area` | `EPS` | Skip by relative area. |
| `--beam_cutoff_importance` | `EPS` | Skip by relative `|J|^2*area`. |
| `--trace_cutoff` | `EPS` | Set trace pruning cutoffs. |
| `--trace_cutoff_j` | `EPS` | Trace prune by `|J|^2`. |
| `--trace_cutoff_area` | `EPS` | Trace prune by area. |
| `--trace_cutoff_importance` | `EPS` | Trace prune by `|J|^2*area`. |
| `--trace_max_beams` | `N` | Аварийно завершить расчёт после `N` узлов дерева в одной ориентации; `0` отключает предел. |
| `--gpu-trace-prefilter` | none | Явно включить автоматический пакетный профиль CUDA-трассировки. |
| `--no-gpu-trace-prefilter` | none | Отключить автоматическое CUDA-упорядочивание и отбор кандидатов. |
| `--trace_prefilter` | none | Enable CPU projected-AABB prefilter. |
| `--no_trace_prefilter` | none | Disable CPU prefilter. |
| `--trace_prefilter_margin` | `M` | AABB prefilter margin. |
| `--trace_prefilter_stats` | none | Print prefilter counters. |
| `-r` | `RATIO` | Beam area restriction ratio. |

### Multi-size и multi-GPU

| Флаг | Аргументы | Описание |
|---|---:|---|
| `--multigrid` | `Dmin Dmax N` | Log-spaced scan по максимальному размеру частицы. |
| `--multikeq` | `Kmin Kmax N` | Log-spaced scan по equivalent size parameter `k_eq`. |
| `--multikeq_list` | `FILE` | Точные значения `k_eq`, один на строку. |
| `--multigrid_parallel` | `N` | Запускать scan points как child processes; `0` выбирает число jobs автоматически. |
| `--multigrid_threads` | `N` | OpenMP threads на один child process. |
| `--gpu_devices` | `LIST` | CUDA devices, назначаются child processes по кругу. |
| `--multikeq_shared_batches` | none | Группировать близкие `k_eq` и переиспользовать трассировку от максимального `k_eq` в batch. |
| `--multikeq_batch_ratio` | `R` | Максимальное `kmax/kmin` внутри shared batch; default `1.05`. |

### Output, diagnostics, legacy

| Флаг | Аргументы | Описание |
|---|---:|---|
| `-o` | `NAME` | Output path/name. |
| `--close` | none | Exit after computation. |
| `--log` | `SEC` | Progress log interval. |
| `--checkpoint` | none | Save/resume long runs. |
| `--save_betas` | none | Save per-beta Mueller files. |
| `--full_only` | none | Write only full Mueller output; default. |
| `--noshadow_output` | none | Also write `_noshadow`. |
| `--shadow` | none | Legacy flag, currently no effect. |
| `--shadow_off` | none | Disable shadow beam. |
| `--incoh` | none | Incoherent per-beam Mueller sum. |
| `--jones` | none | Write Jones matrices where supported. |
| `--karczewski` | none | Use Karczewski polarization matrix. |
| `--legacy_sign` | none | Old Fresnel sign convention. |
| `--ot_phase_avg` | none | Average OT extinction over phase period. |
| `--ot_phase_shift` | `F` | Diagnostic OT phase shift. |
| `--ot_ping` | `D` | Legacy OT phase rotation. |
| `--filter` | `DEG` | Restrict output to backscatter cone. |
| `--point` | none | Legacy backscatter point mode. |
| `--tr` | `FILE` | Load trajectory file. |
| `--all` | none | Calculate all loaded trajectories. |
| `--gr` | none | Output trajectory groups. |
| `--forced_convex` | none | Force convex processing. |
| `--forced_nonconvex` | none | Force nonconvex processing. |
| `--help`, `-h` | none | Short help. |
| `--help-debug` | none | Full debug/experimental help. |

## Переменные окружения

Все активные переменные `MBS_*` записываются в итоговый журнал. Переменные,
которые меняют выборку, геометрию, отсечения или численный результат, по
умолчанию запрещены. Флаг `--allow-experimental-environment` явно разрешает их
и сохраняет имена и значения в журнале; в рабочем расчёте их предпочтительно
снять. Для блоков ориентаций программа начинает с 16 пробных ориентаций,
измеряет фактическую память вложенных пучков и траекторий и подбирает следующий
блок по текущей доступной ОЗУ и ограничениям ниже.

| Переменная | Смысл |
|---|---|
| `CUDA_VISIBLE_DEVICES` | Стандартный выбор видимых CUDA devices. |
| `OMP_NUM_THREADS` | OpenMP threads по умолчанию. |
| `MBS_GPU_ALLOW_FALLBACK=1` | Разрешить старый fallback GPU->CPU после CUDA ошибки. Только debug. |
| `MBS_GPU_MEM_FRACTION` | Доля свободной GPU памяти для buffers. |
| `MBS_GPU_NO_ATOMICS` | Выбор no-atomic (`1`) или atomic (`0`) path. |
| `MBS_GPU_FUSED_MUELLER` | Fused diffraction-to-Mueller kernels. |
| `MBS_GPU_STAGE_MUELLER` | Staged per-orientation Mueller reduction. |
| `MBS_GPU_NO_VERTEX_CACHE` | Отключить cached/packed vertex path. |
| `MBS_GPU_COMPACT_BEAMS=0` | Отключить компактную упаковку пучков с не более чем 8 вершинами. |
| `MBS_GPU_WARP_BEAMS=0` | Вернуть последовательный CUDA-обход пучков одним потоком для диагностики. |
| `MBS_GPU_WARP_GRID_3D=0` | Отключить только специализацию warp-пути для трёхмерной сетки. |
| `MBS_GPU_TIMING` | Печатать CUDA timing breakdown. |
| `MBS_GPU_BLOCK` | Размер блока CUDA: 64 (по умолчанию), 128 или 256. |
| `MBS_GPU_TRACE_OPENMP=0/1` | Автоматическая CUDA-трассировка с несколькими OpenMP-потоками; по умолчанию `1`. |
| `MBS_GPU_TRACE_BATCH_BEAMS=N` | Переопределить пакет `ceil(512/threads)` на поток; большие значения могут исчерпать стек CUDA. |
| `MBS_GPU_TRACE_FULL_SORT=0/1` | Точное топологическое упорядочивание CUDA; по умолчанию `1`. |
| `MBS_GPU_TRACE_SORT_CACHE=0/1` | Кэш точной сортировки; по умолчанию `1`. |
| `MBS_GPU_TRACE_LARGE_SKIP_FLAGS=0/1` | Консервативные флаги пропуска больших списков; по умолчанию `1`. |
| `MBS_GPU_TRACE_LARGE_THREADS=N` | Размер блока большого ядра; по умолчанию `64`, допустимо 64/128/256. |
| `MBS_GPU_TRACE_NONBLOCKING_STREAM=0/1` | Неблокирующие потоки workers; по умолчанию `1`; `0` включает последовательный диагностический режим. |
| `MBS_GPU_TRACE_PROCESS_LOCK=0/1` | Межпроцессная блокировка тяжёлых пакетов трассировки на физической GPU; по умолчанию `1`; `0` только для диагностики. |
| `MBS_GPU_TRACE_TIMING=1` | Печатать времена передачи и ядер CUDA-трассировки. |
| `MBS_GPU_TRACE_VERIFY_SORT=1` | Сверять CUDA-упорядочивание с CPU. |
| `MBS_ORIENTATION_TIMING` | Печатать разбиение CPU-времени поворота, трассировки и подготовки пучков. |
| `MBS_HOST_MEM_FRACTION` | Доля текущей доступной ОЗУ для измеряемых блоков подготовленных ориентаций. |
| `MBS_HOST_MEM_RESERVE_MB` | Резерв host RAM в MB. |
| `MBS_HOST_MEM_BUDGET_MB` | Жесткий host RAM budget. |
| `MBS_OLDAUTO_GRID_MEM_SAFETY` | Коэффициент запаса для оценки host RAM больших сеток `theta x phi` в oldauto; default `1.25`. |
| `MBS_OLDAUTO_BYTES_PER_GAMMA_MB` | Оценка памяти prepared-beams на один gamma для oldauto/random memory guard. |
| `MBS_OLDAUTO_GAMMA_CHUNK` | Default oldauto gamma chunk до автоматического memory clamp. |
| `MBS_OLDAUTO_BETA_MIDPOINT=0/1` | Midpoint beta rings в oldauto/random quadrature. |
| `MBS_OLDAUTO_GAMMA_STAGGER=1` | Сдвигать gamma samples между beta rings, чтобы уменьшить aliasing узких событий. |
| `MBS_FFT_PHI_FACTOR` | Override reduced phi factor для FFT. |
| `MBS_FFT_THETA_FACTOR` | Override theta batching factor для FFT. |
| `MBS_FFT_CHECK=1` | FFT diagnostic checks. |
| `MBS_FFT_ADAPTIVE_PHI=1` | Adaptive reduced-phi behavior в FFT backend. |
| `MBS_GPU_MULTI=0` | Отключить automatic multi-orientation GPU batching. |
| `MBS_GPU_MULTI_MAX=N` | Ограничить GPU multi batching. |
| `MBS_GPU_GROUPS=1` | Распределить группы theta переменной latitude-phi сетки по видимым GPU. |
| `MBS_GPU_MULTI_K_FULL=1` | Experimental fused multi-`k_eq` diffraction. |
| `MBS_GPU_BEAM_STATS=1` | Печатать GPU beam packing/count diagnostics. |
| `MBS_SHARED_BETA_GROUP=N` | Override beta grouping для shared batch. |
| `MBS_SHARED_ORIENT_CHUNK=N` | Override orientation chunk для shared batch. |
| `MBS_PARALLEL_MEM_FRACTION` | Memory fraction scheduler для parallel multi-size. |
| `MBS_PARALLEL_SHARED_KEQ=1` | Environment equivalent для shared `k_eq` batching. |
| `MBS_PARALLEL_KEQ_BATCH_RATIO=R` | Default batch ratio для shared `k_eq`. |

## Выходные файлы

Основной `.dat`:

```text
ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44
```

| Колонка | Смысл |
|---|---|
| `ScAngle` | Угол рассеяния theta в градусах. |
| `2pi*dcos` | Вес кольца телесного угла для theta row. |
| `Mij` | Элементы Mueller после суммирования пучков и усреднения ориентаций. |

Дополнительные outputs:

| Output | Как включается | Смысл |
|---|---|---|
| `_noshadow` | `--noshadow_output` | Mueller без shadow/external beam. |
| `_betas/` | `--save_betas` | Per-beta diagnostics. |
| checkpoint files | `--checkpoint` | Resume для длинных orientation-grid runs. |
| Jones output | `--jones` | Raw Jones matrices в поддержанных PO modes. |

Предупреждения про optical-theorem/integral mismatch являются диагностикой. Для non-absorbing runs физическое поглощение фиксируется нулем; integral characteristic не заменяет выходную Mueller матрицу.
