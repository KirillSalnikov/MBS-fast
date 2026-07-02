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
make -C gpu -j          # GPU float_fast по умолчанию
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
# GPU binary по умолчанию: float + fast math
PATH=/usr/local/cuda/bin:$PATH make -C gpu -j

# Явные варианты
make -C gpu float       -j   # gpu/bin/mbs_po_gpu_float
make -C gpu float_fast  -j   # gpu/bin/mbs_po_gpu_float_fast
make -C gpu double_fast -j   # gpu/bin/mbs_po_gpu_double_fast
```

| Target | Binary | CUDA точность | Fast math | Когда использовать |
|---|---|---:|---:|---|
| `make -C gpu` | `gpu/bin/mbs_po_gpu_float_fast` | FP32 | да | Быстрые production сканы |
| `make -C gpu float` | `gpu/bin/mbs_po_gpu_float` | FP32 | нет | Проверка FP32 без fast math |
| `make -C gpu double_fast` | `gpu/bin/mbs_po_gpu_double_fast` | FP64 | да | Контрольные GPU расчеты |
| `make -C gpu double_debug` | `gpu/bin/mbs_po_gpu_double_debug` | FP64 | да | `--help-debug` и диагностика |

В split GPU build CUDA backend включен по умолчанию. Флаг `--gpu` можно писать, но он не обязателен. Флаг `--cpu` принудительно запускает CPU backend внутри GPU-capable бинарника.

Архитектура GPU определяется через `nvidia-smi`. При необходимости задавайте вручную:

```bash
make -C gpu double_fast -j GPU_ARCH=86   # Ampere
make -C gpu float_fast  -j GPU_ARCH=89   # Ada
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
| `--orientfile` | `FILE` | Пользовательские ориентации | Одна пара beta/gamma на строку. |
| `--sobol` | `N` | Quasi-random average | Хорош для сходимости сканов. |
| `--sobol_seed` | `N S` | Sobol/Owen с seed | Повторяемые проверки сходимости. |
| `--sobol_ring` | `Nb Ng` | Sobol beta + uniform gamma rings | Гибридная сетка. |
| `--so3_quat` | `N` | Full SO(3) quaternions | Не опирается на beta/gamma symmetry domain. |
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
| `--mirror_gamma` | none | Половина gamma domain + зеркалирование Mueller. |
| `--sym` | `Sb Sg` | Override symmetry: beta range `pi/Sb`, gamma range `2*pi/Sg`. |
| `--b` | `B1 B2` | Диапазон beta для `--random`, градусы. |
| `--g` | `G1 G2` | Диапазон gamma для `--random`, градусы. |
| `--maxorient` | `N` | Верхняя граница adaptive orientations. |
| `--chunk` | `N` | Chunk size по ориентациям/gamma. |
| `--coh_orient` | none | Legacy coherent-across-orientations. |
| `--pole` | none | Одна gamma на точных beta poles. |
| `--owen_avg` | `K` | Усреднение `K` Owen seeds в `--autofull`. |
| `--owen_seeds` | `S...` | Явный список Owen seeds. |

## Сетка рассеяния

| Флаг | Аргументы | Описание |
|---|---:|---|
| `--grid` | `T1 T2 Nphi Nth` | Uniform theta range от `T1` до `T2`; выводит `Nth + 1` theta строк. |
| `--grid` | `R Nphi Nth` | Backscatter cone радиуса `R` градусов. |
| `--tgrid` | `FILE` | Неравномерная theta сетка, градусы, одна строка на theta. |
| `--auto_tgrid` | `EPS` | Adaptive theta grid через bisection. |
| `--auto_phi` | none | Автоматический выбор `Nphi`. |
| `--nphi` | `N` | Override phi count, максимальный приоритет. |
| `--filter` | `DEG` | Ограничить output backscatter cone. |
| `--point` | none | Legacy backscatter point mode. |

Приоритет:

| Величина | Приоритет |
|---|---|
| Theta | `--tgrid` > `--grid` > `--auto_tgrid` > `--auto` > default |
| Phi | `--nphi` > `--grid` > `--auto_phi` > `--auto` > default |

## GPU и multi-GPU

### За счет чего GPU быстрее

Самая дорогая часть PO - многократный расчет дифракционного edge-интеграла для большого числа направлений, пучков и ориентаций. CPU готовит пучки, а GPU параллелит дифракцию и, в быстрых режимах, сразу накопление Mueller.

| Операция | CPU path | GPU path | Гранулярность |
|---|---|---|---|
| Трассировка через грани | OpenMP по ориентациям/gamma blocks | В основном CPU; `--gpu_trace` только prefilter кандидатов | Orientation/chunk |
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
| Fused Mueller | `diffraction_grid_mueller_kernel`, `*_full_kernel`, `*_full8_kernel`, `*_mixed8_kernel` | Kernel считает Jones локально и сразу переводит/добавляет Mueller. | Быстрый production путь. |
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
| `MBS_GPU_TIMING=1` | Печатать breakdown count/pack/copy/kernels/d2h/add. |
| `MBS_GPU_BLOCK=N` | Override CUDA block size. |

### Multi-size и multi-GPU

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
| `MBS_HOST_MEM_FRACTION` | Доля total host RAM для oldauto/random chunking. |
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
| `--so3_quat` | `N` | Full SO(3) quaternions. |
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
| `--orientfile` | `FILE` | Load beta/gamma from file. |
| `--b` | `B1 B2` | Beta range for `--random`. |
| `--g` | `G1 G2` | Gamma range for `--random`. |
| `--maxorient` | `N` | Max adaptive orientations. |
| `--chunk` | `N` | Orientation/gamma chunk size. |
| `--coh_orient` | none | Legacy coherent orientation mode. |
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
| `--trace_max_beams` | `N` | Abort orientation after `N` beam nodes; `0` disables. |
| `--gpu_trace` | none | Experimental CUDA candidate prefilter. |
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
| `MBS_GPU_TIMING` | Печатать CUDA timing breakdown. |
| `MBS_GPU_BLOCK` | Override CUDA block size. |
| `MBS_HOST_MEM_FRACTION` | Доля host RAM для oldauto/random chunking. |
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
