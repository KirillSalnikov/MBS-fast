# Fullauto: краткая памятка

[Подробное руководство с математикой](FULLAUTO_ADAPTIVE_RU.html).
Текущая инструкция к подготовке 2.0; EPS=0.02 означает 2%.

## CPU: минимум параметров

```bash
make -C cpu -j2
cpu/bin/mbs_po_mpi --method po --backend cpu \
  --particle 1 5 2 --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --output results/column_5um
```

Тип 1 — шестигранный столбик, 5 — осевая длина, 2 — диаметр основания.
Все длины в тех же единицах, что длина волны.
Для файла: `--particle-file examples/cube.particle --resize-dmax-um 5`.

## GPU: конус назад

```bash
make -C gpu fp64 GPU_ARCH=70 -j2
CUDA_VISIBLE_DEVICES=0,1 gpu/bin/mbs_po_gpu_double \
  --method po --backend cuda --particle 1 100 50 \
  --ri 1.3116 0 --wavelength-um 0.532 \
  --fullauto 0.02 --fullauto-theta-range 170 180 \
  --fullauto-gpus 0,1 --output results/column_back
```

Вперёд замените диапазон на `--fullauto-theta-range 0 25` и выберите другую
папку результата. GPU_ARCH=70 для V100, 86 для RTX 3080 Ti, 89 для RTX 4090.

## Что можно закрепить

| Флаг | Эффект |
|---|---|
| `--max-reflections 8` | Фиксировать глубину, выключить её автопоиск. |
| `--phi-points 4` | Фиксировать четыре азимута. |
| `--fullauto-initial 512` | Старт 512 ориентаций на seed; дальнейший рост остаётся. |
| `--threads 4` | Фиксировать число потоков. |
| `--fullauto-analytic-controls off` | Отключить аналитические контроли. |
| `--cutoff-profile off` | Без отсечения. |
| `--fullauto-rotation-symmetry off --fullauto-mirror off` | Не сокращать область по симметриям. |
| `--max-orientations 1048576` | Ограничить N; предел не является сходимостью. |

По умолчанию: восемь Sobol/Owen seed, подбор глубины 6…18, автоматические φ,
проверяемые аналитические контроли, coherent PO balanced с парной поправкой off,
интерполяционный допуск EPS/4. Для 2% это 0,5%. Численная точность внутри PO
не равна физической точности относительно решения Максвелла.

## Смотреть и восстанавливать

```bash
python3 -m json.tool results/column_5um/fullauto_status.json
python3 -m json.tool results/column_5um/intermediate_current.json
python3 -m json.tool results/column_5um/execution_progress.json
```

`converged` — включённые проверки пройдены. Точность относится к последнему
опубликованному раунду; активность GPU сама её не обновляет.
`mueller_fullauto.dat` — полная сетка; `mueller_evaluated.dat` — прямые углы;
`mueller_fullauto_annotated.csv` — происхождение строк и ошибки.

После отключения повторите **исходную команду тем же бинарником** в ту же папку.
Сохраняйте папку целиком. Завершённые блоки используются повторно. Новая версия
требует подтверждённой совместимости, а не ручной замены хэша.

Полная команда с закреплёнными параметрами, определение каждого флага,
усреднение Sobol, моменты Френеля, остаток интерполяции и GO разобраны
в [подробном руководстве](FULLAUTO_ADAPTIVE_RU.html).
