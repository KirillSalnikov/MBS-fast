#!/usr/bin/env python3
"""Build offline native-fullauto documentation using completed audit evidence."""
import argparse
import base64
import html
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def figure(fig, caption):
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', dpi=155, bbox_inches='tight')
    plt.close(fig)
    return '<figure><img alt="'+html.escape(caption)+'" src="data:image/png;base64,'+base64.b64encode(buffer.getvalue()).decode()+'"><figcaption>'+caption+'</figcaption></figure>'


def body(validation):
    cases = []
    plots = ''
    for name in ['cpu_final_one_percent', 'cpu_backscatter_depth12']:
        root = validation/name
        if not (root/'independent_fullauto_audit.json').exists():
            continue
        audit = json.loads((root/'independent_fullauto_audit.json').read_text())
        rows = np.loadtxt(root/'python_vs_native.csv', delimiter=',', skiprows=1)
        state = json.loads((root/'fullauto_status.json').read_text())
        cases.append((name, audit, state))
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.7))
        axes[0].plot(rows[:, 0], rows[:, 1], label='native fullauto')
        axes[0].plot(rows[:, 0], rows[:, 2], '--', label='direct full theta grid')
        axes[0].set_ylabel('M11')
        axes[1].plot(rows[:, 0], -rows[:, 3]/rows[:, 1], label='native fullauto')
        axes[1].plot(rows[:, 0], -rows[:, 4]/rows[:, 2], '--', label='direct full theta grid')
        axes[1].set_ylabel('-M12 / M11')
        axes[2].plot(rows[:, 0], rows[:, 5]*100, label='max interpolation residual, 16 entries')
        axes[2].plot(rows[:, 0], rows[:, 6]*100, label='max pointwise 95% interval')
        axes[2].axhline(audit['epsilon']*25, color='black', linestyle=':', label='interpolation budget')
        axes[2].axhline(audit['epsilon']*75, color='grey', linestyle=':', label='statistical budget')
        axes[2].set_ylabel('% of M11')
        for ax in axes:
            ax.set_xlabel('theta, deg')
            ax.grid(alpha=.2)
            ax.legend(fontsize=7)
        plots += figure(fig, f'{name}: реальные матрицы M11 и M12 и независимая проверка всех 16 элементов. Пунктир — отдельный прямой расчёт всех запрошенных углов при тех же seeds и бюджетах.')
    assert cases, 'No completed independent audits'
    records = ''
    for name, audit, state in cases:
        dense = audit['independent_dense']
        records += '<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in (
            name, state['target_depth'], audit['phi'], str(audit['final_counts_per_seed']),
            f"{audit['requested_theta']} → {audit['evaluated_theta']} / {audit['support_theta']}",
            f"{100*audit['statistical95_scaled_M11']:.6f}%",
            f"{100*audit['guard_residual95_scaled_M11']:.6f}%",
            f"{100*dense['matrix_error_scaled_M11']:.6f}%",
            f"{100*dense['residual_plus95_scaled_M11']:.6f}%"))+'</tr>'
    max_difference = max(a['matrix_difference_scaled_M11'] for _, a, _ in cases)
    return '''
<h2>Нативный режим <code>-fullauto EPS</code></h2>
<p><code>EPS</code> — относительный допуск: <code>0.01</code> означает 1%. Запуск C++ сам выполняет пилоты, выбирает число азимутов φ, фиксирует контрольные коэффициенты, распределяет ориентации между уровнями отражений и уточняет θ. Для расчёта Python не нужен. NumPy и SciPy применяются только при независимой проверке и построении отчёта.</p>
<p><code>-fullauto</code> и <code>--fullauto</code> равнозначны. Старый <code>--autofull</code> сохраняет прежний алгоритм; теперь это отдельный режим.</p>
<h3>Команда для GPU</h3>
<pre>./gpu/bin/mbs_po_gpu_double --method po --backend cuda \\
  --particle 1 446.7 138.9 --ri 1.3116 0 --wavelength-um 0.532 \\
  --max-reflections 12 --beam-cutoff 0.001 \\
  --fullauto-theta-range 0 25 -fullauto 0.01 \\
  --fullauto-gpus 0,1,2,3 --output results/fullauto_H446</pre>
<p>Это инструкция запуска, а не уже завершённый тест нового флага на H=446,7 мкм. <code>--particle 1 H D</code> задаёт шестигранную призму с высотой H и диаметром D, в мкм. <code>--ri</code> задаёт действительную и мнимую части показателя преломления. Длина волны также в мкм. В примере задана конечная глубина 12 и отсечение пучков 0,001. Диапазон θ задаётся без числа узлов: проверочная сетка генерируется по размеру и длине волны, а опорные узлы выбираются адаптивно.</p>
<p>Внутренний подбор φ заменяет число NPHI в <code>--scattering-grid</code>. Чтобы закрепить φ, задайте отдельно <code>--phi-points N</code>. <code>--fullauto-gpus</code> перечисляет локальные видимые GPU и учитывает <code>CUDA_VISIBLE_DEVICES</code>. Программа распределяет независимые задания по очереди: освободившаяся карта получает следующее. Между серверами задания распределяет внешний планировщик.</p>
<h3>Бюджет точности</h3>
<p>Интерполяции всегда выделяется <b>EPS/4</b>, статистической оценке — <b>3 EPS/4</b> без зеркального сокращения. При зеркальном сокращении EPS/256 резервируется для парного аудита, статистике остаётся 191 EPS/256. Для 1% интерполяция всегда ограничена 0,25%. Более мягкий допуск интерполяции установить нельзя. Для каждого угла проверяются все 16 элементов матрицы Мюллера; матрица преобразует четыре компоненты вектора Стокса падающего света в компоненты рассеянного света. Все абсолютные ошибки делятся на среднее M11 в том же угле. Это позволяет проверять M12 и элементы, проходящие через ноль.</p>
<p>По восьми независимым перемешиваниям Sobol программа вычисляет интервал <code>e = t₇ · s / (√8 · M11)</code>, где <code>t₇=2,364624251</code>, а s — выборочное стандартное отклонение восьми оценок с делителем 7. Это оценённый поточечный интервал 95%. Он не является одновременной гарантией 95% для всех углов и не оценивает отличие физической оптики от ADDA.</p>
<h3>Пилоты, Френель и ориентации</h3>
<p>На CUDA пять коротких запусков сравнивают блоки 64/128/256, конвейер обработки ориентаций и варианты организации пучков. Выбирается самый быстрый вариант, совпавший с базовым по всем элементам в пределах 10⁻⁷ от M11. На CPU этот подбор пропускается.</p>
<p>Затем независимые обучающие и проверочные пилоты сравнивают φ=1,2,4,8,16,32,64. Для каждого варианта оценивается произведение измеренного времени проверочных заданий на квадрат максимального интервала. Минимальное произведение задаёт выбор φ. Это прогноз стоимости сходимости, а не доказательство отсутствия смещения азимутальной квадратуры. Отдельная проверка с более плотной φ полезна при смене формы или физического режима.</p>
<p>Аналитические средние отражения от граней и круговой модели тени используются как контрольные переменные для M11. Формула: <code>Ycv = Y + bR(μR−CR) + bS(μS−CS)</code>. Y — исходная интенсивность; CR и CS — выборочные оценки контрольных компонентов; μR и μS — их известные средние; bR и bS — коэффициенты. Моменты Френеля — интегралы комбинаций коэффициентов отражения s- и p-поляризации по распределению углов падения. Они входят в μR. Само когерентное поле, включая интерференцию внутренних путей, остаётся в численном остатке.</p>
<p>Коэффициенты получают из ковариации отдельной обучающей выборки с регуляризацией 0,1 и границами −4…4. Отдельная проверочная выборка принимает их лишь при дисперсии менее 0,7 от варианта bR=bS=1. Рабочие seeds не участвуют в этом выборе. На новых θ коэффициенты детерминированно интерполируются и остаются фиксированными относительно рабочей выборки.</p>
<p>Обучение: 1009,1013,1019,1021,1031,1033,1039,1049. Проверка: 2003,2011,2017,2027,2029,2039,2053,2063. Расчёт: 11,23,37,53,71,89,107,131. Пересечение запрещено. <code>--sobol-seed N SEED</code> во внутренних командах означает N ориентаций и номер перемешивания. N в отчёте относится к одному seed, таких seeds восемь.</p>
<p><code>--haar-alpha</code> добавляет равномерный полный угол α к мере Хаара вращений. В этой схеме малое число φ дополняется независимым α. Оно не получается делением прежней регулярной сетки на число симметрий. Для призмы стандартный генератор сохраняет фундаментальный диапазон γ, который задаёт описание частицы. Для произвольной частицы используются её объявленные диапазоны; fullauto новых симметрий не предполагает.</p>
<h3>Уровни внутренних отражений</h3>
<p>Для конечной глубины n уровни: <code>max(1,n−4), max(1,n−2), n</code>, повторы удаляются. При n=12 оценка имеет вид <code>M8 + (M10−M8) + (M12−M10)</code>. В каждой разности используются те же ориентации, поэтому поправка обычно имеет меньшую дисперсию. Разные уровни могут использовать разные N. Среднее суммируется отдельно для каждого seed; ковариация между уровнями сохраняется при оценке общего интервала.</p>
<p>Для распределения бюджета используются измеренная стоимость cℓ одной ориентации и <code>Vℓ = Nℓ Var(оценки уровня)</code>. Предложение бюджета пропорционально <code>√(Vℓ/cℓ) · Σ√(Vj cj)</code> и обратно пропорционально квадрату допуска. Проверка берёт худший элемент и угол; бюджет округляется до степени двойки, за раунд растёт не более чем в четыре раза. Закон дисперсии 1/N служит прогнозом, фактическую остановку решает измеренный интервал. Нужны две последовательные успешные проверки и отдельная проверка уточнения основного уровня.</p>
<p>Глубина n и отсечение остаются параметрами конечной физической модели. Этот режим не доказывает сходимость по n→∞ и не подбирает допустимое смещение от cutoff. Отсечение действует во внутренних заданиях и может сокращать трассировку; его эффект на точность нужно проверять отдельно. Это действующий механизм, а не неиспользуемый флаг.</p>
<h3>Адаптивная θ и интерполяция</h3>
<p>Запрошенная сетка определяет строки итогового файла. Можно задать <code>--theta-grid-file</code>: числовые возрастающие углы от 0 до 180°, без заголовка. Либо <code>--scattering-grid TH1 TH2 NPHI NTH</code>. Без них для стандартной призмы или файла частицы строится оценочная сетка с шагом порядка λ/(16L), где L — консервативная граница диаметра частицы. Для остальных встроенных форм требуется явная сетка.</p>
<p>Программа выбирает стартовые опорные узлы. В каждом интервале добавляет прямые расчёты около ¼, ½ и ¾ его длины на запрошенной сетке. Кубический сплайн с граничным условием not-a-knot строится отдельно для каждого seed: первые два и последние два сегмента продолжают один кубический полином. Линейность этого сплайна по значениям сохраняет совместную статистику интерполированных элементов.</p>
<p>В контрольных узлах проверяется <code>(|среднее остатка| + t₇·sостатка/√8)/M11 ≤ EPS/4</code>. Не прошедшие узлы становятся опорными; контрольные узлы пересоздаются. Интервал проверяется и на всех восстановленных выходных углах. Отрицательное восстановленное M11 вызывает уточнение. Проверка конечного набора контрольных узлов не является доказательством максимальной ошибки между всеми узлами. В этом отчёте отдельные плотные расчёты дополнительно проверяют каждый запрошенный угол.</p>
<p>При изменении θ текущая версия запускает задания на новой общей сетке. Она переиспользует прежние угловые строки при совпадении N, seed, глубины, φ, симметрии и режима вычислений. Для новых θ выполняются отдельные задания. M11 пересчитывается из исходного Y и сохранённых контрольных компонентов относительно согласованных аналитических средних. Поэтому сокращение числа конечных θ не равно такому же ускорению всего адаптивного расчёта.</p>
<h3>Память, ограничения и возобновление</h3>
<p>N не означает, что все ориентации одновременно хранятся на GPU. Рабочий процесс обрабатывает ограниченные пакеты; <code>--orientation-chunk</code> задаёт их размер. Контроллер использует заданный размер пакета; штатный механизм рабочего процесса дополнительно ограничивает его по памяти. Меньшая θ-сетка уменьшает размер угловых рабочих массивов. На карте работает одно задание контроллера. Число CPU-потоков по умолчанию ограничено 16 на задание с учётом числа карт; <code>--threads</code> задаёт его явно.</p>
<table><tr><th>Флаг</th><th>Назначение / значение по умолчанию</th></tr>
<tr><td>--fullauto-pilot N</td><td>32768 ориентаций на каждый обучающий и проверочный seed</td></tr>
<tr><td>--fullauto-kernel-count N</td><td>8192 ориентации для подбора CUDA</td></tr>
<tr><td>--fullauto-initial N</td><td>131072 на seed основного уровня</td></tr>
<tr><td>--fullauto-min-correction N</td><td>8192 на seed каждой поправки</td></tr>
<tr><td>--fullauto-theta-start N</td><td>33 опорных узла, включая границы</td></tr>
<tr><td>--fullauto-max-rounds N</td><td>32 раунда; достижение лимита не означает сходимость</td></tr>
<tr><td>--max-orientations N</td><td>67108864 на seed и уровень</td></tr>
<tr><td>--max-theta-points N</td><td>Лимит реально вычисляемых узлов, включая контрольные</td></tr>
<tr><td>--max-phi-points N</td><td>Ограничивает проверяемые варианты φ; fullauto не требует кратности 6</td></tr>
<tr><td>--adaptive-config FILE</td><td>Может задать начальные значения и лимиты; его tolerance не меняют бюджеты EPS/4 и 3EPS/4</td></tr>
</table>
<p>Режим рассчитан на одну форму и один размер с когерентной PO. Внешние коэффициенты контроля, --incoh, --multisize, FFT и запуск нескольких MPI-рангов не поддерживаются. Для нескольких серверов запускайте локальные контроллеры на разных размерах.</p>
<p>Повторите ту же команду с той же папкой, чтобы возобновить расчёт. Контроллер проверяет исполняемый файл, аргументы, файлы геометрии и сетки, окружение MBS, GPU и имя сервера. Несовпадение требует новой папки. Совпавшие задания читаются из кэша с проверкой контрольных сумм DAT и контрольных компонентов. Блокировка предотвращает два контроллера в одной папке.</p>
<p>Код выхода <b>0</b> означает converged; <b>3</b> — requires_further_sampling при исчерпанном лимите; <b>2</b> — ошибка ввода, задания или кэша. После остановки проверяйте <code>fullauto_status.json</code>. В нём указаны фактически вычисленные бюджеты, а не ещё не выполненное предложение следующего раунда.</p>
<h3>Файлы результата</h3>
<p><code>mueller_fullauto.dat</code> содержит все запрошенные углы, а <code>mueller_evaluated.dat</code> — непосредственно вычисленные. Оба имеют 18 столбцов: θ, вес углового интервала и 16 элементов матрицы. <code>all_mueller_confidence.csv</code> содержит интервалы всех элементов; <code>theta_validation.csv</code> — ошибки в контрольных узлах. <code>theta_support.csv</code>, <code>theta_evaluated.csv</code> и <code>requested_theta.csv</code> задают три сетки.</p>
<p><code>calibration.json</code> записывает выбор φ и CUDA; <code>weights_*.tsv</code> — фиксированные коэффициенты; <code>refinements.tsv</code> — историю уточнений. В <code>jobs/</code> сохранены команды <code>command.args</code>, логи <code>stdout.log</code>, матрицы, аналитические компоненты и контрольные суммы. Кэш средних <code>mean_*.cache</code> сохраняется между заданиями.</p>
<h3>Настройки новых оптимизаций и произвольные формы</h3>
<p><code>--fullauto-mirror auto</code> — стандартный режим: проверяет отражение геометрии относительно плоскости xz и явно рассчитанные пары (β,Γ−γ,−α). Только затем использует половину γ и восстановление P M(−φ) P, P=diag(1,1,−1,−1). <code>on</code> требует успешной проверки; <code>off</code> отключает сокращение. Флаг <code>--mirror-gamma</code> равнозначен требованию on. Изменение знаков применяется ко всем элементам, а не только к M11.</p>
<p><code>--fullauto-angular-cache on</code> включён по умолчанию. <code>off</code> позволяет независимо повторить расчёт без объединения строк. Собранные результаты помечены <code>angular_reuse.json</code>; <code>angular_sources.tsv</code> содержит источники с контрольными суммами. Их command.args — рецепт полного прямого расчёта, а stdout.log явно сообщает, что запуск был заменён сборкой. Все действительные запуски новых строк сохранены отдельными заданиями.</p>
<p>Аналитические средние вычисляются один раз на всей проверочной сетке через <code>--analytic-mean-reference-grid</code>, в том числе при числе узлов больше4097. Уточнение θ читает тот же кэш: старые и новые строки используют одинаковые μR и μS. Расход памяти на сохранённые средние растёт линейно с числом строк; радиальная таблица зависит от максимального угла и размера грани. В прежней версии большие сетки требовали повторной подготовки средних для каждого нового набора θ. В распределении бюджета используются прогнозы стоимости полной сетки по измеренным прямым заданиям; фактическая остановка всегда определяется матрицами и интервалами.</p>
<p><code>--fullauto-analytic-controls auto</code> применяет контрольные компоненты при допустимой памяти и сходимости средних. При превышении радиальной таблицы или непригодной скорректированной оценке автоматически переходит к полному выборочному PO и сбрасывает проверки сходимости. <code>on</code> требует контрольных компонентов, <code>off</code> отключает их. Это меняет способ оценки и дисперсию, сохраняя полный когерентный физический расчёт, заданные глубину и cutoff.</p>
<p>Для произвольной формы: <code>--particle-file shape.particle --symmetry 1 1 --fullauto 0.01</code>. Полная область β=0…180°, γ=0…360°, α=0…360° используется без неподтверждённых вращательных симметрий. Зеркальный режим остаётся условным после проверки. Неплоские грани и другие ошибки сетки отклоняются; программа не подправляет геометрию молча.</p>
<p>Равномерные φ вместе со случайным α образуют квадратуру азимутального интеграла со случайным сдвигом. <code>--alpha-points</code> — синоним <code>--phi-points</code>; это плотность квадратуры, а не отдельное число значений третьей координаты Sobol. В ровно0° и180° поворот поляризационного базиса дополнительно усредняется аналитически функциями PoleMueller.</p>
<p><code>--fullauto-theta-range TH1 TH2</code> задаёт только область углов, без числа точек. В автоматическом режиме до производства проводится независимый плотный пилот по всем точкам сгенерированной проверочной сетки. Он выбирает опорные узлы с пилотным допуском EPS/8. Число ориентаций на seed: min(pilot,8192), меняется через <code>--fullauto-theta-audit-count</code>; 0 отключает дополнительный пилот. Рабочие проверки интерполяции остаются EPS/4. Проверки конечной сетки не доказывают ошибку в непрерывном промежутке между её точками.</p>
<h3>Сравнение с независимым Python-аудитом</h3>
<p>Аудит читает исходные DAT, восстанавливает парные разности уровней и итог по восьми seeds. Он вызывает прежние Python-функции confidence, fit_control_weights, estimate_cost_score и allocation; сплайн строит отдельно через SciPy CubicSpline. Совпали выбранное φ, контрольные коэффициенты и предложение распределения N. Сравнение использует одинаковые задания и измеренные времена: старый Python-контроллер не имел этой адаптивной θ-сетки.</p>
<table><tr><th>Тест</th><th>n</th><th>φ</th><th>N по уровням на seed</th><th>θ: запрос → прямые / опорные</th><th>Статистика</th><th>Контроль интерполяции</th><th>Плотный расчёт: отклонение</th><th>Остаток +95%</th></tr>'''+records+'''</table>
<p>Первая проверка: H=0,2 мкм, D=0,1 мкм, λ=1 мкм, m=1,31, θ=0…25°. Это проверка реализации, а не физической точности PO для малой частицы. Проверка depth12, если приведена в таблице: H=10 мкм, D=3,1 мкм, λ=0,532 мкм, m=1,3116, θ=170…180°, cutoff=0,001.</p>
<p>Максимальное расхождение восстановленных средних Python и C++, нормированное на M11: '''+f'{max_difference:.3g}'+'''. Отдельно проверены возобновление без новых заданий, отказ при повреждённом кэше и выход 3 при недостигнутой точности. CUDA FP64 собрана и проверена на epyc1 и epyc2: восемь независимых seeds, все 16 элементов и допуск 1%. Новые большие расчёты запущены отдельной кампанией. Ускорение полной GPU-кампании ещё не измерено.</p>
<pre>bash tests/run_fullauto_tests.sh
python3 tests/audit_fullauto.py PATH_TO_COMPLETED_OUTPUT --dense</pre>
'''+plots


def optimization_body(mean_validation, arbitrary_validation):
    content = '<h2>Ускорение подготовки средних и проверка произвольных форм</h2>'
    content += '<p>Аналитические средние μR и μS — средние контрольных компонентов отражения и тени по ориентациям. Их подготовка распараллелена по θ с помощью OpenMP. Число потоков задаёт <code>--threads N</code>. Каждая строка сохраняет порядок суммирования, а решение о сходимости принимается после последовательного объединения строк. Этот приём применим к произвольным граням и не использует симметрию призмы.</p>'
    if mean_validation:
        records = json.loads((mean_validation/'comparison.json').read_text())
        content += '<table><tr><th>Контрольная грань</th><th>Строк θ</th><th>1 поток, с</th><th>8 потоков, с</th><th>Ускорение подготовки</th><th>Побитное совпадение</th></tr>'
        for row in records:
            assert row['bitwise_identical']
            content += '<tr>'+''.join('<td>'+html.escape(str(value))+'</td>' for value in (
                row['case'], row['rows'], f"{row['serial_seconds']:.4f}",
                f"{row['parallel_seconds']:.4f}", f"{row['speedup']:.3f}×", 'да'))+'</tr>'
        content += '</table><p>Это ускорение только холодной подготовки средних. Полное время включает пилоты, трассировку, дифракционный интеграл и аудит; ускорение всей программы из этой таблицы не следует.</p>'
    if arbitrary_validation:
        records = json.loads((arbitrary_validation/'comparison.json').read_text())
        content += '<p>Следующие проверки выполнены новой версией на CPU, на файлах частиц с <code>--symmetry 1 1</code>, допуском 1%, λ=1 мкм и m=1,3116. θ выбирается автоматически по диапазону 0…180°. В проверке φ закреплено равным 4. Выпуклый тест — повёрнутый параллелепипед 1×0,8×0,6 мкм; вогнутый — <code>examples/particles/concave_hexagonal.particle</code>. Это проверка реализации и точности усреднения относительно прямого PO.</p>'
        content += '<table><tr><th>Форма</th><th>n</th><th>Время контроллера, с</th><th>θ: запрос / прямые</th><th>Статистика</th><th>Плотный остаток +95%</th><th>Зеркальное сокращение</th><th>Контрольные компоненты</th></tr>'
        plots = ''
        for row in records:
            state, audit = row['state'], row['audit']
            assert audit['independent_dense']['all_requested_nodes_pass']
            content += '<tr>'+''.join('<td>'+html.escape(str(value))+'</td>' for value in (
                row['case'], state['target_depth'], f"{row['controller_seconds']:.2f}",
                f"{state['requested_theta']} / {state['evaluated_theta']}",
                f"{100*state['max_pointwise95_scaled_M11']:.6f}%",
                f"{100*audit['independent_dense']['residual_plus95_scaled_M11']:.6f}%",
                'да' if state['mirror_gamma'] else 'нет',
                'да' if state['analytic_controls'] else 'нет'))+'</tr>'
            rows = np.loadtxt(arbitrary_validation/row['case']/'python_vs_native.csv', delimiter=',', skiprows=1)
            fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.7))
            axes[0].plot(rows[:,0], rows[:,1], label='fullauto')
            axes[0].plot(rows[:,0], rows[:,2], '--', label='direct theta grid')
            axes[0].set_ylabel('M11')
            axes[1].plot(rows[:,0], rows[:,3]/rows[:,1], label='fullauto')
            axes[1].plot(rows[:,0], rows[:,4]/rows[:,2], '--', label='direct theta grid')
            axes[1].set_ylabel('M12 / M11')
            axes[2].plot(rows[:,0], 100*rows[:,5], label='max difference, all 16 entries')
            axes[2].plot(rows[:,0], 100*rows[:,6], label='max statistical 95% interval')
            axes[2].axhline(.25, linestyle=':', color='black', label='interpolation budget')
            axes[2].set_ylabel('% of M11')
            for ax in axes:
                ax.set_xlabel('theta, deg'); ax.grid(alpha=.2); ax.legend(fontsize=7)
            plots += figure(fig, row['case']+': реально вычисленные матрицы новой версии и независимый прямой расчёт всех запрошенных θ.')
        content += '</table><p>Зеркальное сокращение разрешается только после геометрической и численной проверки пар по всем 16 элементам. Отказ проверки включает полный диапазон γ. Время в таблице не включает последующий независимый плотный аудит. Сравнение с ADDA в этих тестах не выполнялось.</p>'+plots
    return content


def automatic_phi_body(validation):
    records = json.loads((validation/'comparison.json').read_text())
    content = '<h2>Окончательная версия: произвольные формы и автоматический φ</h2>'
    content += '<p>Эти отдельные тесты выполнены на CPU тем же исходным кодом контроллера, который работает в текущей GPU-кампании. Допуск EPS=1%, λ=1 мкм, все три угла ориентации покрывают полную область: β=0…180°, γ=0…360°, α=0…360°. Число φ и опорные θ программа выбирает автоматически. Тетраэдр имеет неравные рёбра и смещённую вершину; второй тест использует повёрнутый параллелепипед с поглощением. Ни у одной из этих форм геометрическая проверка не разрешила зеркальное сокращение.</p>'
    content += '<table><tr><th>Форма</th><th>m</th><th>n</th><th>Выбранное φ</th><th>θ: запрос / прямые</th><th>Статистика, % M11</th><th>Плотный остаток +95%, % M11</th><th>Расхождение с прежней версией / M11</th></tr>'
    plots = ''
    for row in records:
        state, report = row['state'], row['audit']
        assert row['automatic_phi'] and report['passed']
        assert not state['mirror_gamma'] and report['independent_dense']['all_requested_nodes_pass']
        previous = row['previous_implementation']
        assert previous['audit']['passed']
        content += '<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in (
            row['case'], row['refractive_index'], state['target_depth'], state['phi'],
            f"{state['requested_theta']} / {state['evaluated_theta']}",
            f"{100*report['statistical95_scaled_M11']:.6f}",
            f"{100*report['independent_dense']['residual_plus95_scaled_M11']:.6f}",
            f"{previous['max_matrix_difference_scaled_M11']:.3g}"))+'</tr>'
        rows = np.loadtxt(validation/row['result_directory']/'python_vs_native.csv', delimiter=',', skiprows=1)
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.7))
        for j, ylabel in [(0, 'M11'), (1, 'M12 / M11')]:
            native = rows[:,1] if j == 0 else rows[:,3]/rows[:,1]
            direct = rows[:,2] if j == 0 else rows[:,4]/rows[:,2]
            axes[j].plot(rows[:,0], native, label='fullauto')
            axes[j].plot(rows[:,0], direct, '--', label='direct theta grid')
            axes[j].set_ylabel(ylabel)
        axes[2].plot(rows[:,0], 100*rows[:,5], label='max difference, all 16 entries')
        axes[2].plot(rows[:,0], 100*rows[:,6], label='max statistical 95% interval')
        axes[2].axhline(.25, color='black', linestyle=':', label='interpolation budget')
        axes[2].set_ylabel('% of M11')
        for ax in axes:
            ax.set_xlabel('theta, deg'); ax.grid(alpha=.2); ax.legend(fontsize=7)
        plots += figure(fig, row['case']+': реальные матрицы при автоматическом выборе φ, прямой расчёт всех θ и ошибки по всем 16 элементам.')
    content += '</table><p>Прежняя версия проверена при том же выбранном φ; совпадение числа ориентаций сохранено в JSON. Сравнение времени этих двух запусков не является измерением ускорения: новая версия дополнительно выбирала φ. Независимый NumPy/SciPy-аудит восстановил средние, интервалы, коэффициенты контроля и распределение работы, затем отдельно рассчитал каждую точку полной сетки θ. Плотный остаток проверен против EPS/4=0,25%. Проверка относится к численной реализации PO; она не задаёт физическую погрешность относительно ADDA.</p>'
    content += '<pre>python3 tests/test_fullauto_arbitrary.py --binary bin/mbs_po \\\n  --case asymmetric --phi-points 0 --threads 4 \\\n  --reference-binary /path/to/previous/mbs_po \\\n  --output results/arbitrary_auto_phi</pre>'
    return content+plots


def build(validation, output, mean_validation=None, arbitrary_validation=None, shared_mean_validation=None, large_mean_validation=None, automatic_phi_validation=None):
    content = body(validation)
    if mean_validation or arbitrary_validation:
        content += optimization_body(mean_validation, arbitrary_validation)
    if automatic_phi_validation:
        content += automatic_phi_body(automatic_phi_validation)
    if shared_mean_validation:
        evidence = json.loads((shared_mean_validation/'comparison.json').read_text())
        assert evidence['same_counts']
        assert all(c['audit']['independent_dense']['all_requested_nodes_pass'] for c in evidence['cases'])
        residual = max(c['audit']['independent_dense']['residual_plus95_scaled_M11'] for c in evidence['cases'])
        content += '<h3>Общий кэш на большой проверочной сетке</h3><p>Отдельный тест сравнил прежний кэш для каждого набора θ и новый единый кэш на4099 запрошенных углах. Оба расчёта использовали H=0,2 мкм, D=0,1 мкм, λ=1 мкм, m=1,31, n=1, φ=4, EPS=1% и одинаковые seeds и N. Это тест реализации большого кэша на явной густой сетке; рабочая кампания задаёт диапазон без числа θ.</p>'
        content += '<p>Максимальное расхождение всех16 элементов, нормированное на M11: <b>'+f"{evidence['max_matrix_difference_scaled_M11']:.3g}"+'</b>. Независимый прямой расчёт всех4099 углов дал остаток с оценкой95% <b>'+f'{100*residual:.6f}%'+'</b>, ниже бюджета интерполяции0,25%. Все задания новой версии с контрольными компонентами используют один <code>means_shared.cache</code>; повторная подготовка для новых подмножеств θ отсутствует. Время этого малого теста не характеризует ускорение больших частиц: подготовка полной таблицы выполняется заранее, а её преимущество возникает при повторном использовании.</p>'
    if large_mean_validation:
        equality=json.loads((large_mean_validation/'large_mean_equivalence_epyc2.json').read_text())
        timing=json.loads((large_mean_validation/'large_mean_timing_epyc2.json').read_text())
        assert equality['bitwise_equal'] and equality['new_sha256']==equality['old_sha256']
        content += '<h3>Измерение на большой частице текущей кампании</h3><p>На epyc2 для призмы H=501,2 мкм, D=155,8 мкм, λ=0,532 мкм и m=1,3116 сравнили подготовку одного и того же массива из6889 аналитических средних. Предыдущая версия считала строки последовательно; новая использует32 CPU-потока. Контрольная сумма всего кэша совпала: средние и оценка их внутренней сходимости совпадают побитно.</p>'
        content += '<table><tr><th>Этап</th><th>Прежняя версия</th><th>Новая версия</th><th>Ускорение</th></tr><tr><td>Подготовка средних, с</td><td>'+f"{timing['old']['seconds']:.2f}"+'</td><td>'+f"{timing['new']['seconds']:.2f}"+'</td><td>'+f"{timing['cold_preparation_speedup']:.2f}×"+'</td></tr></table><p>Время взято из строк <code>Analytic facet average ... setup=</code> двух действительных запусков. Это время подготовки средних, а не полное время рассеяния и не окончательное ускорение всей кампании. Дополнительная экономия возникает потому, что новая версия использует этот массив повторно при уточнении θ.</p>'
    style = 'body{font:17px/1.65 system-ui,sans-serif;color:#172536;background:#f3f5f9}main{max-width:1180px;margin:auto;padding:28px;background:white}h2,h3{color:#183454}pre{background:#13263b;color:#e5effb;padding:18px;overflow:auto;font-size:14px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:9px;border-bottom:1px solid #ccd5df;text-align:left}figure{margin:28px 0}figure img{width:100%}figcaption{font-size:14px;color:#506276}code{background:#edf2f8;padding:2px 4px}'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('<!doctype html><html lang="ru"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MBS-fast: native fullauto и независимый аудит</title><style>'+style+'</style></head><body><main>'+content+'</main></body></html>')
    return content


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--mean-validation', type=Path)
    parser.add_argument('--arbitrary-validation', type=Path)
    parser.add_argument('--shared-mean-validation', type=Path)
    parser.add_argument('--large-mean-validation', type=Path)
    parser.add_argument('--automatic-phi-validation', type=Path)
    args = parser.parse_args()
    build(args.validation, args.output, args.mean_validation, args.arbitrary_validation, args.shared_mean_validation, args.large_mean_validation, args.automatic_phi_validation)
