#!/usr/bin/env python3
"""Standalone scientific plots and HTML from actual native tuning test data."""
import argparse,base64,html,json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();root=a.root.resolve()
    plots=root/'plots';plots.mkdir(exist_ok=True);parts=[]
    def figure(name,caption):
        path=plots/(name+'.png');plt.savefig(path,dpi=170,bbox_inches='tight');plt.close()
        data=base64.b64encode(path.read_bytes()).decode()
        parts.append('<figure><img src="data:image/png;base64,'+data+'"><figcaption>'+html.escape(caption)+'</figcaption></figure>')
    labels={'prism':'Hexagonal prism','asymmetric':'Asymmetric tetrahedron','absorbing':'Absorbing rotated box','concave':'Concave prism'}
    fig,axes=plt.subplots(4,3,figsize=(13,12))
    for row,name in enumerate(labels):
        for mode,color in [('off','black'),('auto','#087fbb')]:
            matrix=np.loadtxt(root/'cpu_amortized'/(name+'_'+mode)/'mueller_fullauto.dat',skiprows=1)
            theta=matrix[:,0];m=matrix[:,2:]
            values=[m[:,0],m[:,1]/m[:,0],m[:,5]/m[:,0]]
            for col,y in enumerate(values):axes[row,col].plot(theta,y,color=color,label=mode,linewidth=1.8)
        for col in range(3):axes[row,col].grid(alpha=.25);axes[row,col].set_xlabel('theta, deg');axes[row,col].set_ylabel(['M11','M12 / M11','M22 / M11'][col])
        axes[row,0].set_title(labels[name]);axes[row,0].legend()
    fig.tight_layout();figure('mueller_four_shapes','Реальные расчёты четырёх форм: auto и контроль off. Оба запуска используют тот же код PO; это не сравнение с ADDA. λ=1 мкм, заданный допуск 1%.')
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for ax,name in zip(axes.flat,labels):
        folder=root/'cpu_amortized'/(name+'_auto');raw=np.loadtxt(folder/'all_mueller_confidence.csv',delimiter=',',skiprows=1)
        ratios=np.loadtxt(folder/'normalized_mueller_confidence.csv',delimiter=',',skiprows=1)
        state=json.loads((folder/'fullauto_status.json').read_text())
        ax.plot(raw[:,0],raw[:,1:].max(1)*100,label='all Mij / mean M11')
        ax.plot(ratios[:,0],ratios[:,1:].max(1)*100,label='ratios Mij / M11')
        ax.axhline(state['statistical_budget']*100,color='red',ls='--',label='statistical budget')
        ax.set_title(labels[name]);ax.set_xlabel('theta, deg');ax.set_ylabel('pointwise 95% estimate, %');ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.tight_layout();figure('precision','Максимальная оценка статистической ошибки по 16 элементам для каждого θ. Отношения имеют абсолютный допуск; все интервалы точечные, не одновременная гарантия диапазона.')
    case=root/'cutoff_back30_final';state=json.loads((case/'fullauto_status.json').read_text());saved=json.loads((case/'automatic_cutoff_validation.json').read_text())
    first=1600000000+state['round']*1024;seeds=[first+17*s for s in range(8)];raw={'off':{},saved['profile']:{}};angles=None
    for folder in (case/'jobs').iterdir():
        if not (folder/'complete.txt').exists():continue
        args=(folder/'command.args').read_text().splitlines()
        def option(name):return args[args.index(name)+1] if name in args else None
        if '--sobol-seed' not in args or option('--max-reflections')!='30':continue
        i=args.index('--sobol-seed');count,seed=map(int,args[i+1:i+3]);profile=option('--cutoff-profile')
        if seed in seeds and count==saved['count_per_seed'] and profile in raw:
            data=np.loadtxt(folder/'result/result.dat',skiprows=1);angles=data[:,0];raw[profile][seed]=data[:,2:]
    means={key:np.array([value[s] for s in seeds]).mean(0) for key,value in raw.items()}
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for profile,color in [('off','black'),(saved['profile'],'#c65f10')]:
        m=means[profile]
        for ax,value,title in zip(axes,[m[:,0],m[:,1]/m[:,0],m[:,5]/m[:,0]],['M11','M12 / M11','M22 / M11']):
            ax.plot(angles,value,'o-',color=color,label=profile);ax.set_title(title);ax.set_xlabel('theta, deg');ax.grid(alpha=.25);ax.legend()
    fig.tight_layout();figure('paired_cutoff_backscatter','Окрестность 180°: n=30, правильная призма H=0,2 мкм, D=0,1 мкм, λ=1 мкм, m=1,3116. Два полных когерентных расчёта используют одинаковые 128 ориентаций на каждый из восьми новых seed. Разность с учётом статистики — '+format(saved['bias95']*100,'.6f')+'%, бюджет — 0,25%.')
    calibration=json.loads((case/'cutoff_calibration.json').read_text());rows=calibration['candidates']
    plt.figure(figsize=(8,4));plt.bar([r['profile'] for r in rows],[r['seconds'] for r in rows],color=['#444','#8db','#4aa','#c65f10']);plt.ylabel('represented multilevel native seconds');plt.grid(axis='y',alpha=.25)
    figure('cutoff_cost','Измеренная стоимость работы уровней 26, 28, 30 на пилоте. В off учтена совместная трассировка; при отсечении ветвей она отключена. Это стоимость части расчёта, не ускорение всей программы.')
    cpu=json.loads((root/'cpu_amortized/comparison.json').read_text());gpu_file=root/'gpu_complete_comparison.json';gpu=json.loads(gpu_file.read_text()) if gpu_file.exists() else []
    text='''<!doctype html><html lang="ru"><meta charset="utf-8"><title>Fullauto: автоматические флаги и проверки</title><style>body{font:17px system-ui;line-height:1.5;max-width:1200px;margin:30px auto;padding:20px;color:#20303c}img{max-width:100%}table{border-collapse:collapse;width:100%}td,th{padding:8px;border:1px solid #ccd}pre{background:#edf3f7;padding:18px;overflow:auto}figure{margin:30px 0}figcaption{color:#435}</style>
    <h1>Fullauto: автоматические флаги и проверка результата</h1>
    <p>EPS — допуск из --fullauto. N — число ориентаций в каждой независимой выборке seed; программа использует восемь таких выборок. n — максимальное число внутренних столкновений. θ — полярный угол рассеяния; φ — число азимутальных направлений наблюдения. α — третий угол Эйлера, поворот частицы вокруг направления падающего света. Mij — элементы матрицы Мюллера; M11 задаёт масштаб интенсивности, M12 и M22 описывают поляризационный отклик.</p>
    <pre>bin/mbs_po --method po --backend cuda \\
  --particle 1 10 6 --ri 1.3116 0 --wavelength-um .532 \\
  --max-reflections 12 --fullauto .01 --fullauto-gpus auto \\
  --output results/new_auto</pre>
    <p>Три новые политики по умолчанию: fullauto-resources auto, fullauto-rotation-symmetry auto, fullauto-cutoffs auto. Явные потоки, пакеты, блоки, симметрия и отсечения имеют приоритет. Без fullauto-gpus auto сохраняется GPU 0. CUDA_VISIBLE_DEVICES ограничивает разрешённые карты.</p>
    <h2>Что делает программа</h2><p>Она измеряет нативные фазы — трассировку, подготовку пучков и вычисление рассеянного поля. Затем проверяет, окупятся ли дополнительные пилоты потоков, пакетов и отсечений. Для короткой задачи дорогой перебор пропускается. Размер сохраняемого блока оценивается по уже измеренному времени запуска и вычислений. Рабочий пакет дополнительно ограничивается текущей RAM и VRAM. RAM разделяется между параллельными работниками; конвейер и два набора глубин учитываются в памяти.</p>
    <p>Размеры рабочих пакетов и сохраняемых блоков различаются. Пакет находится в памяти. Сохраняемый блок — часть последовательности Sobol, после которой готовый результат можно восстановить. Большой сохраняемый блок уменьшает число запусков и файловых операций; маленький уменьшает работу, потерянную при аварии.</p>
    <h2>Вращательная симметрия</h2><p>Программа сравнивает полные ориентированные грани после поворотов вокруг z частицы и полуповорота вокруг x. Проверяет фиксированные оси и несколько малых порядков вращения; это не поиск всей произвольной группы симметрий. Агрегаты и специальные флаги видимости исключены. Заголовок файла не служит доказательством. Затем явно трассируются эквивалентные повороты и сравниваются все 16 элементов и отношения.</p>
    <pre>F(α,β,γ) = F(α,β,γ + 2π/q)       при подтверждённом вращении порядка q
F(α,β,γ) = F(α+π,π−β,−γ)          при подтверждённом полуповороте вокруг x</pre>
    <p>У правильной призмы проверены q=6 и полуповорот: область γ сокращается в 6 раз, область cos β — в 2 раза. α остаётся полным. Это уменьшает область интеграла; число ориентаций до заданного допуска программа всё равно определяет по статистике. У несимметричной формы область не сокращается.</p>
    <h2>Отсечения и допуск</h2><p>Без явных порогов сравниваются off, safe, balanced, fast. off задаёт эталон при том же n. Пилот ограничен 2048 ориентациями на seed. Измеренное время учитывает потерю совместной трассировки при отсечении ветвей. Выбранный профиль должен проходить EPS/8 и давать хотя бы 10% выигрыша в пилоте. На короткой задаче перебор может быть пропущен.</p>
    <p>Перед завершением выполняется новый парный аудит: одинаковые ориентации для выбранного профиля и off, отдельные seed. Проверяются все элементы и отношения средних Mij/M11. Для малой разности сначала берётся небольшая выборка, которая удваивается при необходимости. При отказе профиль заменяется на off, и адаптация продолжается. n сохраняется; отдельная fullauto-physics-audit on проверяет чувствительность к n+2 и отсечениям. Автоматического подбора n сейчас нет.</p>
    <p>Бюджеты: интерполяция EPS/4; численная подготовка аналитических средних EPS/64; проверенное зеркало EPS/256; проверенные вращения EPS/256; настройка ресурсов/GPU EPS/1024; выбранные автоматические отсечения EPS/8. Неиспользуемые резервы остаются статистике. Это точечные статистические оценки и внутренние проверки, не строгая глобальная гарантия и не ошибка PO относительно ADDA.</p>
    <h2>Реальные сравнения</h2>'''
    text+=''.join(parts)
    text+='<h2>Время полного запуска</h2><p>В таблице время включает подбор и контроль до завершения, без последующего независимого Python-аудита. Малые задачи чувствительны к накладным расходам; эти однократные тесты не задают универсальное ускорение больших частиц. CPU-серия выполнялась при других локальных проверках. GPU-пары выполнены последовательно на соответствующем сервере.</p><table><tr><th>Серия</th><th>Форма</th><th>auto, с</th><th>off, с</th><th>off/auto</th><th>Разность всех Mij / M11</th></tr>'
    for series,records in [('CPU',cpu),('GPU',gpu)]:
        for r in records:
            if not (r.get('mode')=='auto' or r.get('case','').endswith('_auto')):continue
            name=r['case'].removesuffix('_auto');other=next(v for v in records if (v.get('case')==name and v.get('mode')=='off') or v.get('case')==name+'_off')
            text+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [series,name,round(r['seconds'],3),round(other['seconds'],3),round(r['auto_vs_reference_speedup'],3),format(r['all16_difference_scaled_M11'],'.4g')])+'</tr>'
    text+='''</table><h2>Файлы и воспроизведение</h2><p>В архиве сохранены все выданные матрицы, логи, команды, промежуточные попытки, замороженные бинарные файлы и исходники. Решения записываются в resource_calibration.json, rotation_verification.json, cutoff_calibration.json и automatic_cutoff_validation.json. execution_tuning.tsv восстанавливает точный выбор без нового подбора. Принудительное завершение и повторный запуск сохранили 86 готовых работ; повреждение итоговой матрицы было отвергнуто.</p><pre>python3 tests/test_fullauto_tuning.py --binary bin/mbs_po --output tuning_checks
python3 tests/audit_fullauto_tuning.py RESULT_DIRECTORY</pre></html>'''
    (root/'Fullauto_флаги_ускорения_и_проверки.html').write_text(text)

if __name__=='__main__':main()
