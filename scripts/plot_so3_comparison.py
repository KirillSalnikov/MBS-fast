#!/usr/bin/env python3
"""Standalone HTML and PNG figures from completed scattering ensembles."""
import argparse
import base64
import html
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter

LABELS={'sobol':'Sobol','lebedev':'Лебедев × α','wigner':'Вигнер (Gauss × Fourier)'}
COLORS={'sobol':'#2067b2','lebedev':'#ca7324','wigner':'#2a9260'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--backend',choices=['cpu','gpu'],default='cpu')
    a=p.parse_args();root=a.root.resolve();report=json.loads((root/('comparison_'+a.backend+'.json')).read_text())
    records=report['completed_ensembles'];figures=[]
    out=root/('plots_'+a.backend);out.mkdir(exist_ok=True)
    names={'prism_1um':'Призма H=D=1 мкм','asymmetric_1um':'Асимметричный тетраэдр (~1 мкм)',
           'prism_H501':'Призма H=501,2; D=155,8 мкм','prism_H1000':'Призма H=1000; D=220,1 мкм'}
    rows=[]
    for case in sorted({r['case'] for r in records}):
        subset=[r for r in records if r['case']==case and r['level']!='reference']
        fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
        for method in LABELS:
            values=sorted([r for r in subset if r['method']==method],key=lambda r:r['orientations_per_replica'])
            if not values:continue
            y=[r['statistical95_scaled_M11']*100 for r in values]
            for ax,key in zip(axes,['orientations_per_replica','total_seconds']):
                ax.loglog([r[key] for r in values],y,'o-',color=COLORS[method],label=LABELS[method])
        for ax in axes:
            ax.axhline(1,color='black',linestyle='--',linewidth=1,label='1%')
            ax.set_ylabel('Максимальная оценка 95%, % от M11');ax.grid(True,which='both',alpha=.2);ax.legend(fontsize=8)
        axes[0].set_xlabel('Ориентаций на один повтор');axes[1].set_xlabel('Суммарное время восьми повторов, с')
        if subset:
            times=[r['total_seconds'] for r in subset]
            axes[1].xaxis.set_major_locator(FixedLocator(np.geomspace(min(times),max(times),4)))
            axes[1].xaxis.set_major_formatter(FuncFormatter(lambda x,pos:'{:.3g}'.format(x)))
            axes[1].xaxis.set_minor_formatter(NullFormatter())
        fig.suptitle(names.get(case,case)+'; '+a.backend.upper())
        path=out/(case+'_convergence.png');fig.savefig(path,dpi=150);plt.close(fig);figures.append(path)
        for r in subset:
            difference=r.get('difference_to_reference_scaled_M11')
            rows.append('<tr><td>{}</td><td>{}</td><td>{}</td><td>{:.3f}</td><td>{}</td><td>{:.2f}</td></tr>'.format(
                names.get(case,case),LABELS[r['method']],r['orientations_per_replica'],100*r['statistical95_scaled_M11'],
                '—' if difference is None else '{:.3f}'.format(100*difference),r['total_seconds']))
        matrices={}
        for method in LABELS:
            level='31';items=[]
            for marker in (root/'results'/a.backend).glob('*/complete.json'):
                j=json.loads(marker.read_text())['job']
                if j['case']==case and j['method']==method and j['level']==level:
                    items.append(np.loadtxt(marker.parent/'matrix/matrix.dat',skiprows=1))
            if len(items)==8:matrices[method]=np.mean(items,axis=0)
        refs=[]
        for marker in (root/'results'/a.backend).glob('*/complete.json'):
            j=json.loads(marker.read_text())['job']
            if j['case']==case and j['level']=='reference':refs.append(np.loadtxt(marker.parent/'matrix/matrix.dat',skiprows=1))
        if len(refs)==8 and matrices:
            ref=np.mean(refs,axis=0);fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
            for col,title,ax in zip([0,1,5],['M11','M12/M11','M22/M11'],axes):
                ax.plot(ref[:,0],ref[:,2] if col==0 else ref[:,col+2]/ref[:,2],'k.--',label='Sobol reference')
                for method,data in matrices.items():
                    ax.plot(data[:,0],data[:,2] if col==0 else data[:,col+2]/data[:,2],'o-',
                            markersize=3,color=COLORS[method],label=LABELS[method])
                if col==0:ax.set_yscale('log')
                ax.set_title(title);ax.set_xlabel('θ, градусы');ax.grid(alpha=.2);ax.legend(fontsize=7)
            fig.suptitle(names.get(case,case)+'; квадратуры степени 31, только рассчитанные углы')
            path=out/(case+'_mueller.png');fig.savefig(path,dpi=150);plt.close(fig);figures.append(path)
    page='''<!doctype html><html lang="ru"><meta charset="utf-8"><title>Сравнение усреднения SO(3)</title>
<style>body{font:16px system-ui;max-width:1200px;margin:30px auto;padding:0 20px;color:#233}table{border-collapse:collapse}td,th{border:1px solid #ccd;padding:8px}img{max-width:100%}p{line-height:1.5}</style>
<h1>Реальное сравнение Sobol, Лебедева и Вигнера</h1>
<p>Восемь независимых повторов: Owen-scramble для Sobol и глобальные равномерные повороты квадратуры для двух других методов. Полная мера Haar SO(3); симметрия не сокращалась. α входит в каждую ориентацию, φ=4 для всех методов. λ=0,532 мкм, m=1,3116+0i, выходная отсечка 0,001; глубина 6 для малых форм и 12 для больших призм. Для малых форм рассчитаны 13 углов, включая 170–180°; для больших — 11 углов в 0–25°.</p>
<p>«Оценка 95%» — t₇·s/√8 для среднего восьми повторов, максимум по рассчитанным θ и всем 16 элементам, нормированный на M11. Это оценка статистики в проверенных точках. Плотность θ, глубина и отсечка в этом опыте не уточнялись. Reference — отдельный Sobol 65536×8; его конечная погрешность также приведена в JSON. Соединительные линии на графиках не являются новыми расчётами.</p>
<p>Вигнер здесь означает выделение среднего (постоянной моды D⁰₀₀) квадратурой Gauss–Fourier. Полный спектр коэффициентов Dʲₘₙ не восстанавливается. Алгебраическая точность проверена отдельно до степени 131; такая точность не означает точность 1% реального рассеяния.</p>
'''
    page+='<p>Время: '+a.backend.upper()+', сумма времени процессов, включая подготовку и вывод. CPU и GPU необходимо сравнивать отдельно.</p>'
    page+='<table><tr><th>Форма</th><th>Метод</th><th>N на повтор</th><th>Оценка 95%, %</th><th>Разница с reference, %</th><th>Время, с</th></tr>'+''.join(rows)+'</table>'
    for path in figures:
        encoded=base64.b64encode(path.read_bytes()).decode('ascii')
        page+='<p><img src="data:image/png;base64,'+encoded+'" alt="Рассчитанное сравнение"></p>'
    page+='<p><a href="comparison_'+a.backend+'.json">Полные метрики и ограничения</a>; <a href="manifest.json">Все параметры запусков</a>; исходные матрицы, команды и логи находятся в results/'+a.backend+'/.</p></html>'
    target=root/('Сравнение_усреднения_'+a.backend+'.html');target.write_text(page)
    print(target)


if __name__=='__main__':main()
