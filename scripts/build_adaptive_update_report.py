#!/usr/bin/env python3
"""Embed measured fullauto acceptance figures in a standalone HTML report."""
import argparse,base64,json,html
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);a=p.parse_args();root=a.results.resolve();figures=root/'figures';figures.mkdir(exist_ok=True)
    names={'prism':'Призма','asymmetric':'Асимметричный тетраэдр','absorbing':'Поглощающее тело','concave':'Вогнутая частица'};images=[]
    def save(fig,name):
        path=figures/name;fig.savefig(path,dpi=140,bbox_inches='tight');plt.close(fig);images.append(path)
    cpu=json.loads((root/'cpu/comparison.json').read_text());rows=[]
    for record in cpu:
        rows.append('<tr><td>{}</td><td>{}</td><td>{:.3f}</td><td>{:.3f}</td><td>{:.3f}</td><td>{}/{}</td></tr>'.format(names[record['case']],record['local'],100*record['statistical95_scaled_M11'],100*record['ratio95'],100*record['independent_dense']['residual_plus95_scaled_M11'],record['evaluated_theta'],record['requested_theta']))
    for case in names:
        folder=root/'cpu'/(case+'_auto');data=np.loadtxt(folder/'mueller_fullauto.dat',skiprows=1);direct=np.loadtxt(folder/'mueller_evaluated.dat',skiprows=1);ci=np.loadtxt(folder/'all_mueller_confidence.csv',delimiter=',',skiprows=1);ratio=np.loadtxt(folder/'normalized_mueller_confidence.csv',delimiter=',',skiprows=1)
        fig,axes=plt.subplots(1,3,figsize=(13,3.8),layout='constrained')
        for ax,col,label in zip(axes,[0,1,5],['M11','M12/M11','M22/M11']):
            value=data[:,2] if col==0 else data[:,col+2]/data[:,2];points=direct[:,2] if col==0 else direct[:,col+2]/direct[:,2]
            ax.plot(data[:,0],value,label='Адаптивная оценка');ax.scatter(direct[:,0],points,s=12,label='Прямые θ');ax.set_xlabel('θ, градусы');ax.set_title(label);ax.grid(alpha=.25);ax.legend(fontsize=8)
        fig.suptitle(names[case]+'; реально рассчитанное усреднение, EPS=1%');save(fig,case+'_matrices.png')
        fig,axes=plt.subplots(1,2,figsize=(11,3.7),layout='constrained')
        axes[0].plot(ci[:,0],ci[:,1:].max(1)*100,label='Все 16 элементов / M11');axes[0].plot(ratio[:,0],ratio[:,1:].max(1)*100,label='Все 16 отношений')
        axes[0].axhline(.75,color='k',linestyle='--',label='верхний статистический бюджет');axes[0].set_ylabel('Оценка 95%, %');axes[0].legend(fontsize=8)
        history=np.loadtxt(folder/'refinements.tsv',ndmin=2);axes[1].semilogy(history[:,0],history[:,4]*100,'o-',label='Статистика');axes[1].semilogy(history[:,0],history[:,6]*100,'o-',label='Интерполяция');axes[1].axhline(.25,color='k',linestyle='--');axes[1].set_xlabel('Раунд');axes[1].legend(fontsize=8)
        for ax in axes:ax.grid(alpha=.25)
        axes[0].set_xlabel('θ, градусы');fig.suptitle(names[case]+': точность и ход уточнения');save(fig,case+'_accuracy.png')
    pairs=json.loads((root/'paired_cpu/comparison.json').read_text()) if (root/'paired_cpu/comparison.json').exists() else []
    gpu_path=root/'remote/paired_gpu_comparison.json'
    if gpu_path.exists():pairs+=json.loads(gpu_path.read_text())
    if pairs:
        fig,ax=plt.subplots(figsize=(9,4),layout='constrained');labels=[names[r['case']]+' '+r['backend'].upper() for r in pairs];ax.bar(np.arange(len(pairs)),[r['paired_speedup'] for r in pairs]);ax.set_xticks(np.arange(len(pairs)),labels,rotation=20,ha='right');ax.axhline(1,color='k',linewidth=1);ax.set_ylabel('Ускорение парного расчёта, ×');ax.set_title('Общая трассировка; когерентные поля считаются отдельно');save(fig,'paired_speed.png')
    cal_path=root/'phi_auto/calibration.json'
    if cal_path.exists():
        cal=json.loads(cal_path.read_text());fig,ax=plt.subplots(figsize=(7,4),layout='constrained');ax.plot([r['phi'] for r in cal['phi_candidates']],[r['score'] for r in cal['phi_candidates']],'o-');ax.set_xlabel('φ');ax.set_ylabel('Время × оценка 95%²');ax.set_title('Рабочая сетка; целевая глубина '+str(cal['calibration_depth'])+'; выбран φ='+str(cal['selected_phi']));ax.grid(alpha=.25);save(fig,'phi_selection.png')
    manual=(ROOT/'docs/FULLAUTO_ADAPTIVE_RU.html').read_text();page=manual.rsplit('</html>',1)[0]
    page+='<h2>Проверенные результаты</h2><p>Таблица относится к ориентационной статистике и интерполяции при фиксированной физической модели. Полная сетка содержит 33 угла в 0–25°. Независимый Python-аудит восстановил все матрицы из исходных worker-файлов и проверил каждую запрошенную строку.</p>'
    page+='<table><tr><th>Форма</th><th>Режим N</th><th>Статистика, %</th><th>Отношения, %</th><th>Интерполяция, %</th><th>Прямые θ</th></tr>'+''.join(rows)+'</table>'
    for path in images:page+='<p><img style="max-width:100%" src="data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode()+'" alt="Реально рассчитанные данные"></p>'
    physics=root/'physics/physics_validation.json'
    if physics.exists():
        result=json.loads(physics.read_text());page+='<h2>Отдельная проверка физической глубины</h2><p>Для малой призмы при λ=1 мкм, H=0,2 мкм, D=0,1 мкм, m=1,3116 сравнение n=4 с n=6 дало остаток с оценкой95% <b>{:.3f}%</b>. Поэтому программа вернула physics_validation_failed при EPS=1%, сохранив матрицу для запрошенного n=4. Это показывает разницу между сходимостью усреднения и чувствительностью к физической глубине.</p>'.format(result['depth_plus2_residual95_scaled_M11']*100)
    page+='<p>Исходные файлы: cpu/, remote/cases/, remote/paired/, paired_cpu/, phi_auto/, physics/. Команды и логи сохранены рядом с результатами. Измеренные ускорения парного расчёта не являются ускорением всей кампании.</p></html>'
    destination=root/'Fullauto_адаптация_и_проверки.html';destination.write_text(page);print(destination)
if __name__=='__main__':main()
