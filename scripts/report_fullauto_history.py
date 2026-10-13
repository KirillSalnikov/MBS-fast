#!/usr/bin/env python3
"""Render saved adaptive decisions and real round matrices as standalone HTML."""
import argparse,base64,html,json
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--output',type=Path);a=p.parse_args();root=a.root.resolve();target=a.output or root/'fullauto_history.html'
 entries=[];invalid=0
 for line in (root/'fullauto_history.jsonl').read_text().splitlines():
  try:entries.append(json.loads(line))
  except ValueError:invalid+=1
 params=json.loads((root/('final_parameters.json' if (root/'final_parameters.json').exists() else 'fullauto_parameters.json')).read_text())
 style='body{font:16px system-ui;max-width:1250px;margin:30px auto;padding:0 20px;color:#233}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd;padding:7px;vertical-align:top}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f3f6f8;padding:12px}img{max-width:100%}details{margin:10px 0}p{line-height:1.5}'
 page=['<!doctype html><html lang="ru"><meta charset="utf-8"><title>Журнал fullauto</title><style>'+style+'</style><h1>Параметры и ход fullauto</h1>',
 '<p>Данные взяты из расчёта. Время работников может перекрываться; его сумма не является временем завершения всей программы. Полная история: fullauto_history.jsonl. Итог: final_parameters.json. Снимки раундов: round_snapshots/.</p>']
 if invalid:page.append('<p>Неполных строк журнала после прерывания: '+str(invalid)+'. Они сохранены в исходном файле.</p>')
 page.append('<h2>Все итоговые параметры</h2><pre>'+html.escape(json.dumps(params,ensure_ascii=False,indent=2))+'</pre>')
 plots=[]
 try:
  import numpy as np
  import matplotlib;matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  decisions=[e['data'] for e in entries if e['event']=='round_decision']
  if decisions:
   fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
   for key,label in [('statistical95','Прямые узлы'),('output_statistical95','После интерполяции')]:axes[0].semilogy([d['round'] for d in decisions],[100*d[key] for d in decisions],'o-',label=label)
   axes[0].axhline(100*params['statistical_budget'],color='k',ls='--',label='Бюджет статистики')
   axes[1].semilogy([d['round'] for d in decisions],[max(1e-12,100*d['interpolation_residual95']) for d in decisions],'o-',label='Контроль интерполяции')
   axes[1].axhline(100*params['interpolation_budget'],color='k',ls='--',label='EPS/4')
   for ax in axes:ax.set_xlabel('Раунд');ax.set_ylabel('% от M11');ax.grid(alpha=.25);ax.legend(fontsize=8)
   plots.append(fig)
  snapshots=sorted((root/'round_snapshots').glob('*/mueller_evaluated.dat'),key=lambda p:int(p.parent.name.split('_')[0][5:]))
  if snapshots:
   fig,axes=plt.subplots(1,3,figsize=(14,4),layout='constrained')
   for path in snapshots:
    m=np.loadtxt(path,skiprows=1,ndmin=2);label=path.parent.name.split('_')[0]
    for ax,col in zip(axes,[2,3,7]):ax.plot(m[:,0],m[:,col] if col==2 else m[:,col]/m[:,2],'.-',label=label,ms=3)
   for ax,title in zip(axes,['M11','M12/M11','M22/M11']):ax.set_title(title);ax.set_xlabel('θ, градусы');ax.grid(alpha=.25);ax.legend(fontsize=7)
   axes[0].set_yscale('log');plots.append(fig)
  for fig in plots:
   import io
   buf=io.BytesIO();fig.savefig(buf,format='png',dpi=140);plt.close(fig);page.append('<img alt="Реальные промежуточные результаты" src="data:image/png;base64,'+base64.b64encode(buf.getvalue()).decode()+'">')
 except ImportError:page.append('<p>Для графиков установите NumPy и Matplotlib. Все значения сохранены ниже.</p>')
 page.append('<h2>Решения в порядке выполнения</h2><p>Для каждого события сохранены время, номер процесса и значения. worker_start ссылается на полную команду; worker_complete — на измеренное время и профили фаз.</p>')
 worker_count=sum(e['event']=='worker_start' for e in entries);page.append('<p>Запусков работников: '+str(worker_count)+'. Их отдельные команды и stdout.log лежат в jobs/.</p>')
 for e in entries:
  if e['event'] in ['worker_start','worker_complete','worker_cache_reuse']:continue
  title='%.2f с · %s · процесс %s'%(e['elapsed_seconds'],e['event'],e['controller_pid'])
  page.append('<details><summary>'+html.escape(title)+'</summary><pre>'+html.escape(json.dumps(e['data'],ensure_ascii=False,indent=2))+'</pre></details>')
 page.append('</html>');target.parent.mkdir(parents=True,exist_ok=True);target.write_text('\n'.join(page));print(target)
if __name__=='__main__':main()
