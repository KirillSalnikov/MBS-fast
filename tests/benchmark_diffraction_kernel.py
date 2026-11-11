#!/usr/bin/env python3
"""Fixed optical work, baseline versus optimized executable; no adaptive stopping."""
import argparse,json,time,subprocess,re,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--candidate',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--orientations',default='256');p.add_argument('--length',type=float,default=20);p.add_argument('--theta-min',default='0');p.add_argument('--theta-max',default='180');p.add_argument('--repeats',type=int,default=3);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);results={}
for case,absorption,incoherent in [('coherent','0',False),('coherent_absorbing','.00042407',False),('incoherent','0',True)]:
 values={k:[] for k in ['baseline','candidate']};matrices={}
 for repeat in range(a.repeats):
  for label in (['baseline','candidate'] if repeat%2==0 else ['candidate','baseline']):
   dest=out/(case+'_'+label+str(repeat));cmd=[str(getattr(a,label).resolve()),'--method','po','--backend','cpu','--particle','2',str(a.length),str(a.length/2),'--ri','1.3116',absorption,'--wavelength-um','.532','--max-reflections','6','--sobol-seed',a.orientations,'11','--haar-alpha','--scattering-grid',a.theta_min,a.theta_max,'8','256','--phi-points','8','--symmetry','1','1','--threads','1','--cutoff-profile','off','--output',str(dest),'--close']
   if incoherent:cmd+=['--incoherent']
   (out/(dest.name+'.argv.json')).write_text(json.dumps(cmd));start=time.monotonic()
   with (out/(dest.name+'.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=120)
   elapsed=time.monotonic()-start;log=(out/(dest.name+'.log')).read_text();phase=float(re.search(r'Phase 2 \(diffraction, OpenMP\): ([\d.]+) s',log)[1]);values[label].append({'wall_seconds':elapsed,'diffraction_seconds':phase})
   rows=[[float(x) for x in line.split()] for line in (dest/(dest.name+'.dat')).read_text().splitlines()[1:] if line.strip()];matrices[label]=rows
 x,y=matrices['baseline'],matrices['candidate'];assert len(x)==len(y)
 diff=max(abs(u[j]-v[j])/max(abs(u[2]),1e-300) for u,v in zip(x,y) for j in range(2,18));assert diff<1e-8,(case,diff)
 med={k:statistics.median(v['wall_seconds'] for v in z) for k,z in values.items()};results[case]={'runs':values,'median_wall_seconds':med,'speedup':med['baseline']/med['candidate'],'max_all16_difference_scaled_M11':diff};(out/'receipt.json').write_text(json.dumps(results,indent=2));print(case,results[case],flush=True)
