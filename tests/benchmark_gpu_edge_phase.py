import os,json,subprocess,time,statistics
from pathlib import Path
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--baseline',type=Path,required=True)
parser.add_argument('--candidate',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
old=args.baseline.resolve();new=args.candidate.resolve()
out=args.output.resolve();out.mkdir(parents=True,exist_ok=False)
env=os.environ.copy();env['MBS_GPU_TIMING']='1'
cases=[('small',5,0,180,'0',False),('whole',20,0,180,'0',False),('forward1000',1000,0,1,'0',False),('back1000',1000,170,180,'0',False),('abs_back1000',1000,170,180,'.00042407',False),('incoherent',20,150,170,'0',True)]
receipt={}
for name,L,lo,hi,imag,incoh in cases:
 record={'runs':{},'commands':{}};matrices={}
 for repeat in range(2):
  for label in (['baseline','candidate'] if repeat%2==0 else ['candidate','baseline']):
   dest=out/(name+'_'+label+str(repeat));binary=old if label=='baseline' else new
   cmd=[str(binary),'--method','po','--backend','gpu','--particle','2',str(L),str(L/2),'--ri','1.3116',imag,'--wavelength-um','.532','--max-reflections','6','--sobol-seed','64','11','--haar-alpha','--scattering-grid',str(lo),str(hi),'8','256','--phi-points','8','--symmetry','1','1','--threads','2','--cutoff-profile','off','--output',str(dest),'--close']
   if incoh:cmd+=['--incoherent']
   record['commands'][label]=cmd
   start=time.monotonic()
   with (out/(dest.name+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,timeout=180)
   elapsed=time.monotonic()-start
   if r.returncode:raise RuntimeError((name,label,r.returncode))
   record['runs'].setdefault(label,[]).append(elapsed)
   matrices[label]=[[float(x) for x in line.split()] for line in (dest/(dest.name+'.dat')).read_text().splitlines()[1:] if line.strip()]
 a,b=matrices['baseline'],matrices['candidate'];assert len(a)==len(b)
 scale=lambda row:max(abs(row[2]),1e-300)
 difference=max(abs(x[j]-y[j])/scale(x) for x,y in zip(a,b) for j in range(2,18))
 ratios=max(abs(x[j]/scale(x)-y[j]/scale(y)) for x,y in zip(a,b) for j in range(3,18))
 record.update(max_all16_scaled_M11=difference,max_ratio_difference=ratios,pass_accuracy=difference<=.001 and ratios<=.001,median_speedup=statistics.median(record['runs']['baseline'])/statistics.median(record['runs']['candidate']))
 receipt[name]=record;(out/'receipt.json').write_text(json.dumps(receipt,indent=2));print(name,difference,ratios,record['median_speedup'],flush=True)
if not all(r['pass_accuracy'] for r in receipt.values()):raise RuntimeError('Accuracy regression: see receipt')
