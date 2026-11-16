from pathlib import Path
import os,json,time,subprocess,math,traceback,hashlib
import argparse
parser=argparse.ArgumentParser(description='ABA all-element CUDA validation of analytic mean geometry reuse')
parser.add_argument('--reference',required=True)
parser.add_argument('--candidate',required=True)
parser.add_argument('--particle-file',required=True)
parser.add_argument('--output',required=True)
args=parser.parse_args()
p=Path(args.output).resolve();p.mkdir(parents=True,exist_ok=False)
old=str(Path(args.reference).resolve());new=str(Path(args.candidate).resolve());shape=str(Path(args.particle_file).resolve())
env=os.environ.copy();env.update(OMP_NUM_THREADS='2',MBS_HOST_MEM_BUDGET_MB='512',MBS_GPU_MEM_FRACTION='.05')
s={'pid':os.getpid(),'state':'running','checks':[],'timing_valid':False,'reason':'caller must isolate GPU and CPU for timing claims','hashes':{k:hashlib.sha256(Path(v).read_bytes()).hexdigest() for k,v in [('old',old),('new',new)]}}
def save():
 q=p/'native_status.tmp';q.write_text(json.dumps(s,indent=2));q.replace(p/'native_status.json')
def matrix(f):
 rows=[]
 for line in f.read_text().splitlines():
  try:r=list(map(float,line.split()))
  except ValueError:continue
  if len(r)==18:rows.append(r)
 assert rows and all(all(math.isfinite(v) for v in r) for r in rows)
 return rows
save()
try:
 cases=[('arb',['--particle-file',shape,'--resize-dmax-um','15.85','--symmetry','1','1'],'circular','.355','1.3243','0'),('bullet',['--particle','2','10','9.85','--symmetry','1','6'],'facets','.4','1.3194','0'),('droxtal',['--particle','4','5','--symmetry','2','6'],'circular','1.55','1.2906','.00042407')]
 theta=p/'native_theta.txt';theta.write_text(''.join(str(i*180/64)+'\n' for i in range(65)))
 for name,args,shadow,w,n,k in cases:
  for incoherent in [False,True]:
   label=name+('_incoherent' if incoherent else '_coherent');data=[];durations=[]
   for tag,exe in [('old',old),('new',new),('old_repeat',old)]:
    out=p/'native_validation_v3'/(label+'_'+tag);out.parent.mkdir(exist_ok=True)
    cmd=[exe,'--method','po','--backend','cuda','--sobol-seed','32','53','--haar-alpha','--phi-points','8','--max-reflections','6','--threads','2','--orientation-chunk','16','--theta-grid-file',str(theta),'--cutoff-profile','balanced','--po-trace-backend','cpu','--orientation-pipeline','--analytic-facet-average','--analytic-shadow-control',shadow,'--analytic-mean-threads','2','--wavelength-um',w,'--ri',n,k,'--close','--output',str(out)]+args+(['--incoherent'] if incoherent else [])
    (out.parent/(out.name+'_cmd.json')).write_text(json.dumps(cmd))
    started=time.monotonic()
    with (out.parent/(out.name+'.log')).open('w') as log:subprocess.run(cmd,env=env,stdout=log,stderr=log,check=True,timeout=180)
    durations.append(time.monotonic()-started)
    data.append(matrix(out/(out.name+'.dat')))
   a,b,repeat=data; repeat_error=max(abs(x-y)/max(abs(r[2]),1e-30) for r,t in zip(a,repeat) for x,y in zip(r[2:],t[2:]));assert repeat_error<3e-5; assert len(a)==len(b) and all(x[:2]==y[:2] for x,y in zip(a,b))
   err=max(abs(x-y)/max(abs(r[2]),1e-30) for r,t in zip(a,b) for x,y in zip(r[2:],t[2:]));rat=max(abs(x/r[2]-y/t[2]) for r,t in zip(a,b) for x,y in zip(r[2:],t[2:]))
   s['checks'].append({'case':label,'all16_scaled_M11':err,'ratio_difference':rat,'identical':a==b,'seconds_ABA':durations,'baseline_repeat_error':repeat_error});save();assert max(err,rat)<3e-5
 s['state']='passed'
except Exception as e:s.update(state='failed',error=str(e));traceback.print_exc()
finally:save()
if s['state']!='passed':raise SystemExit(1)
