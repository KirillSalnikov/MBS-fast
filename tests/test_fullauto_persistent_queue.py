#!/usr/bin/env python3
import argparse,pathlib,subprocess,os,json,time
import numpy as np
parser=argparse.ArgumentParser(description='Compare persistent workers, durable resume and SIGKILL recovery with real frozen fullauto blocks.')
parser.add_argument('--probe',type=pathlib.Path,required=True)
parser.add_argument('--binary',type=pathlib.Path,required=True)
parser.add_argument('--output',type=pathlib.Path,required=True)
settings=parser.parse_args();p=settings.output.resolve();p.mkdir(parents=True)
binary=str(settings.binary.resolve())
base=[str(settings.probe.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1','--ri','1.3116','0','--wavelength-um','1','--max-reflections','4','--threads','1','--fullauto','.02','--fullauto-analytic-controls','off','--fullauto-mirror','off','--fullauto-resources','off','--cutoff-profile','off','--phi-points','2','--scattering-grid','0','180','2','8','--fullauto-initial','256','--fullauto-orientation-block','128','--fullauto-depth-sharing','off','--fullauto-work-queue','on','--fullauto-cpu-workers','2']
report=[]
for mode in ['process','persistent']:
 dest=p/('queue4_'+mode);args=base+['--fullauto-workers',mode,'--output',str(dest)]
 (p/('queue4_'+mode+'.args')).write_text('\n'.join(args)+'\n');started=time.time()
 with (p/('queue4_'+mode+'.log')).open('w') as log:subprocess.run(args,env=dict(os.environ,TASK_MBS_WORKER=binary),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
 events=[json.loads(l) for l in (dest/'fullauto_history.jsonl').read_text().splitlines()];starts=[r for r in events if r['event']=='worker_start']
 report.append({'mode':mode,'seconds':time.time()-started,'jobs':len(starts),'pids':len({r['data']['pid'] for r in starts})})
a=np.loadtxt(p/'queue4_process/frozen_matrix.dat',skiprows=1);b=np.loadtxt(p/'queue4_persistent/frozen_matrix.dat',skiprows=1)
error=float(np.max(abs(a[:,2:]-b[:,2:])/a[:,2,None]));assert error<1e-12
assert report[1]['pids']==2 and report[0]['pids']==report[0]['jobs']
(p/'queue4_validation.json').write_text(json.dumps({'modes':report,'max_all16_scaled_M11':error},indent=2)+'\n');print(report,error)
# Disk cache remains sufficient after all workers exit: no requests on resume.
dest=p/'queue4_persistent';args=base+['--fullauto-workers','persistent','--output',str(dest)]
before=(dest/'fullauto_history.jsonl').read_text().count('"event":"worker_start"')
with (p/'queue4_resume.log').open('w') as log:subprocess.run(args,env=dict(os.environ,TASK_MBS_WORKER=binary,TASK_MBS_WARM_CACHE='1'),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=60)
assert before==(dest/'fullauto_history.jsonl').read_text().count('"event":"worker_start"')
print('Resume: all saved blocks reused, zero workers launched')

import pathlib,subprocess,os,json,time,signal,numpy as np
P=p
args=(P/'queue4_persistent.args').read_text().splitlines();root=P/'queue4_killed';args[args.index('--output')+1]=str(root)
env=dict(os.environ,TASK_MBS_WORKER=binary)
with (P/'queue4_killed.log').open('w') as log:
 worker=subprocess.Popen(args,env=env,stdout=log,stderr=subprocess.STDOUT)
 deadline=time.time()+30;events=[]
 while time.time()<deadline:
  if (root/'fullauto_history.jsonl').exists():
   events=[]
   for line in (root/'fullauto_history.jsonl').read_text().splitlines():
    try:events.append(json.loads(line))
    except ValueError:pass
   if sum(r['event']=='worker_complete' for r in events)>=8:break
  assert worker.poll() is None,'controller finished before interruption'
  time.sleep(.02)
 else:raise RuntimeError('No completed blocks')
 pids={r['data']['pid'] for r in events if r['event']=='worker_start'}
 saved={str(f.parent):f.read_bytes() for f in (root/'jobs').glob('*/complete.txt')}
 assert saved
 worker.kill();worker.wait()
 time.sleep(.4)
 def alive(pid):
  try:return pathlib.Path('/proc/%d/stat'%pid).read_text().split(') ',1)[1][0]!='Z'
  except OSError:return False
 assert all(not alive(pid) for pid in pids),'orphan worker remained active'
 with (P/'queue4_killed_resume.log').open('w') as log:subprocess.run(args,env=dict(env,TASK_MBS_WARM_CACHE='1'),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
 for f,contents in saved.items():assert pathlib.Path(f,'complete.txt').read_bytes()==contents
 x=np.loadtxt(P/'queue4_process/frozen_matrix.dat',skiprows=1);y=np.loadtxt(root/'frozen_matrix.dat',skiprows=1)
 error=float(np.max(abs(x[:,2:]-y[:,2:])/x[:,2,None]));assert error<1e-12
 result={'saved_blocks_preserved':len(saved),'worker_pids_terminated':list(pids),'max_all16_scaled_M11':error,'parent_sigkill_recovery':True}
 (P/'recovery4_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
