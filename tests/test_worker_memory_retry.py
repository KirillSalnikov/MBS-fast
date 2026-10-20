#!/usr/bin/env python3
"""Exercise one-shot resource retry, bounded failure and unchanged Mueller matrices."""
import argparse,json,os,subprocess,shutil
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.output.resolve();r.mkdir(parents=True,exist_ok=False);binary=a.binary.resolve();shutil.copy2(str(Path(__file__).with_name('memory_fault_worker.py')),str(r/'fake_memory_worker.py'));(r/'fake_memory_worker.py').chmod(0o755);base=[str(a.probe.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1','--ri','1.3116','0','--wavelength-um','1','--max-reflections','4','--threads','2','--fullauto','.01','--fullauto-analytic-controls','off','--fullauto-mirror','off','--fullauto-resources','off','--cutoff-profile','off','--phi-points','4','--scattering-grid','170','180','4','4','--fullauto-initial','256','--fullauto-orientation-block','128','--fullauto-depth-sharing','off','--fullauto-workers','persistent','--fullauto-work-queue','on']
report=[];reference=None
for mode in ['none','once','always','physics']:
 out=r/('retry_test_'+mode);env=dict(os.environ,TASK_MBS_WORKER=str(r/'fake_memory_worker.py'),TASK_REAL_WORKER=str(binary),TASK_FAULT=mode,TASK_FAULT_MARKER=str(r/(mode+'.marker')))
 with (r/(mode+'.log')).open('w') as f:p=subprocess.run(base+['--output',str(out)],env=env,stdout=f,stderr=subprocess.STDOUT,timeout=180)
 events=[json.loads(x) for x in (out/'fullauto_history.jsonl').read_text().splitlines()];retries=[x for x in events if x['event']=='worker_memory_retry']
 assert p.returncode==(0 if mode in ['none','once'] else 1),(mode,p.returncode)
 assert len(retries)==(1 if mode in ['once','always'] else 0),(mode,len(retries))
 error=None
 if p.returncode==0:
  matrix=np.loadtxt(out/'frozen_matrix.dat',skiprows=1)
  if reference is None:reference=matrix
  error=float(np.max(abs(matrix[:,2:]-reference[:,2:])/reference[:,2,None]));assert error<1e-12,error
  if mode=='once':assert list((out/'jobs').glob('*/stdout.memory-retry.log'))
 report.append(dict(mode=mode,returncode=p.returncode,retries=len(retries),matrix_error=error))
(r/'retry_test_result.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
