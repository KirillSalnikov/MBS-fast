#!/usr/bin/env python3
"""Compare real paired pilots with a baseline and check common-depth optics.

The reference policy may select a different depth than the adjacent-pair
baseline. Compare unchanged optical matrices only at identical depth/sample
prefixes; never require physically different truncations to be identical.
"""
import argparse,json,os,signal,subprocess,time
import numpy as np
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True)
p.add_argument('--baseline',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=True);records=[]
for name,binary in [('baseline',a.baseline),('candidate',a.binary)]:
 out=root/name
 cmd=[str(binary.resolve()),'--method','po','--backend','cpu','--particle','1','5','2.5',
      '--ri','1.3116','0','--wavelength-um','.532','--threads','1','--fullauto','.04',
      '--fullauto-depth-start','6','--fullauto-depth-max','12','--fullauto-depth-pilot','8',
      '--fullauto-initial','32','--fullauto-min-correction','8','--fullauto-pilot','8',
      '--fullauto-max-rounds','1','--max-orientations','128','--phi-points','1',
      '--scattering-grid','170','180','1','2','--fullauto-theta-start','3',
      '--fullauto-resources','off','--fullauto-cutoffs','off','--cutoff-profile','off',
      '--fullauto-analytic-controls','off','--fullauto-mirror','off','--fullauto-rotation-symmetry','off',
      '--fullauto-final-checks','off','--symmetry','1','1','--output',str(out),'--close']
 (root/(name+'.args')).write_text('\n'.join(cmd)+'\n');start=time.monotonic()
 with (root/(name+'.log')).open('w')as log:
  r=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  try:r.wait(timeout=300)
  except subprocess.TimeoutExpired:
   os.killpg(r.pid,signal.SIGTERM)
   try:r.wait(timeout=5)
   except subprocess.TimeoutExpired:os.killpg(r.pid,signal.SIGKILL);r.wait()
   raise
 assert r.returncode in [0,3,4],(name,r.returncode)
 pilots=[json.loads(f.read_text())for f in out.glob('depth_selection_candidate*.json')]
 pilots.sort(key=lambda x:(x['candidate'],x['count_per_seed']))
 assert pilots
 selection=json.loads((out/'depth_selection.json').read_text())
 if selection['pilot_passed']:assert selection['residual_upper95']<=selection['budget']
 else:assert not (out/'adaptive_checkpoint.tsv').exists()
 if selection.get('policy')=='common-reference-local-pilot-v2':
  assert selection['comparison_depth']==16 and selection['reference_depth']==14
  if selection['pilot_passed']:
   assert selection['reference_upper95']<=selection['reference_budget']
   assert (out/'automatic_depth_policy.tsv').exists()
 counts={}
 for v in pilots:counts.setdefault(str(v['candidate']),[]).append(v['count_per_seed'])
 if name=='candidate':
  assert len(counts)>1,counts
  assert all(v[0]==8 for v in counts.values()),counts
  assert any(v[-1]>8 for v in counts.values()),'pilot did not exercise an expanded sample'
 history=[json.loads(line)for line in (out/'fullauto_history.jsonl').read_text().splitlines()]
 selected_at=next(row['elapsed_seconds']for row in history if row['event']=='physical_depth_selected')
 records.append(dict(mode=name,seconds=time.monotonic()-start,depth_selection_seconds=selected_at,counts=counts,selection=selection,
     depth_worker_seconds=sum(row['data'].get('seconds',0)for row in history if row['event']=='worker_complete'and row['elapsed_seconds']<=selected_at)))
new_policy=records[1]['selection'].get('policy')=='common-reference-local-pilot-v2'
if new_policy:
 matched=[]
 for old in (root/'baseline/depth_pilot_snapshots').glob('n*_N8/mueller_evaluated.dat'):
  depth=old.parent.name.split('_')[0]
  new=root/'candidate/depth_pilot_snapshots'/f'{depth}_reference16_N8/mueller_evaluated.dat'
  if not new.exists():continue
  a0=np.loadtxt(old,skiprows=1);a1=np.loadtxt(new,skiprows=1)
  assert a0.shape==a1.shape and np.allclose(a0[:,:2],a1[:,:2],rtol=0,atol=1e-12)
  error=float(np.max(np.abs(a0[:,2:]-a1[:,2:])/np.maximum(np.abs(a0[:,2:3]),1e-300)))
  assert error<1e-9,(depth,error)
  matched.append(dict(depth=depth,count=8,error_scaled_M11=error))
 assert matched,'No identical depth/orientation prefixes were compared'
 (root/'common_depth_optics.json').write_text(json.dumps(matched,indent=2)+'\n')
else:
 assert records[0]['selection']['target_depth']==records[1]['selection']['target_depth']
 matrices=[Path(json.loads((root/name/'intermediate_current.json').read_text())['matrix']).read_bytes()for name in ['baseline','candidate']]
 assert matrices[0]==matrices[1], 'production matrices changed on the fixed test sample'
(root/'comparison.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records,indent=2))
