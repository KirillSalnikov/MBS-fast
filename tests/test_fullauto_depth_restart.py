#!/usr/bin/env python3
"""Compare real paired depth pilots with/without carrying an enlarged sample.

Also validates that an accepted depth still meets the same confidence bound.
Run with a baseline binary built from the previous FullAutoPhysics.inc.
"""
import argparse,json,subprocess,time
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
 with (root/(name+'.log')).open('w')as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=300)
 assert r.returncode in [0,3,4],(name,r.returncode)
 pilots=[json.loads(f.read_text())for f in out.glob('depth_selection_candidate*.json')]
 pilots.sort(key=lambda x:(x['candidate'],x['count_per_seed']))
 assert pilots
 selection=json.loads((out/'depth_selection.json').read_text())
 if selection['pilot_passed']:assert selection['residual_upper95']<=selection['budget']
 else:assert not (out/'adaptive_checkpoint.tsv').exists()
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
assert records[0]['selection']['target_depth']==records[1]['selection']['target_depth']
matrices=[Path(json.loads((root/name/'intermediate_current.json').read_text())['matrix']).read_bytes()for name in ['baseline','candidate']]
assert matrices[0]==matrices[1], 'production matrices changed on the fixed test sample'
(root/'comparison.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records,indent=2))
