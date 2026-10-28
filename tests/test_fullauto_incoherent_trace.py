#!/usr/bin/env python3
"""Exercise automatic tracing selection for incoherent PO on a real CUDA job."""
import argparse,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=True)
cmd=[str(a.binary.resolve()),'--method','po','--backend','cuda','--incoherent','--particle','1','5','2.5',
 '--ri','1.3116','0','--wavelength-um','.532','--symmetry','1','1','--max-reflections','6',
 '--fullauto','.04','--fullauto-gpus','0','--threads','2','--fullauto-resources','off',
 '--fullauto-analytic-controls','off','--fullauto-mirror','off','--fullauto-rotation-symmetry','off',
 '--phi-points','4','--scattering-grid','160','170','4','8','--fullauto-initial','1024',
 '--fullauto-pilot','64','--fullauto-kernel-count','64','--fullauto-min-correction','64',
 '--fullauto-theta-start','5','--max-orientations','2048','--fullauto-max-rounds','1',
 '--fullauto-final-audit','off','--fullauto-final-checks','off','--output',str(root/'result'),'--close']
(root/'command.args').write_text('\n'.join(cmd)+'\n')
with (root/'stdout.log').open('w')as log:q=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=360)
assert q.returncode in [0,3,4],q.returncode
r=root/'result';trace=json.loads((r/'trace_backend_calibration.json').read_text())
assert {x['backend']for x in trace['candidates']}=={'cpu','resident','hybrid'},trace
assert any(x.get('gpu_executed',False)for x in trace['candidates']),trace
cpu=trace['candidates'][0]['wall_seconds']
if trace['selected']!='cpu':
 selected=next(x for x in trace['candidates']if x['backend']==trace['selected'])
 assert selected['accepted'] and selected['error']<trace['tolerance'] and selected['wall_seconds']<.95*cpu
assert json.loads((r/'depth_selection.json').read_text())['explicit_target_preserved']
(root/'verification.json').write_text(json.dumps(dict(passed=True,exit_code=q.returncode,trace=trace),indent=2)+'\n')
print(json.dumps(trace))
