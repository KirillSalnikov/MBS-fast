#!/usr/bin/env python3
"""A real noisy depth pilot must not silently enter production, even with audits off."""
import argparse,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
r=a.output.resolve();r.mkdir(parents=True,exist_ok=True)
base=[str(a.binary.resolve()),'--method','po','--backend','cpu','--particle','1','20','10','--ri','1.3116','0','--wavelength-um','1','--threads','1','--fullauto','.02','--fullauto-final-checks','off','--fullauto-depth-start','6','--fullauto-depth-max','6','--fullauto-depth-pilot','8','--fullauto-initial','32','--fullauto-min-correction','8','--fullauto-pilot','8','--fullauto-max-rounds','1','--max-orientations','64','--phi-points','1','--scattering-grid','170','180','1','2','--fullauto-theta-start','3','--fullauto-resources','off','--fullauto-cutoffs','off','--cutoff-profile','off','--fullauto-analytic-controls','off','--fullauto-mirror','off','--fullauto-rotation-symmetry','off','--symmetry','1','1','--output',str(r/'result'),'--close']
(r/'command.args').write_text('\n'.join(base)+'\n')
for attempt in [0,1,2]:
 if attempt==2:(r/'result/depth_selection.json').rename(r/'result/depth_selection.before_missing_test.json')
 with (r/('run%d.log'%attempt)).open('w') as log:q=subprocess.run(base,stdout=log,stderr=subprocess.STDOUT,timeout=180)
 assert q.returncode==4,(q.returncode,(r/('run%d.log'%attempt)).read_text()[-2000:])
 status=json.loads((r/'result/fullauto_status.json').read_text());assert status['status']=='depth_selection_unresolved'
 assert not (r/'result/adaptive_checkpoint.tsv').exists()
 selection=json.loads((r/'result/depth_selection.json').read_text());assert not selection['pilot_passed']
 current={str(f):f.stat().st_mtime_ns for f in (r/'result/jobs').glob('*/complete.txt')}
 if attempt:assert current==markers,'Unresolved depth replay repeated completed worker data'
 markers=current
(r/'verification.json').write_text(json.dumps(dict(production_blocked=True,final_checks_off=True,completed_jobs_preserved=len(markers),selection=selection),indent=2)+'\n')
print('Unresolved depth checkpoint and replay PASS')
