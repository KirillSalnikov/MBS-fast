#!/usr/bin/env python3
"""Regression for phase cancellation on an absorbing, triangulated dust shape."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--reference', default=str(ROOT/'gpu/bin/mbs_po_gpu_float'))
parser.add_argument('--candidate', default=str(ROOT/'gpu/bin/mbs_po_gpu_float_reduced'))
args = parser.parse_args()
env = {k:v for k,v in os.environ.items() if not k.startswith(('MBS_', 'OMP_', 'GOMP_'))}
env.update(MBS_GPU_BLOCK='256', MBS_GPU_WARP_BEAMS='0')
with tempfile.TemporaryDirectory(prefix='mbs-reduced-phase-') as tmp:
    work = Path(tmp)
    grid = work/'theta.csv'
    grid.write_text('0\n0.01\n19.1\n86\n101\n170\n171.675\n177.825\n177.875\n180\n')
    base = ['--method','po','--backend','cuda','--pf',str(ROOT/'tests/data/phase_dust15.dat'),
            '--rs','500','--ri','1.2906','.00042407','--wavelength-um','1.55','--absorption',
            '--max-reflections','8','--sobol-seed','64','71','--sobol-offset','1048576',
            '--haar-alpha','--phi-points','8','--theta-grid-file',str(grid),
            '--symmetry','1','1','--cutoff-profile','off','--threads','2','--close']
    for mode in ['coherent','incoherent']:
        matrices=[]
        for label,binary in [('reference',args.reference),('candidate',args.candidate)]:
            dest=work/(mode+'_'+label)
            cmd=[str(Path(binary).resolve())]+base+([] if mode=='coherent' else ['--incoh'])
            cmd+=['--output',str(dest)]
            run=subprocess.run(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,timeout=120)
            if run.returncode:raise RuntimeError(run.stdout)
            matrices.append([[float(x) for x in row.split()] for row in
                             (dest/(dest.name+'.dat')).read_text().splitlines()[1:] if row.strip()])
        a,b=matrices
        assert len(a)==len(b) and len(a)>0
        assert all(len(r)==len(c)==18 and r[:2]==c[:2] and
                   all(math.isfinite(x) for x in c) and r[2]>0 and c[2]>0
                   for r,c in zip(a,b))
        error=max(abs(r[k]-c[k])/r[2] for r,c in zip(a,b) for k in range(2,18))
        ratio=max(abs(r[k]/r[2]-c[k]/c[2]) for r,c in zip(a,b) for k in range(2,18))
        print(json.dumps({'mode':mode,'all16_scaled_M11':error,'ratios':ratio}),flush=True)
        assert max(error,ratio)<3e-5,(mode,error,ratio)
