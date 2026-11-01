#!/usr/bin/env python3
"""Physical transparent-body limit through native PO paths (not a zero-index shortcut)."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('binary', type=Path)
p.add_argument('--backend', choices=['cpu','cuda'], default='cpu')
p.add_argument('--trace-backend', choices=['cpu','resident'], default='cpu')
p.add_argument('--out', type=Path)
a = p.parse_args()
binary = a.binary.resolve()
temporary = tempfile.TemporaryDirectory(prefix='mbs-transparent-') if a.out is None else None
root = Path(temporary.name) if temporary else a.out.resolve()
root.mkdir(parents=True, exist_ok=True)
base = [str(binary), '--method','po','--backend',a.backend,'--particle','1','1.2371791482634837','1',
        '--k-eq','40','--wavelength-um','.532','--max-reflections','14','--threads','2',
        '--scattering-grid','0','180','12','36','--cutoff-profile','off','--close']
results = {}
def run(name, orientation, index='1', extra=(), particle=None):
    folder = root/name
    assert not folder.exists(), f'Use an empty output directory: {folder}'
    if a.backend == 'cuda' and '--fixed-orientation' in orientation:
        rule = root/(name+'_orientation.txt')
        rule.write_text(f'{orientation[1]} {orientation[2]} 0 1\n')
        orientation = ['--so3-rule',str(rule)]
    args = base + ['--refractive-index',index,'0'] + orientation + list(extra) + ['--output',str(folder)]
    if a.backend == 'cuda':
        args += ['--po-trace-backend',a.trace_backend,'--orientation-chunk','8']
        if a.trace_backend == 'resident':
            # Two interactions suffice at zero contrast and keep zero-amplitude
            # branches below the resident queue limit with cutoffs disabled.
            args[args.index('--max-reflections')+1]='2'
    if particle is not None:
        start=args.index('--particle')+1;args[start:start+3]=particle
    (root/(name+'_command.json')).write_text(json.dumps(args,indent=2)+'\n')
    with (root/(name+'.log')).open('w') as log:
        subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=True,
                       env={**os.environ,'OMP_NUM_THREADS':'2','MBS_GPU_MEM_FRACTION':'.4'})
    x = np.loadtxt(folder/(name+'.dat'), skiprows=1)[:,-16:]
    assert np.all(np.isfinite(x)), name
    value = float(np.max(np.abs(x)))
    results[name] = value
    print(name, 'max |Mij| =', value, flush=True)
    return value

threshold = 1e-8 if a.backend == 'cpu' else 2e-7
sobol = ['--sobol-seed','16','11','--haar-alpha']
for name, orientation, extra in [
    ('aligned', ['--fixed-orientation','0','0'], []),
    ('tilted', ['--fixed-orientation','37','19'], []),
    ('averaged', sobol, []),
    ('forced_nonconvex_trace', sobol, ['--geometry','nonconvex']),
]:
    assert run(name,orientation,extra=extra) < threshold, name
assert run('concave',sobol,particle=['10','1.2371791482634837','1','.3']) < threshold
# Removing the incident-field contribution must restore the transmitted beam.
assert run('without_shadow',sobol,extra=['--no-shadow-beam']) > 1.
# Keep the contrast signal above FP32 cancellation noise on CUDA.
contrasts = ('1.00000001','1.0000000001') if a.backend == 'cpu' else ('1.0001','1.000001')
near = run('small_contrast', ['--fixed-orientation','37','19'], index=contrasts[0])
nearer = run('smaller_contrast', ['--fixed-orientation','37','19'], index=contrasts[1])
assert 0 < nearer < near, (near,nearer)
assert nearer < 1e-3*near + threshold
(root/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print('Transparent PO and vanishing-contrast checks passed')
