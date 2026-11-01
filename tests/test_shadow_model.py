#!/usr/bin/env python3
"""Shadow model selection must survive tracing clones and CPU/CUDA preparation."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('binary', type=Path)
p.add_argument('--backend', choices=['cpu', 'cuda'], default='cpu')
p.add_argument('--trace-backend', choices=['cpu', 'resident'], default='cpu')
p.add_argument('--out', type=Path)
a = p.parse_args()
temporary = tempfile.TemporaryDirectory(prefix='mbs-shadow-model-') if a.out is None else None
root = Path(temporary.name) if temporary else a.out.resolve()
root.mkdir(parents=True, exist_ok=True)
base = [str(a.binary.resolve()), '--method', 'po', '--backend', a.backend,
        '--particle', '1', '1.2371791482634837', '1', '--k-eq', '40',
        '--wavelength-um', '.532', '--max-reflections', '4', '--threads', '2',
        '--sobol-seed', '16', '11', '--haar-alpha', '--scattering-grid', '0', '180', '12', '36',
        '--cutoff-profile', 'off', '--close']
if a.backend == 'cuda':
    base += ['--po-trace-backend', a.trace_backend, '--orientation-chunk', '8']


def run(name, model=None, index='1.3116', extra=()):
    folder = root / name
    assert not folder.exists(), folder
    cmd = base + ['--refractive-index', index, '0', '--output', str(folder)] + list(extra)
    if model:
        cmd += ['--po-shadow-model', model]
    (root / (name + '_command.json')).write_text(json.dumps(cmd, indent=2) + '\n')
    with (root / (name + '.log')).open('w') as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT,
                       env={**os.environ, 'OMP_NUM_THREADS': '2', 'MBS_GPU_MEM_FRACTION': '.4'})
    x = np.loadtxt(folder / (name + '.dat'), skiprows=1)[:, -16:]
    assert np.isfinite(x).all()
    return x


def difference(x, y):
    return float(np.max(abs(x-y) / np.maximum(abs(y[:, :1]), 1)))


default = run('default')
projected = run('projected', 'projected')
surface = run('surface', 'surface')
no_projected = run('no_projected', 'projected', extra=['--no-shadow-beam'])
no_surface = run('no_surface', 'surface', extra=['--no-shadow-beam'])
vacuum_projected = run('vacuum_projected', 'projected', index='1')
vacuum_surface = run('vacuum_surface', 'surface', index='1')
tolerance = 2e-7 if a.backend == 'cuda' else 1e-10
report = {
    'default_vs_projected': difference(default, projected),
    'surface_vs_projected': difference(surface, projected),
    'unchanged_particle_beams': difference(no_projected, no_surface),
    'vacuum_projected_max_abs': float(np.max(abs(vacuum_projected))),
    'vacuum_surface_max_abs': float(np.max(abs(vacuum_surface))),
}
assert report['default_vs_projected'] < tolerance, report
assert report['surface_vs_projected'] > 1e-4, report
assert report['unchanged_particle_beams'] < tolerance, report
assert report['vacuum_surface_max_abs'] < tolerance, report
# This is an explicit compatibility tradeoff, not a transparent-limit fix.
assert report['vacuum_projected_max_abs'] > 1e-4, report
(root / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
print('Shadow selection, unchanged particle beams and explicit vacuum tradeoff passed')
