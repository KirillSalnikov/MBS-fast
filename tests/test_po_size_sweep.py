#!/usr/bin/env python3
"""Check shared PO size scaling against independent traces, including m=1."""
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
p.add_argument('--shadow-model', choices=['projected', 'surface'], default='surface')
p.add_argument('--out', type=Path)
a = p.parse_args()
temporary = tempfile.TemporaryDirectory(prefix='mbs-po-sizes-') if a.out is None else None
root = Path(temporary.name) if temporary else a.out.resolve()
root.mkdir(parents=True, exist_ok=True)
sizes = root/'sizes.txt'
sizes.write_text('20\n40\n')
base = [str(a.binary.resolve()), '--method', 'po', '--backend', a.backend,
        '--po-shadow-model', a.shadow_model,
        '--particle', '1', '1.2371791482634837', '1', '--wavelength-um', '.532',
        '--max-reflections', '4', '--euler-grid', '3', '3',
        '--scattering-grid', '0', '180', '12', '18', '--threads', '2',
        '--cutoff-profile', 'off', '--close']
results = {}


def run(name, index, flags, extra_env=None, geometry=()):
    folder = root/name
    # The executable renames existing output directories: never read old data.
    assert not folder.exists(), f'Use an empty output directory: {folder}'
    args = base + ['--refractive-index', *index] + list(geometry) + flags
    if extra_env:
        args += ['--allow-experimental-environment']
    args += ['--output', str(folder)]
    env = {**os.environ, 'OMP_NUM_THREADS': '2', 'MBS_GPU_MEM_FRACTION': '.4',
           **(extra_env or {})}
    (root/(name+'_command.json')).write_text(json.dumps(args, indent=2)+'\n')
    with (root/(name+'.log')).open('w') as log:
        subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=True, env=env)
    return folder


threshold = 1e-8 if a.backend == 'cpu' else 2e-7
for case, index, geometry in [
    ('transparent', ('1', '0'), ()),
    ('transparent_nonconvex', ('1', '0'), ('--geometry', 'nonconvex')),
    ('lossless', ('1.3116', '0'), ()),
    ('absorbing', ('1.3116', '.02'), ()),
]:
    single = {}
    for k in (20, 40):
        out = run(f'{case}_single{k}', index, ['--k-eq', str(k)], geometry=geometry)
        single[k] = np.loadtxt(out/(out.name+'.dat'), skiprows=1)
    variants = [('default', {})]
    if a.backend == 'cuda':
        # Exercise scaling on the CPU and the fused device path as well as
        # the default GPU pack-time scaling. These reconstruct phases separately.
        variants += [('explicit', {'MBS_SHARED_GPU_EXPLICIT_SCALE': '1'}),
                     ('fused', {'MBS_GPU_MULTI_K_FULL': '1'})]
    for variant, env in variants:
        out = run(f'{case}_{variant}', index, ['--k-eq-list', str(sizes)],
                  extra_env=env, geometry=geometry)
        entry = {}
        for k in (20, 40):
            x = np.loadtxt(out/(out.name+f'_keq{k}.dat'), skiprows=1)
            ref = single[k]
            assert np.all(np.isfinite(x))
            np.testing.assert_allclose(x[:, :2], ref[:, :2], rtol=0, atol=1e-12)
            error = float(np.max(abs(x[:, 2:]-ref[:, 2:])/np.maximum(ref[:, 2, None], 1.)))
            assert error < (1e-8 if a.backend == 'cpu' else 5e-5), (case, variant, k, error)
            maximum = float(np.max(abs(x[:, 2:])))
            if index == ('1', '0') and a.shadow_model == 'surface':
                assert maximum < threshold, (case, variant, k, maximum)
            entry[str(k)] = dict(max_difference_over_max_M11_1=error,
                                 maximum_abs_mueller=maximum)
        results[case+'_'+variant] = entry
        print(case, variant, entry, flush=True)
(root/'results.json').write_text(json.dumps(results, indent=2)+'\n')
print('Shared PO size scaling checks passed')
