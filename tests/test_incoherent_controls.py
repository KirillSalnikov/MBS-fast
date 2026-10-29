#!/usr/bin/env python3
"""Physical incoherent target, facet CV identity, raw diagnostics and prefixes."""
import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

repo = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--backend', choices=['cpu', 'cuda'], default='cpu')
p.add_argument('--tolerance', type=float, default=2e-10)
a = p.parse_args()
root = a.output.resolve()
root.mkdir(parents=True, exist_ok=True)
theta = np.array([0., 17., 90., 170., 180.])
np.savetxt(root/'theta.txt', theta, fmt='%.17g')
weights = np.tile([.75, -.25], (len(theta), 1))
np.savetxt(root/'weights.tsv', np.column_stack([theta, weights]), fmt='%.17g',
           header='theta_deg reflection_weight shadow_weight', comments='')
base = [str(a.binary.resolve()), '--method', 'po', '--backend', a.backend,
        '--wavelength-um', '.532', '--max-reflections', '8', '--cutoff-profile', 'balanced',
        '--sobol-seed', '128', '19', '--haar-alpha', '--theta-grid-file', str(root/'theta.txt'),
        '--phi-points', '4', '--threads', '2', '--orientation-chunk', '64', '--full-only', '--close']
controls = ['--analytic-facet-average', '--analytic-shadow-control', 'circular',
            '--analytic-mean-cache', str(root/'means.cache'), '--analytic-control-weights', str(root/'weights.tsv')]
def run(name, shape, extra):
    output = root/name/'result'
    output.parent.mkdir(parents=True, exist_ok=True)
    command = base+shape+extra+['--output', str(output)]
    (output.parent/'command.json').write_text(json.dumps(command, indent=2)+'\n')
    with (output.parent/'run.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
    return output, np.loadtxt(output/'result.dat', skiprows=1)
records = []
for name, shape in [
    ('prism', ['--particle', '1', '2', '1', '--ri', '1.3116', '0']),
    ('asymmetric_absorbing', ['--particle-file', str(repo/'tests/data/resident_asymmetric.particle'),
                             '--resize-dmax-um', '1', '--ri', '1.3243', '.0004'])]:
    raw_path, raw = run(name+'_raw', shape, ['--incoherent'])
    out, actual = run(name+'_controls', shape, ['--incoherent']+controls+['--analytic-facet-samples'])
    summary = np.genfromtxt(out/'result_analytic_facets.tsv', names=True, delimiter='\t')
    samples = np.genfromtxt(out/'result_analytic_facet_samples.tsv', names=True, delimiter='\t')
    expected = raw[:,2]+weights[:,0]*(summary['reflection_mean']-summary['sampled_reflection'])+weights[:,1]*(summary['shadow_mean']-summary['sampled_shadow'])
    scale = np.maximum(abs(raw[:,2]), 1e-100)
    error = max(np.max(abs(actual[:,2]-expected)/scale), np.max(abs(actual[:,3:]-raw[:,3:])/scale[:,None]))
    assert error < a.tolerance, (name, 'control identity', error)
    assert np.any(summary['sampled_reflection']>0), 'CPU control samples were skipped in incoherent mode'
    for i, angle in enumerate(theta):
        rows = samples[np.isclose(samples['theta_deg'], angle, atol=1e-10)]
        for sample, column in [('raw_M11','raw_M11'), ('reflection_control','sampled_reflection'), ('shadow_control','sampled_shadow')]:
            weighted = np.sum(rows['weight']*rows[sample])
            assert abs(weighted-summary[column][i])/scale[i] < a.tolerance, (name, sample, angle)
    # C itself is the same scalar function for either target, never a
    # coherent approximation substituted for the incoherent physical matrix.
    coh, unused = run(name+'_coherent_control', shape, controls)
    c = np.genfromtxt(coh/'result_analytic_facets.tsv', names=True, delimiter='\t')
    for field in ['reflection_mean','shadow_mean','sampled_reflection','sampled_shadow']:
        np.testing.assert_allclose(summary[field], c[field], rtol=1e-12, atol=1e-13)
    if name=='prism':
        paired, upper = run(name+'_paired_prefix', shape, ['--incoherent']+controls+['--paired-depth','6','--orientation-prefix-output','64'])
        # A standalone command for the lower depth / first 64 orientations.
        saved = base[:]
        base[base.index('--max-reflections')+1]='6'
        lower, reference = run(name+'_lower', shape, ['--incoherent']+controls)
        np.testing.assert_allclose(np.loadtxt(paired/'result_depth6.dat',skiprows=1),reference,rtol=a.tolerance,atol=1e-12)
        base[base.index('--sobol-seed')+1]='64'
        prefix, reference = run(name+'_lower_prefix', shape, ['--incoherent']+controls)
        np.testing.assert_allclose(np.loadtxt(paired/'result_prefix64_depth6.dat',skiprows=1),reference,rtol=a.tolerance,atol=1e-12)
        base[:] = saved
    records.append({'case': name, 'max_all16_scaled_M11': float(error)})
    print('PASS', name, 'incoherent target + intensity controls', error, flush=True)
(root/'results.json').write_text(json.dumps(records, indent=2)+'\n')
