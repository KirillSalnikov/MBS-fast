#!/usr/bin/env python3
"""Independently check paired Haar mirror pilots and converged on/off runs."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from audit_fullauto import matrix, option

T95 = 2.3646242510102993


def run(root):
    on, off = root/'fullauto_on', root/'fullauto_off'
    state = json.loads((on/'fullauto_status.json').read_text())
    baseline = json.loads((off/'fullauto_status.json').read_text())
    assert state['status'] == baseline['status'] == 'converged'
    assert state['mirror_gamma'] and not baseline['mirror_gamma']
    assert state['seeds'] == baseline['seeds'] and state['phi'] == baseline['phi']
    jobs = []
    for path in (on/'jobs').glob('*/command.args'):
        if not (path.parent/'complete.txt').exists():
            continue
        args = path.read_text().splitlines()
        n, seed = map(int, option(args, '--sobol-seed', 2))
        jobs.append(dict(args=args, count=n, seed=seed, depth=int(option(args, '--max-reflections')),
                         phi=int(option(args, '--phi-points')), data=matrix(path.parent/'result/result.dat')))

    checks = []
    for path in sorted(on.glob('mirror_verification_*.json')):
        saved = json.loads(path.read_text())
        key = path.stem.removeprefix('mirror_verification_')
        grid = f'theta_{key}.csv'

        def samples(depth, paired):
            result = []
            for seed in saved['seeds']:
                found = [j for j in jobs if j['count'] == saved['base_count_per_seed']
                         and j['seed'] == seed and j['depth'] == depth and j['phi'] == saved['phi']
                         and Path(option(j['args'], '--theta-grid-file')).name == grid
                         and (('--haar-mirror-audit' in j['args']) == paired)
                         and (('--mirror-gamma' in j['args']) != paired)]
                assert len(found) == 1, (depth, seed, paired, len(found))
                result.append(found[0]['data'][:, 2:])
            return np.stack(result)

        reference = samples(state['target_depth'], True).mean(axis=0)[:, 0]
        worst = 0
        for level in saved['levels']:
            delta = samples(level['depth'], False)-samples(level['depth'], True)
            bound = (abs(delta.mean(axis=0))+T95*delta.std(axis=0, ddof=1)/np.sqrt(8))/reference[:, None]
            value = float(bound.max())
            assert abs(value-level['max_residual_plus95_scaled_M11']) < 1e-12
            worst = max(worst, value)
        ml_bound = (2*len(saved['levels'])-1)*worst
        assert abs(ml_bound-saved['mlmc_residual_bound_plus95_scaled_M11']) < 1e-12
        assert ml_bound <= saved['budget']
        checks.append(dict(theta_count=saved['theta_count'], mlmc_residual_bound=ml_bound,
                           paired_native_ratio=saved['paired_native_seconds']/saved['reduced_native_seconds']))

    a, b = matrix(on/'mueller_fullauto.dat'), matrix(off/'mueller_fullauto.dat')
    assert np.array_equal(a[:, 0], b[:, 0])
    error = abs(a[:, 2:]-b[:, 2:])/b[:, 2, None]
    c1 = np.loadtxt(on/'all_mueller_confidence.csv', delimiter=',', skiprows=1)[:, 1:]
    c2 = np.loadtxt(off/'all_mueller_confidence.csv', delimiter=',', skiprows=1)[:, 1:]
    allowed = c1*a[:, 2, None]/b[:, 2, None]+c2+state['interpolation_budget']+baseline['interpolation_budget']
    assert np.all(error <= allowed), float((error-allowed).max())
    timing_on = json.loads((root/'fullauto_on_timing.json').read_text())['wall_seconds']
    timing_off = json.loads((root/'fullauto_off_timing.json').read_text())['wall_seconds']
    report = dict(passed=True, paired_checks=checks, on_seconds=timing_on, off_seconds=timing_off,
                  end_to_end_speedup=timing_off/timing_on, max_on_off_difference_scaled_M11=float(error.max()),
                  on_off_agreement_within_sum_of_pointwise95_and_interpolation_budgets=True,
                  on_counts=state['counts_per_seed'], off_counts=baseline['counts_per_seed'],
                  on_evaluated_theta=state['evaluated_theta'], off_evaluated_theta=baseline['evaluated_theta'])
    (root/'mirror_comparison_audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    run(parser.parse_args().root)
