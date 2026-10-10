#!/usr/bin/env python3
"""Independent NumPy/SciPy audit of native fullauto and saved worker matrices.

Python is a validation dependency, never a fullauto runtime dependency.
--dense also executes the saved native commands on the entire requested grid.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
import math
from pathlib import Path
import subprocess

import numpy as np
from scipy.interpolate import CubicSpline


def option(args, name, count=1):
    ids = [i for i, value in enumerate(args) if value == name]
    if not ids:
        return None
    i = ids[-1]
    return args[i+1] if count == 1 else args[i+1:i+1+count]


def replace(args, name, value):
    args = args.copy()
    i = max(i for i, item in enumerate(args) if item == name)
    args[i+1] = str(value)
    return args


def fnv(data):
    result = 14695981039346656037
    for byte in data:
        result = ((result ^ byte) * 1099511628211) & ((1 << 64)-1)
    return format(result, 'x')


def matrix(path):
    result = np.loadtxt(path, skiprows=1, ndmin=2)
    assert result.shape[1] == 18 and np.isfinite(result).all(), path
    return result


def run(root, reference, dense=False, dense_workers=1):
    spec = importlib.util.spec_from_file_location('python_campaign_reference', reference)
    py = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(py)
    state = json.loads((root/'fullauto_status.json').read_text())
    assert state['status'] == 'converged', state
    calibration = json.loads((root/'calibration.json').read_text())
    evaluated = np.loadtxt(root/'theta_evaluated.csv', ndmin=1)
    support = np.loadtxt(root/'theta_support.csv', ndmin=1)
    requested = np.loadtxt(root/'requested_theta.csv', ndmin=1)
    jobs = []
    for folder in sorted((root/'jobs').iterdir()):
        if not (folder/'complete.txt').exists():
            continue
        args = (folder/'command.args').read_text().splitlines()
        metadata = (folder/'complete.txt').read_text().split()
        data_path = folder/'result/result.dat'
        summary_path = folder/'result/result_analytic_facets.tsv'
        assert fnv(data_path.read_bytes()) == metadata[2], folder
        if metadata[3] != 'none':
            assert fnv(summary_path.read_bytes()) == metadata[3], folder
        log = (folder/'stdout.log').read_text()
        for label in ('Hard tree-limit hits:', 'Orientations still incomplete:'):
            assert int(log.split(label, 1)[1].split()[0]) == 0, folder
        n, seed = map(int, option(args, '--sobol-seed', 2))
        jobs.append(dict(folder=folder, args=args, count=n, seed=seed,
                         depth=int(option(args, '--max-reflections')),
                         phi=int(option(args, '--phi-points')),
                         weighted=option(args, '--analytic-control-weights') is not None,
                         seconds=float(metadata[1]), data=matrix(data_path),
                         summary=np.genfromtxt(summary_path, names=True, delimiter='\t')
                         if summary_path.exists() else None))

    def select(seeds, count, depth, phi, weighted, angles):
        selected = []
        for seed in seeds:
            found = [j for j in jobs if j['seed'] == seed and j['count'] == count
                     and j['depth'] == depth and j['phi'] == phi and j['weighted'] == weighted
                     and ('--analytic-facet-average' in j['args']) == state.get('analytic_controls', bool(calibration['weights']))
                     and ('--mirror-gamma' in j['args']) == state.get('mirror_gamma', False)
                     and '--haar-mirror-audit' not in j['args']
                     and len(j['data']) == len(angles)
                     and np.allclose(j['data'][:, 0], angles, rtol=0, atol=1e-10)]
            assert len(found) == 1, (seed, count, depth, phi, len(found))
            selected.append(found[0])
        return selected

    has_controls=bool(calibration['weights'])
    if has_controls:
        weight_path=root/Path(calibration['weights']).name
        if not weight_path.exists():weight_path=Path(calibration['weights'])
        weights = np.loadtxt(weight_path, skiprows=1, ndmin=2)
    else:
        pilot_job=next(j for j in jobs if j['seed']==py.VALIDATION_SEEDS[0] and j['count']==calibration['pilot_count'] and j['phi']==calibration['phi_candidates'][0]['phi'] and '--analytic-facet-average' not in j['args'])
        weights=np.column_stack([pilot_job['data'][:,0],np.ones((len(pilot_job['data']),2))])
    weight_difference=0.
    calibration_checks = []
    target = state['target_depth']
    for candidate in calibration['phi_candidates']:
        phi = candidate['phi']
        train = select(py.TRAIN_SEEDS, calibration['pilot_count'], max(1, target-4), phi, False, weights[:, 0])
        valid = select(py.VALIDATION_SEEDS, calibration['pilot_count'], max(1, target-4), phi, False, weights[:, 0])
        td = np.stack([j['data'] for j in train])
        vd = np.stack([j['data'] for j in valid])
        if has_controls:
            coefficients, info = py.fit_control_weights(td, vd, [j['summary'] for j in train], [j['summary'] for j in valid])
            adjusted = py.reweight(vd, [j['summary'] for j in valid], coefficients)
        else:
            coefficients=weights[:,1:];info=dict(rows_fitted=0,controls=False);adjusted=vd
        score = py.estimate_cost_score(adjusted, sum(j['seconds'] for j in valid))
        assert np.isclose(score, candidate['score'], rtol=1e-11, atol=1e-18)
        if phi == state['phi']:
            weight_difference = float(np.max(abs(coefficients-weights[:, 1:])))
            assert weight_difference < 1e-11, weight_difference
        calibration_checks.append(dict(phi=phi, python_score=score, native_score=candidate['score'], **info))
    assert min(calibration_checks, key=lambda x: x['python_score'])['phi'] == state['phi']

    depths = sorted(set([max(1, target-4), max(1, target-2), target]))
    seeds = state['seeds']
    counts = state['counts_per_seed']
    level_data, costs, raw_levels = [], [], []
    for level, (depth, count) in enumerate(zip(depths, counts)):
        upper = select(seeds, count, depth, state['phi'], has_controls, evaluated)
        data = np.stack([j['data'] for j in upper])
        cost = np.mean([j['seconds'] for j in upper])/count
        raw = (upper, None)
        if level:
            lower = select(seeds, count, depths[level-1], state['phi'], has_controls, evaluated)
            data[:, :, 2:] -= np.stack([j['data'] for j in lower])[:, :, 2:]
            cost += np.mean([j['seconds'] for j in lower])/count
            raw = (upper, lower)
        level_data.append(data)
        costs.append(cost)
        raw_levels.append(raw)
    combined = level_data[0].copy()
    for level in level_data[1:]:
        combined[:, :, 2:] += level[:, :, 2:]
    direct_metrics, direct_mean, _ = py.confidence(combined, tolerance=state['statistical_budget'])
    indices = np.searchsorted(evaluated, support)
    assert np.allclose(evaluated[indices], support, atol=1e-10, rtol=0)
    spline = CubicSpline(support, combined[:, indices, 2:], axis=1, bc_type='not-a-knot')
    reconstructed = np.zeros((8, len(requested), 18))
    reconstructed[:, :, 0] = requested
    reconstructed[:, :, 2:] = spline(requested)
    metrics, mean, ci = py.confidence(reconstructed, tolerance=state['statistical_budget'])
    native = matrix(root/'mueller_fullauto.dat')
    mean_difference = float(np.max(abs(mean[:, 2:]-native[:, 2:])/mean[:, 2, None]))
    native_ci = np.loadtxt(root/'all_mueller_confidence.csv', skiprows=1, delimiter=',')[:, 1:]
    ci_difference = float(np.max(abs(ci-native_ci)))
    delta = spline(evaluated)-combined[:, :, 2:]
    guards = (abs(delta.mean(axis=0))+py.T95_7*delta.std(axis=0, ddof=1)/math.sqrt(8))/direct_mean[:, 2, None]
    guard_max = float(guards.max())
    saved_guards = np.loadtxt(root/'theta_validation.csv', skiprows=1, delimiter=',')[:, 3]
    guard_difference = float(np.max(abs(guards.max(axis=1)-saved_guards)))
    direct_difference = float(np.max(abs(direct_mean[:, 2:]-matrix(root/'mueller_evaluated.dat')[:, 2:])/direct_mean[:, 2, None]))
    assert max(mean_difference, ci_difference, guard_difference, direct_difference) < 1e-12
    assert abs(max(metrics['max_mueller95_scaled_M11'], direct_metrics['max_mueller95_scaled_M11'])-state['max_pointwise95_scaled_M11']) < 1e-12
    assert abs(guard_max-state['max_guard_residual95_scaled_M11']) < 1e-12
    eps = state['statistical_budget']+state['interpolation_budget']+state.get('mirror_budget', 0)
    assert np.isclose(state['interpolation_budget'], eps/4, atol=1e-16, rtol=0)
    assert metrics['confidence_pass'] and direct_metrics['confidence_pass']
    assert guard_max <= eps/4 and state['pass_streak'] >= 2
    assert state.get('max_mirror_residual95_scaled_M11', 0) <= state.get('mirror_budget', 0)
    report = dict(passed=True, mirror_gamma=state.get('mirror_gamma', False), epsilon=eps, requested_theta=len(requested), evaluated_theta=len(evaluated),
                  support_theta=len(support), final_counts_per_seed=counts, phi=state['phi'],
                  matrix_difference_scaled_M11=mean_difference, confidence_absolute_difference=ci_difference,
                  guard_absolute_difference=guard_difference, evaluated_matrix_difference_scaled_M11=direct_difference,
                  statistical95_scaled_M11=state['max_pointwise95_scaled_M11'], guard_residual95_scaled_M11=guard_max,
                  control_weight_absolute_difference=weight_difference, calibration_checks=calibration_checks,
                  python_next_allocation_if_requested=py.allocation(level_data, costs, counts, direct_mean[:, 2], state['statistical_budget']*.85),
                  complete_jobs_audited=len(jobs), python_reference=str(reference.resolve()))

    probe = Path(__file__).resolve().parent/'.build/fullauto_probe'
    if probe.exists():
        allocation_dir = root/'independent_allocation'
        allocation_dir.mkdir(exist_ok=True)
        lines = [f"{state['statistical_budget']*.85:.17g} {len(level_data)}"]
        for index, (level, cost, count) in enumerate(zip(level_data, costs, counts)):
            lines.append(f'{count} {cost:.17g}')
            for seed in range(8):
                path = allocation_dir/f'level{index}_seed{seed}.dat'
                np.savetxt(path, level[seed], fmt='%.17g', header='theta area sixteen Mueller entries', comments='')
                lines.append(str(path))
        manifest = allocation_dir/'manifest.txt'
        manifest.write_text('\n'.join(lines)+'\n')
        result = subprocess.run([str(probe), '--allocation', str(manifest)], text=True, capture_output=True, check=True)
        native_allocation = [int(n) for n in result.stdout.split()]
        assert native_allocation == report['python_next_allocation_if_requested']
        report['native_next_allocation_if_requested'] = native_allocation

    if dense:
        destination = root/'independent_dense'
        destination.mkdir(exist_ok=True)
        dense_weights = destination/'weights.tsv'
        np.savetxt(dense_weights, np.column_stack([requested]+[np.interp(requested, weights[:, 0], weights[:, j]) for j in (1, 2)]),
                   fmt='%.17g', comments='', header='theta_deg reflection_weight shadow_weight')
        dense_levels = []
        for level, (upper, lower) in enumerate(raw_levels):
            def compute_seed(seed_index):
                pair = []
                for label, group in [('upper', upper), ('lower', lower)]:
                    if group is None:
                        continue
                    item = group[seed_index]
                    folder = destination/f'level{level}_{label}_seed{item["seed"]}'
                    args = replace(item['args'], '--theta-grid-file', root/'requested_theta.csv')
                    if has_controls:
                        args = replace(args, '--analytic-control-weights', dense_weights)
                        args = replace(args, '--analytic-mean-cache', destination/'means.cache')
                    args = replace(args, '--output', folder/'result')
                    if not (folder/'result/result.dat').exists():
                        folder.mkdir(exist_ok=True)
                        (folder/'command.args').write_text('\n'.join(args)+'\n')
                        with (folder/'stdout.log').open('w') as log:
                            subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, check=True)
                    log = (folder/'stdout.log').read_text()
                    for label in ('Hard tree-limit hits:', 'Orientations still incomplete:'):
                        assert int(log.split(label, 1)[1].split()[0]) == 0, folder
                    full = matrix(folder/'result/result.dat')
                    assert np.allclose(full[:, 0], requested, rtol=0, atol=1e-10)
                    pair.append(full[:, 2:])
                return pair[0] if len(pair) == 1 else pair[0]-pair[1]
            with ThreadPoolExecutor(max_workers=dense_workers) as pool:
                values = list(pool.map(compute_seed, range(8)))
            dense_levels.append(np.stack(values))
        dense_samples = sum(dense_levels)
        dense_mean = dense_samples.mean(axis=0)
        actual_error = abs(mean[:, 2:]-dense_mean)/dense_mean[:, 0, None]
        dense_delta = reconstructed[:, :, 2:]-dense_samples
        dense_bound = (abs(dense_delta.mean(axis=0))+py.T95_7*dense_delta.std(axis=0, ddof=1)/math.sqrt(8))/dense_mean[:, 0, None]
        report['independent_dense'] = dict(matrix_error_scaled_M11=float(actual_error.max()),
                                         residual_plus95_scaled_M11=float(dense_bound.max()),
                                         all_requested_nodes_pass=bool(np.all(dense_bound <= eps/4)))
        assert report['independent_dense']['all_requested_nodes_pass'], report['independent_dense']
        np.savetxt(root/'python_vs_native.csv', np.column_stack([requested, mean[:, 2], dense_mean[:, 0],
                     mean[:, 3], dense_mean[:, 1], actual_error.max(axis=1), ci.max(axis=1)]), delimiter=',',
                   header='theta_deg,native_M11,python_dense_M11,native_M12,python_dense_M12,max_interp_error_scaled_M11,max_stat95_scaled_M11', comments='')
    (root/'independent_fullauto_audit.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    parser.add_argument('--reference', type=Path, default=Path(__file__).resolve().parents[1]/'scripts/gpu_campaign.py')
    parser.add_argument('--dense', action='store_true')
    parser.add_argument('--dense-workers', type=int, default=1)
    settings = parser.parse_args()
    assert settings.dense_workers > 0
    run(settings.output.resolve(), settings.reference, settings.dense, settings.dense_workers)
