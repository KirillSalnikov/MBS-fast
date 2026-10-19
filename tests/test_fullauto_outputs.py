#!/usr/bin/env python3
"""Real PO/GO exports: preserve direct checks, match budgets, verify recovery hashes."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import numpy as np


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def matrix(table):
    return np.array([[float(r['M%d%d' % (i, j)]) for i in range(1, 5) for j in range(1, 5)] for r in table])


def run(command, path, ok=True):
    with path.open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=240)
    if ok:
        assert result.returncode == 0, path.read_text()[-4000:]
    else:
        assert result.returncode != 0 and 'completed output changed' in path.read_text(), path.read_text()[-4000:]


def verify(folder, method, final=True):
    full = rows(folder/'mueller_fullauto_annotated.csv')
    direct = rows(folder/'mueller_evaluated_annotated.csv')
    raw = np.loadtxt(folder/'mueller_evaluated.dat', skiprows=1, ndmin=2)
    np.testing.assert_array_equal(matrix(direct), raw[:, 2:])
    if final:
        output = np.loadtxt(folder/'mueller_fullauto.dat', skiprows=1, ndmin=2)
        np.testing.assert_array_equal(matrix(full), output[:, 2:])
    guide = json.loads((folder/'matrix_output_schema.json').read_text())
    assert guide['full_rows'] == len(full) and guide['direct_rows'] == len(direct)
    assert all(r['value_source'] in ('direct', 'direct_bin') for r in direct)
    assert all(r['max_combined95'] == '' and r['max_interpolation95_or_allowance'] == '' for r in direct)
    if method == 'po':
        lookup = {float(r['theta_deg']): r for r in direct}
        for r in full:
            angle = float(r['theta_deg'])
            assert int(r['has_direct_calculation']) == (angle in lookup)
            assert (r['value_source'] == 'direct') == (r['point_type'] == 'support')
            assert int(r['is_support']) == (r['point_type'] == 'support')
            if r['point_type'] == 'support':
                np.testing.assert_array_equal(matrix([r]), matrix([lookup[angle]]))
            elif r['point_type'] == 'interpolated':
                assert r['interpolation_estimate_kind'] == 'unmeasured_allowance'
            else:
                assert r['interpolation_estimate_kind'] in ('paired_guard95', 'common_prefix_guard95')
        assert {r['point_type'] for r in full} == {'support', 'guard', 'interpolated'}
        budget = rows(folder/('pointwise_budget_production.csv' if final else 'pointwise_budget.csv'))
        np.testing.assert_allclose([float(r['max_combined95']) for r in full],
                                   [float(r['output_total95']) for r in budget], rtol=0, atol=1e-15)
        if final:
            confidence = rows(folder/'all_mueller_confidence.csv')
            expected = [max(float(v) for k, v in r.items() if k != 'theta_deg') for r in confidence]
            np.testing.assert_allclose([float(r['max_statistical95_scaled_M11']) for r in full], expected, rtol=0, atol=1e-15)
    else:
        assert full == direct
        assert all(r['point_type'] == 'finite_bin' and r['is_support'] == '0' for r in full)
        assert (folder/'mueller_forward.dat').exists() and (folder/'mueller_backward.dat').exists()
    return {'full_rows': len(full), 'direct_rows': len(direct), 'types': sorted({r['point_type'] for r in full})}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve();root.mkdir(parents=True, exist_ok=True)
    report = {}
    for method in ('po', 'go'):
        folder = root/method
        command = [str(args.binary.resolve()), '--method', method, '--backend', 'cpu', '--particle', '1', '.2', '.1',
                   '--ri', '1.3116', '0', '--wavelength-um', '1', '--max-reflections', '2', '--threads', '2',
                   '--fullauto', '.05', '--fullauto-final-checks', 'off', '--fullauto-initial', '256',
                   '--max-orientations', '65536', '--fullauto-max-rounds', '16', '--output', str(folder), '--close']
        if method == 'po':
            command += ['--incoherent', '--fullauto-resources', 'off', '--fullauto-cutoffs', 'off', '--cutoff-profile', 'off',
                        '--fullauto-mirror', 'off', '--fullauto-rotation-symmetry', 'off', '--symmetry', '1', '1',
                        '--fullauto-analytic-controls', 'off', '--phi-points', '4', '--fullauto-pilot', '32',
                        '--fullauto-min-correction', '64', '--fullauto-theta-start', '3', '--scattering-grid', '0', '25', '4', '32']
        else:
            command += ['--scattering-grid', '0', '180', '1', '6']
        (root/(method+'.args')).write_text('\n'.join(command)+'\n')
        run(command, root/(method+'.log'))
        assert json.loads((folder/'fullauto_status.json').read_text())['status'] == 'converged'
        report[method] = verify(folder, method)
        snapshot = Path(json.loads((folder/'intermediate_current.json').read_text())['snapshot_folder'])
        report[method]['snapshot'] = verify(snapshot, method, final=False)
        before = {name: (folder/name).read_bytes() for name in ['mueller_fullauto.dat', 'mueller_evaluated.dat',
                  'mueller_fullauto_annotated.csv', 'mueller_evaluated_annotated.csv', 'matrix_output_schema.json']}
        history = (folder/'fullauto_history.jsonl').read_text().count('"event":"worker_start"')
        run(command, root/(method+'_replay.log'))
        assert 'no calculation repeated' in (root/(method+'_replay.log')).read_text()
        assert all((folder/n).read_bytes() == value for n, value in before.items())
        assert history == (folder/'fullauto_history.jsonl').read_text().count('"event":"worker_start"')
        for name in ['mueller_fullauto_annotated.csv', 'mueller_evaluated_annotated.csv', 'matrix_output_schema.json']:
            (folder/name).write_bytes(before[name]+b'\n')
            try:
                run(command, root/(method+'_corrupt_'+name+'.log'), ok=False)
            finally:
                (folder/name).write_bytes(before[name])
        run(command, root/(method+'_restored_replay.log'))
        assert json.loads((folder/'fullauto_status.json').read_text())['status'] == 'converged'
        report[method]['completed_replay_unchanged'] = True
        report[method]['corrupted_annotations_rejected'] = True
    (root/'verification.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
