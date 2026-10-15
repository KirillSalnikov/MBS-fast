#!/usr/bin/env python3
"""Native reference checks for final integrals, prefixes, poles and recovery."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess

import numpy as np

T95 = 2.3646242510102993


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def run(command, log, env):
    log.with_suffix('.args').write_text('\n'.join(map(str, command))+'\n')
    with log.open('w') as stream:
        subprocess.run(list(map(str, command)), stdout=stream, stderr=subprocess.STDOUT,
                       env=env, check=True, timeout=300)


def matrix(path):
    return np.loadtxt(path, skiprows=1, ndmin=2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--legacy-binary', type=Path)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TASK_MBS_WORKER=str(args.binary.resolve()))
    reports = []
    cases = [
        ('partial', 'po', '0', '170', '180', True),
        ('absorbing_partial', 'po', '.02', '0', '25', True),
        ('absorbing_full', 'po', '.02', '0', '180', True),
        ('incoherent', 'po', '0', '0', '180', False),
        ('go_absorbing', 'go', '.02', '0', '180', True),
    ]
    if args.legacy_binary:
        cases.append(('legacy_partial', 'po', '0', '170', '180', True))
    for name, method, imag, start, end, coherent in cases:
        worker = args.legacy_binary.resolve() if name.startswith('legacy') else args.binary.resolve()
        env['TASK_MBS_WORKER'] = str(worker)
        dest = out/name
        common = ['--method', method, '--backend', 'cpu', '--particle', '1', '.2', '.1',
                  '--ri', '1.3116', imag, '--wavelength-um', '1', '--max-reflections', '4',
                  '--threads', '1', '--symmetry', '1', '1', '--cutoff-profile', 'off', '--close']
        if not coherent:
            common += ['--incoherent']
        phi = '1' if method == 'go' else '4'
        command = [args.probe.resolve()]+common+['--fullauto', '.02', '--phi-points', phi,
                   '--fullauto-initial', '256',
                   '--fullauto-orientation-block', '64', '--fullauto-mirror', 'off',
                   '--fullauto-analytic-controls', 'off', '--fullauto-resources', 'off',
                   '--fullauto-workers', 'persistent',
                   '--scattering-grid', start, end, phi, '8', '--output', dest]
        if method == 'po':
            command += ['--fullauto-min-correction', '64']
        run(command, out/(name+'.log'), env)
        summary = table(dest/'mueller_fullauto_integrals.tsv')[0]
        reps = table(dest/'integral_replicas.tsv')
        mean = matrix(dest/'frozen_matrix.dat')
        bins = len(mean)-2 if method == 'go' else len(mean)
        angular = float(np.dot(mean[:bins, 1], mean[:bins, 2]))
        if method == 'go':
            angular += mean[-2:, 2].sum()
        np.testing.assert_allclose(float(summary['angular_scattering_integral']), angular, rtol=2e-13)
        for index, rep in enumerate(reps):
            values = matrix(dest/('replica'+str(index)+'.dat'))
            expected = np.dot(values[:bins, 1], values[:bins, 2])
            if method == 'go':
                expected += values[-2:, 2].sum()
            np.testing.assert_allclose(float(rep['angular_scattering_integral']), expected, rtol=2e-13)
        # Independent standalone full-prefix ray traces, not a replay of the
        # controller's prefix/provenance aggregation implementation.
        seed = int(reps[0]['seed'])
        levels = [(0, 1, 256, 1), (1, 1, 128, 2), (1, -1, 128, 1),
                  (2, 1, 64, 4), (2, -1, 64, 2)] if method == 'po' else [(0, 1, 256, 4)]
        expected_ext = 0
        expected_area = 0
        for level, sign, count, depth in levels:
            ref = out/(name+'_ref_'+str(level)+'_'+str(sign))
            reference_common = common[:]
            reference_common[reference_common.index('--max-reflections')+1] = str(depth)
            replay = [worker]+reference_common+[
                '--sobol-seed', str(count), str(seed), '--scattering-grid', start, end, phi, '8', '--output', ref]
            if method == 'po':
                replay += ['--haar-alpha']
            else:
                replay += ['--go-fullauto-worker']
            run(replay, out/(ref.name+'.log'), env)
            row = table(ref/(ref.name+'_integrals.tsv'))[0]
            if not name.startswith('legacy') or level == 2:
                expected_ext += sign*float(row['Cext'])
            if not level:
                expected_area = float(row.get('projected_area', float(row['Cext'])/float(row['Qext'])))
        if name.startswith('legacy'):
            # Legacy lower-depth scalar output is known to be incorrect; only
            # a direct upper/target sample is admissible for extinction.
            ref = out/(name+'_target_reference')
            ref_common = common[:]
            ref_common[ref_common.index('--max-reflections')+1] = '4'
            run([worker]+ref_common+['--sobol-seed', '64', str(seed), '--haar-alpha',
                '--scattering-grid', start, end, '4', '8', '--output', ref], out/(ref.name+'.log'), env)
            expected_ext = float(table(ref/(ref.name+'_integrals.tsv'))[0]['Cext'])
            assert json.loads((dest/'integral_summary.json').read_text())['legacy_paired_extinction_fallback']
        np.testing.assert_allclose(float(reps[0]['Cext']), expected_ext, rtol=2e-12)
        np.testing.assert_allclose(float(reps[0]['projected_area']), expected_area, rtol=2e-12)
        if imag != '0' and start == '0' and end == '180':
            np.testing.assert_allclose(float(summary['Csca']), angular, rtol=2e-13)
            np.testing.assert_allclose(float(summary['Cabs']), float(summary['Cext'])-angular, rtol=2e-12)
        elif imag != '0':
            assert np.isnan(float(summary['Csca'])) and np.isnan(float(summary['Cabs']))
            assert summary['status'] == 'partial_angular_integral'
        else:
            assert summary['Cext'] == summary['Csca'] and float(summary['Cabs']) == 0
        for row in table(dest/'integral_confidence.tsv'):
            if row['quantity'].startswith('Q'):
                # Central efficiency is a ratio of means, never a mean of ratios.
                wanted = 'C'+row['quantity'][1:]
                np.testing.assert_allclose(float(row['mean']), float(summary[wanted])/float(summary['projected_area']), rtol=2e-12)
                continue
            values = np.array([float(rep[row['quantity']]) for rep in reps])
            np.testing.assert_allclose(float(row['mean']), values.mean(), rtol=2e-12)
            np.testing.assert_allclose(float(row['upper95']), values.mean()+T95*values.std(ddof=1)/np.sqrt(8), rtol=2e-12)
        markers = {p: p.read_bytes() for p in (dest/'jobs').glob('*/complete.txt')}
        previous = (dest/'mueller_fullauto_integrals.tsv').read_bytes()
        run(command, out/(name+'_resume.log'), env)
        assert all(p.read_bytes() == value for p, value in markers.items())
        assert previous == (dest/'mueller_fullauto_integrals.tsv').read_bytes()
        reports.append(dict(case=name, status=summary['status'], angular_integral=angular,
                            native_full_prefix_reference=True, resumed_blocks_preserved=len(markers)))
    (out/'results.json').write_text(json.dumps(reports, indent=2)+'\n')
    print(json.dumps(reports, indent=2))


if __name__ == '__main__':
    main()
