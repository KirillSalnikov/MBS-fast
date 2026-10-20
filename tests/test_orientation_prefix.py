#!/usr/bin/env python3
"""Validate embedded balanced prefixes against separate native Sobol passes.

Includes all 16 Mueller elements and Mij/M11, offset, controls, absorption,
paired-depth extinction, and the incoherent worker path. Production outputs
are never read or modified.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np


def matrix_error(actual, reference):
    assert actual.shape == reference.shape
    assert np.array_equal(actual[:, :2], reference[:, :2])
    scale = np.maximum(abs(reference[:, 2]), 1e-30)
    elements = float(np.max(abs(actual[:, 2:] - reference[:, 2:]) / scale[:, None]))
    ratios = float(np.max(abs(actual[:, 2:] / actual[:, 2, None]
                              - reference[:, 2:] / reference[:, 2, None])))
    return max(elements, ratios)


def integrals(path):
    with path.open() as stream:
        row = next(csv.DictReader(stream, delimiter='\t'))
    return {key: float(row[key]) for key in ('projected_area', 'Cext')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--backend', default='cpu')
    parser.add_argument('--tolerance', type=float)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--count', type=int, default=512)
    parser.add_argument('--phi', type=int, default=4)
    parser.add_argument('--case', default='all')
    parser.add_argument('--pipeline', action='store_true')
    parser.add_argument('--experimental', action='store_true')
    parser.add_argument('--trace-backend', choices=['cpu', 'resident', 'hybrid'])
    parser.add_argument('--mirror', action='store_true')
    parser.add_argument('--arbitrary', action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    binary = str(args.binary.resolve())
    count, prefix = args.count, max(64, args.count // 8)
    tolerance = args.tolerance if args.tolerance is not None else (1e-6 if args.backend == 'gpu' else 1e-10)
    grid = output / 'angles.csv'
    grid.write_text('0\n170\n175\n180\n')
    weights = output / 'weights.tsv'
    weights.write_text('theta_deg\treflection_weight\tshadow_weight\n'
                       '0\t0.37\t1.18\n170\t0.41\t0.9\n175\t0.32\t1.07\n180\t0.68\t1.2\n')
    base = [binary, '--method', 'po', '--backend', args.backend,
            '--particle', '1', '2', '1', '--wavelength-um', '.532',
            '--symmetry', '1', '1', '--threads', str(args.threads),
            '--haar-alpha', '--sobol-offset', '512', '--max-reflections', '8',
            '--theta-grid-file', str(grid), '--phi-points', str(args.phi),
            '--cutoff-profile', 'balanced', '--profile-phases', '--close']
    if args.pipeline:
        base.append('--orientation-pipeline')
    if args.experimental:
        base.append('--allow-experimental-environment')
    if args.trace_backend:
        base += ['--po-trace-backend', args.trace_backend]
    if args.mirror:
        base.append('--mirror-gamma')
    if args.arbitrary:
        from test_fullauto_arbitrary import asymmetric_tetrahedron
        particle = output / 'asymmetric.particle'
        asymmetric_tetrahedron(particle)
        shape = base.index('--particle')
        base[shape:shape + 4] = ['--particle-file', str(particle)]
    records = []
    cases = [('coherent', ['--ri', '1.3116', '0']),
             ('incoherent', ['--ri', '1.3116', '0', '--incoherent']),
             ('paired', ['--ri', '1.3116', '0', '--paired-depth', '6']),
             ('controls_absorbing', ['--ri', '1.3116', '.001', '--paired-depth', '6',
                 '--analytic-facet-average', '--analytic-shadow-control', 'circular',
                 '--analytic-control-weights', str(weights),
                 '--analytic-mean-cache', str(output / 'means.cache')])]
    for name, extra in cases:
        if args.case != 'all' and args.case != name:
            continue
        paths, timings = {}, {}
        for mode, n, flags in [('prefix', prefix, []), ('full', count, []),
                                ('shared', count, ['--orientation-prefix-output', str(prefix)])]:
            directory = output / (name + '_' + mode)
            command = base + extra + ['--sobol-seed', str(n), '11'] + flags + ['--output', str(directory)]
            (output / (directory.name + '.args')).write_text('\n'.join(command) + '\n')
            started = time.monotonic()
            with (output / (directory.name + '.log')).open('w') as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=240)
            timings[mode] = time.monotonic() - started
            paths[mode] = directory / directory.name
        comparisons = []
        for depth in (['', '_depth6'] if '--paired-depth' in extra else ['']):
            for mode, suffix in [('full', ''), ('prefix', '_prefix' + str(prefix))]:
                shared, separate = str(paths['shared']) + suffix + depth, str(paths[mode]) + depth
                error = matrix_error(np.loadtxt(shared + '.dat', skiprows=1),
                                     np.loadtxt(separate + '.dat', skiprows=1))
                assert error <= tolerance, (name, mode, depth, error, tolerance)
                a, b = integrals(Path(shared + '_integrals.tsv')), integrals(Path(separate + '_integrals.tsv'))
                integral_error = max(abs(a[k] - b[k]) / max(abs(b[k]), 1e-30) for k in b)
                assert integral_error < 1e-10, (name, mode, depth, a, b)
                if '--analytic-facet-average' in extra:
                    x, y = np.loadtxt(shared + '_analytic_facets.tsv', skiprows=1), np.loadtxt(separate + '_analytic_facets.tsv', skiprows=1)
                    # Setup timing is deliberately not a numerical observable.
                    columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11]
                    summary_error = float(np.max(abs(x[:, columns] - y[:, columns]) / np.maximum(abs(y[:, 1]), 1e-30)[:, None]))
                    assert summary_error <= tolerance, (name, summary_error)
                comparisons.append(dict(mode=mode, depth=depth, all16_and_ratios=error,
                                        scalar_integrals=integral_error))
        record = dict(case=name, seconds=timings, native_prefix_speedup=(timings['full'] + timings['prefix']) / timings['shared'], comparisons=comparisons)
        records.append(record)
        (output / 'comparison.json').write_text(json.dumps(dict(backend=args.backend, count=count, prefix=prefix, phi=args.phi, pipeline=args.pipeline, arbitrary=args.arbitrary, results=records), indent=2) + '\n')
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    main()
