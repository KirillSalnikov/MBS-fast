#!/usr/bin/env python3
"""Opt-in CUDA integration probe. Run only on a coordinated test GPU.

Compares all 16 matrix elements, paired depth/prefix outputs and integrals.
No production process is stopped and no timing is interpreted as a speedup.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time


def numeric_rows(path):
    result = []
    for line in path.read_text().splitlines():
        try:
            row = [float(value) for value in line.split()]
        except ValueError:
            continue
        if len(row) >= 17:
            assert all(math.isfinite(value) for value in row), path
            result.append(row)
    assert result, path
    return result


def compare_integrals(left_text, right_text):
    """Physical quantities use their own scale; error estimates are fractions.

    A relative tolerance *of a small relative-error estimate* needlessly
    multiplies two small tolerances. Compare those estimates in absolute
    dimensionless units, at the same 1e-6 scale as the matrix probe.
    """
    readers = [csv.DictReader(io.StringIO(text), delimiter='\t')
               for text in (left_text, right_text)]
    physical = {'Cext', 'Csca', 'Cabs', 'Qext', 'Qsca', 'Qabs',
                'angular_scattering_integral', 'projected_area'}
    estimates = {'Csca_error_estimate_rel', 'partial_integral_error_estimate_rel'}
    exact = {'method', 'label', 'status', 'full_angular_coverage'}
    assert readers[0].fieldnames == readers[1].fieldnames, 'different integral schema'
    assert set(readers[0].fieldnames or ()) == physical | estimates | exact, 'unexpected integral schema'
    rows = [list(reader) for reader in readers]
    assert rows[0] and len(rows[0]) == len(rows[1]), 'different integral row counts'
    differences = []
    for left, right in zip(*rows):
        row_diff = {}
        for field in exact:
            assert left[field] == right[field], (field, left[field], right[field])
        for field in physical | estimates:
            u, v = float(left[field]), float(right[field])
            if field in estimates and math.isnan(u) and math.isnan(v):
                row_diff[field] = {'available': False}
                continue
            assert math.isfinite(u) and math.isfinite(v), (field, u, v)
            delta = abs(u-v)
            if field in estimates:
                assert delta <= 1e-6, (field, u, v, delta)
                row_diff[field] = {'absolute_difference': delta, 'absolute_tolerance': 1e-6}
            else:
                # Preserve the previous physical-value tolerance, including
                # its zero floor; +0 and -0 are physically identical.
                assert math.isclose(u, v, rel_tol=1e-6, abs_tol=1e-12), (field, u, v)
                row_diff[field] = {'absolute_difference': delta,
                                   'relative_difference': delta/max(abs(u), abs(v)) if u or v else 0.}
        differences.append(row_diff)
    return differences


def compare(left, right):
    a = {str(p.relative_to(left)): p for p in left.rglob('*.dat')}
    b = {str(p.relative_to(right)): p for p in right.rglob('*.dat')}
    assert a and a.keys() == b.keys(), 'different output matrices'
    reports = {}
    for name in a:
        x, y = numeric_rows(a[name]), numeric_rows(b[name])
        assert len(x) == len(y), name
        worst = 0.
        for u, v in zip(x, y):
            assert len(u) == len(v) and u[:-16] == v[:-16], name
            scale = max(abs(u[-16]), abs(v[-16]), 1e-300)
            worst = max(worst, max(abs(p-q) / scale for p, q in zip(u[-16:], v[-16:])))
        assert worst <= 1e-6, (name, worst)
        reports[name] = {'rows': len(x), 'max_all16_difference_scaled_M11': worst}
    aa = {str(p.relative_to(left)): p for p in left.rglob('*integrals.tsv')}
    bb = {str(p.relative_to(right)): p for p in right.rglob('*integrals.tsv')}
    assert aa and aa.keys() == bb.keys(), 'different integral outputs'
    integral_reports = {name: compare_integrals(aa[name].read_text(), bb[name].read_text()) for name in aa}
    return {'matrices': reports, 'integral_files': list(aa), 'integral_comparisons': integral_reports}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--timeout', type=float, default=120.)
    options = parser.parse_args()
    assert 0 < options.timeout <= 120
    binary = options.binary.resolve(strict=True)
    out = options.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    summary = {'binary': str(binary), 'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
               'production_changed': False, 'cases': [], 'comparisons': [], 'passed': False}
    try:
        for incoherent in (False, True):
            paths = []
            for overlap in (0, 1):
                case = out / (('incoherent' if incoherent else 'coherent') + '_overlap' + str(overlap))
                case.mkdir(); result = case / 'output'; result.mkdir()
                theta = case / 'theta.txt'
                theta.write_text(''.join(str(150 + i*2.5)+'\n' for i in range(9)))
                command = [str(binary), '--method', 'po', '--particle', '1', '5', '2.5',
                           '--refractive-index', '1.3116', '0', '--wavelength-um', '.532',
                           '--backend', 'cuda', '--po-shadow-model', 'projected', '--symmetry', '2', '6',
                           '--sobol-seed', '64', '53', '--threads', '4', '--close', '--haar-alpha',
                           '--theta-grid-file', str(theta), '--phi-points', '8', '--mirror-gamma',
                           '--profile-phases', '--analytic-mean-threads', '4', '--analytic-facet-average',
                           '--analytic-shadow-control', 'circular', '--orientation-pipeline',
                           '--orientation-chunk', '16', '--max-reflections', '6', '--paired-depth', '4',
                           '--cutoff-profile', 'balanced', '--output', str(result / 'result'),
                           '--orientation-prefix-output', '32', '--no-text-reports']
                if incoherent:
                    command.append('--incoherent')
                settings = {'MBS_GPU_OVERLAP_ANALYTIC_PACK': str(overlap), 'MBS_GPU_TIMING': '1',
                            'MBS_HOST_MEM_BUDGET_MB': '1024', 'MBS_GPU_MEM_FRACTION': '.10',
                            'MBS_GPU_BLOCK': '64', 'MBS_GPU_WARP_BEAMS': 'auto',
                            'MBS_GPU_UPLOAD_RESERVE_GROWTH': '0', 'MBS_GPU_RECYCLE_PREPARED': '0'}
                # Explicitly remove ambient experimental variables from this isolated probe.
                env = {k: v for k, v in os.environ.items() if not k.startswith('MBS_')}
                env.update(settings)
                record = {'command': command, 'environment': settings, 'incoherent': incoherent, 'overlap': overlap}
                (case / 'command.json').write_text(json.dumps(record, indent=2)+'\n')
                started = time.monotonic()
                with (case / 'stdout.log').open('w') as log:
                    process = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                    try:
                        rc = process.wait(timeout=options.timeout)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL); process.wait(timeout=5)
                        raise
                record.update(return_code=rc, wall_seconds=time.monotonic()-started)
                summary['cases'].append(record)
                assert rc == 0, (case, rc)
                log = (case / 'stdout.log').read_text()
                record['overlapped_chunks'] = log.count('Analytic packing overlap:')
                assert bool(record['overlapped_chunks']) == bool(overlap), 'requested overlap path not exercised'
                paths.append(result)
            summary['comparisons'].append({'incoherent': incoherent, **compare(*paths)})
        summary['passed'] = True
    finally:
        (out / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
