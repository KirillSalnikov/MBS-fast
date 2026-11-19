#!/usr/bin/env python3
"""Exercise concurrent CPU/GPU tracing with facet controls, paired depth and a prefix."""
import argparse
import csv
import math
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', required=True, type=Path)
parser.add_argument('--output', required=True, type=Path)
args = parser.parse_args()
root = args.output.resolve()
root.mkdir(parents=True, exist_ok=True)
orientations = 4096
common = [str(args.binary.resolve()), '--method', 'po', '--backend', 'cuda',
          '--particle', '2', '10', '9.85', '--ri', '1.3243', '0', '--wavelength-um', '.355',
          '--symmetry', '1', '6', '--haar-alpha', '--sobol-seed', str(orientations), '2029',
          '--max-reflections', '18', '--paired-depth', '16', '--orientation-prefix-output', '2048',
          '--scattering-grid', '170', '180', '4', '6', '--threads', '4', '--orientation-pipeline',
          '--analytic-facet-average', '--analytic-shadow-control', 'circular',
          '--cutoff-profile', 'balanced', '--no-text-reports', '--full-only', '--close']
for backend in ('cpu', 'packet-hybrid'):
    output = root / backend / 'result'
    output.parent.mkdir(exist_ok=False)
    with (output.parent / 'run.log').open('w') as log:
        subprocess.run(common + ['--po-trace-backend', backend, '--output', str(output)],
                       check=True, stdout=log, stderr=subprocess.STDOUT, timeout=180)

log = (root / 'packet-hybrid' / 'run.log').read_text()
assert 'producer_CUDA=0' in log and 'GPU convex resident direct: accepted=' in log
assert 'complete chunk fallback' not in log
with (root / 'packet-hybrid/result/result_hybrid_ranges_rank0.tsv').open() as source:
    ranges = sorted(csv.DictReader(source, delimiter='\t'), key=lambda row: int(row['start']))
cursor = 0
for row in ranges:
    assert int(row['start']) == cursor and int(row['count']) > 0
    cursor += int(row['count'])
assert cursor == orientations and {r['producer'] for r in ranges} == {'CPU', 'GPU'}

for product in ('result', 'result_depth16', 'result_prefix2048', 'result_prefix2048_depth16'):
    rows = []
    for backend in ('cpu', 'packet-hybrid'):
        path = root / backend / 'result' / (product + '.dat')
        rows.append([list(map(float, line.split())) for line in path.read_text().splitlines()[1:]])
    assert len(rows[0]) == len(rows[1]) == 7
    for a, b in zip(*rows):
        assert len(a) == len(b) == 18 and a[:2] == b[:2] and min(a[2], b[2]) > 0
        assert all(math.isfinite(v) for v in a + b)
        assert max(abs(x-y) / a[2] for x, y in zip(a[2:], b[2:])) < 1e-6
        assert max(abs(x/a[2]-y/b[2]) for x, y in zip(a[3:], b[3:])) < 1e-6
    scalars = []
    for backend in ('cpu', 'packet-hybrid'):
        path = root / backend / 'result' / (product + '_integrals.tsv')
        with path.open() as source:
            scalars.append(next(csv.DictReader(source, delimiter='\t')))
    for key in ('Cext', 'Csca', 'Cabs', 'Qext', 'Qsca', 'Qabs', 'projected_area', 'angular_scattering_integral'):
        a, b = [float(v[key]) for v in scalars]
        assert math.isclose(a, b, rel_tol=1e-6, abs_tol=1e-12), (product, key, a, b)
print('PASS disjoint coverage, serial facet controls, upper/lower matrices, prefixes and integrals')
