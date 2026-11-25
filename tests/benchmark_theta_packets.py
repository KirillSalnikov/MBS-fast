#!/usr/bin/env python3
"""Time four 2-angle PO calls and one 8-angle call on the same native sampling.

Pass the solver, an eight-line theta file and native solver arguments after --.
Set the GPU and CPU affinity as for a production worker. Each timing includes
process startup. The caller supplies identical Sobol and optical settings for
both layouts. Results are written to stdout as JSON.
"""
import argparse
import json
import math
from pathlib import Path
import statistics
import subprocess
import tempfile
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path)
    parser.add_argument('angles', type=Path)
    parser.add_argument('native_args', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    native = args.native_args
    if native[:1] == ['--']:
        native = native[1:]
    if '--sobol-seed' not in native:
        parser.error('supply --sobol-seed COUNT SEED to fix the orientation sample')
    if any(flag in native for flag in ('--theta-grid-file', '--output', '-o')):
        parser.error('the benchmark supplies the angle grid and output directory')
    angles = args.angles.read_text().split()
    theta = list(map(float, angles))
    if len(theta) != 8 or theta != sorted(set(theta)) or not all(0 <= t <= 180 for t in theta):
        parser.error('expected eight increasing scattering angles between 0 and 180 degrees')

    runs = []
    reference = None
    max_difference = 0.0
    with tempfile.TemporaryDirectory(prefix='mbs-theta-benchmark-') as directory:
        work = Path(directory)
        grid = work / 'theta.txt'
        for trial, width in enumerate((2, 8, 8, 2)):
            combined = []
            seconds = 0.0
            for start in range(0, 8, width):
                grid.write_text('\n'.join(angles[start:start + width]) + '\n')
                dest = work / f'trial{trial}_row{start}'
                command = [str(binary), *native, '--theta-grid-file', str(grid), '--output', str(dest)]
                begin = time.monotonic()
                run = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                seconds += time.monotonic() - begin
                if run.returncode:
                    raise RuntimeError(run.stdout[-4096:].decode(errors='replace'))
                with (dest / (dest.name + '.dat')).open() as matrix:
                    next(matrix)
                    rows = [list(map(float, line.split())) for line in matrix if line.strip()]
                if len(rows) != width or any(len(row) != 18 or row[2] <= 0
                                             or not all(map(math.isfinite, row)) for row in rows):
                    raise ValueError('unexpected PO matrix')
                combined.extend(rows)
            if reference is None:
                reference = combined
            for want, got in zip(reference, combined):
                if want[0] != got[0]:
                    raise ValueError('scattering angles differ across layouts')
                max_difference = max(max_difference, max(abs(x-y) for x,y in zip(want[2:], got[2:])) / want[2])
            runs.append({'width': width, 'seconds': seconds})
    speedup = statistics.mean(r['seconds'] for r in runs if r['width'] == 2) / statistics.mean(
        r['seconds'] for r in runs if r['width'] == 8)
    print(json.dumps(dict(runs=runs, speedup=speedup, max_difference_scaled_M11=max_difference), indent=2))
    if max_difference > 1e-8:
        raise SystemExit('matrix comparison exceeded 1e-8 of M11')


if __name__ == '__main__':
    main()
