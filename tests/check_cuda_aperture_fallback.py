#!/usr/bin/env python3
"""Check paired CUDA/CPU regression runs containing an unsupported aperture.

Run the same particle, orientation/angle grids and cutoffs with --backend
cuda and --backend cpu, retaining run.log and result/result.dat in separate
directories. Then use --gpu-run DIR --cpu-run DIR. See
docs/GPU_APERTURE_FALLBACK_20260910.md for the production reproducer.
No third-party Python packages or private particle fixtures are required
to check existing outputs. Run this check after rebuilding the CUDA binary.
"""
import argparse
import math
from pathlib import Path
import re
import shlex


def read_run(folder):
    log = (folder / 'run.log').read_text()
    assert 'ERROR: calculation failed' not in log, folder
    grid = re.search(r'Orientation grid:\s*(\d+)\s*x\s*(\d+)\s*=\s*(\d+)', log)
    assert grid, 'Missing orientation grid: ' + str(folder)
    command = next(line[len('Command: '):] for line in log.splitlines()
                   if line.startswith('Command: '))
    with (folder / 'result/result.dat').open() as stream:
        header = next(stream).strip()
        rows = [list(map(float, line.split())) for line in stream if line.strip()]
    assert len(rows) == 181, folder
    for theta, row in enumerate(rows):
        assert len(row) == 18 and all(math.isfinite(x) for x in row), folder
        assert abs(row[0] - theta) < 1e-9 and row[2] > 0, folder
    return log, grid.groups(), shlex.split(command), header, rows


def physics(command):
    # Binary, backend, CPU affinity, tracing prefilter and output location may
    # differ; the physical inputs and sampling must be exactly the same.
    arities = {'--method': 1, '--particle-file': 1, '--k-eq': 1,
               '--wavelength-um': 1, '--refractive-index': 2,
               '--max-reflections': 1, '--orientation-diffraction-sampling': 1,
               '--scattering-grid': 4, '--beam-cutoff-jones': 1,
               '--beam-cutoff-area': 1, '--trace-cutoff-importance': 1,
               '--trace-max-beams': 1, '--trace-limit-retries': 1}
    result = {}
    for option, count in arities.items():
        pos = command.index(option)
        result[option] = command[pos + 1:pos + 1 + count]
    assert '--latitude-phi-grid' in command
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gpu-run', type=Path, required=True)
    parser.add_argument('--cpu-run', type=Path, required=True)
    parser.add_argument('--tolerance', type=float, default=1e-6)
    args = parser.parse_args()
    assert 0 < args.tolerance < 1
    gpu_log, gpu_grid, gpu_cmd, gpu_header, gpu = read_run(args.gpu_run)
    _, cpu_grid, cpu_cmd, cpu_header, cpu = read_run(args.cpu_run)
    assert gpu_cmd[gpu_cmd.index('--backend') + 1] == 'cuda'
    assert cpu_cmd[cpu_cmd.index('--backend') + 1] == 'cpu'
    assert gpu_grid == cpu_grid and gpu_header == cpu_header
    assert physics(gpu_cmd) == physics(cpu_cmd), 'Physics/sampling mismatch'
    assert 'GPU aperture fallback:' in gpu_log, 'Fallback was not exercised'
    vertices = [int(n) for n in re.findall(r'source-vertices=(\d+)', gpu_log)]
    assert any(n >= 32 for n in vertices), 'No oversized source contour tested'
    error = max(abs(g[j] - c[j]) / c[2]
                for g, c in zip(gpu, cpu) for j in range(2, 18))
    assert error <= args.tolerance, 'max |delta Mij|/M11 = %.17g' % error
    print('PASS: all 16 Mueller elements at 181 angles; max |delta Mij|/M11 = %.17g' % error)


if __name__ == '__main__':
    main()
