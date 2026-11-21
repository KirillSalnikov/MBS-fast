"""Check CUDA facet upload ordering on rotated joined rosettes.

Run on one visible GPU with four available CPU cores:
    python3 tests/check_nonconvex_trace_stream.py /path/to/mbs_po_gpu
"""
import csv
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def solve(binary, root, name, mode, env):
    target = root / name
    command = [
        binary, '--method', 'po', '--backend', 'cuda',
        '--particle', '3', '2.9585957810762125', '4.57505088811751',
        '--rosette-joined', '--ri', '1.3194', '0', '--wavelength-um', '0.4',
        '--symmetry', '2', '2', '--haar-alpha', '--sobol-seed', '8', '6043',
        '--theta-grid-file', str(root / 'theta'), '--phi-points', '4',
        '--max-reflections', '12', '--cutoff-profile', 'off',
        '--po-shadow-model', 'projected', '--threads', '4', '--po-prepare-threads', '4',
        '--allow-experimental-environment', '--profile-phases', mode, '--output', str(target), '--close',
    ]
    log = subprocess.check_output(command, cwd=str(root), env=env,
                                  stderr=subprocess.STDOUT, universal_newlines=True)
    # The submission notice precedes memory allocation. Timing records are
    # emitted only after the kernel and the result transfer have completed.
    records = [dict(item.split('=', 1) for item in line.split() if '=' in item)
               for line in log.splitlines() if line.startswith('GPU trace timing:')]
    timed_kernel = any(int(r['items']) > 0 and int(r['candidates']) > 0
                       and float(r['small_ms']) + float(r['large_ms']) > 0 for r in records)
    if (mode == '--gpu-trace-prefilter' and not timed_kernel
            or mode == '--no-gpu-trace-prefilter' and records):
        raise RuntimeError('completed CUDA trace batches do not match requested backend')
    print('%s: %d completed CUDA trace batches' % (name, len(records)), flush=True)
    with (target / (name + '_phase_timings.tsv')).open() as stream:
        times = {r['phase']: float(r['seconds']) for r in csv.DictReader(stream, delimiter='\t')}
    if times.get('gpu_diffraction_including_packing', 0) <= 0:
        raise RuntimeError('CUDA diffraction did not execute')
    with (target / (name + '.dat')).open() as stream:
        next(stream)
        return [[float(x) for x in line.split()] for line in stream if line.strip()]


def check(binary):
    angles = (0, 5, 32, 78, 142, 160, 170, 176, 180)
    env = {k: v for k, v in os.environ.items() if not k.startswith('MBS_GPU_TRACE_')}
    env.update(OMP_NUM_THREADS='4', OMP_THREAD_LIMIT='4', OMP_DYNAMIC='FALSE',
               MBS_GPU_TRACE_NONBLOCKING_STREAM='1', MBS_GPU_TRACE_BATCH_BEAMS='128',
               MBS_GPU_TRACE_TIMING='1')
    with tempfile.TemporaryDirectory(prefix='mbs-trace-stream-') as directory:
        root = Path(directory)
        (root / 'theta').write_text(''.join(str(a) + '\n' for a in angles))
        reference = solve(binary, root, 'cpu_facets', '--no-gpu-trace-prefilter', env)
        for repeat in range(3):
            actual = solve(binary, root, 'gpu_facets_' + str(repeat), '--gpu-trace-prefilter', env)
            if len(reference) != len(angles) or len(actual) != len(angles):
                raise RuntimeError('incomplete angular matrix')
            for angle, x, y in zip(angles, actual, reference):
                if (len(x) != 18 or len(y) != 18 or x[:2] != y[:2] or abs(y[0]-angle) > 1e-10
                        or not all(map(math.isfinite, x + y)) or min(x[2], y[2]) <= 0):
                    raise RuntimeError('invalid matrix row at theta=%g' % angle)
                for index, (test, expected) in enumerate(zip(x[2:], y[2:])):
                    if abs(test - expected) / y[2] > 1e-5:
                        raise RuntimeError('trace mismatch at theta=%g M%d%d repetition %d'
                                           % (angle, index//4 + 1, index%4 + 1, repeat+1))


if __name__ == '__main__':
    check(str(Path(sys.argv[1]).resolve()))
