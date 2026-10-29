#!/usr/bin/env python3
"""Compare raw GPU paths, polygons, Jones amplitudes and phases to CPU tracing.

Uses the existing diagnostic comparator before diffraction. A CPU fallback or
an inactive comparator cannot pass. This is intentionally not a timing test.
"""
import argparse
import json
import os
import re
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True)
    repo = Path(__file__).resolve().parents[1]
    grid = root / 'theta.txt'
    grid.write_text('0\n170\n180\n')
    shapes = [
        ('column', ['--particle', '1', '398.1', '138.9', '--symmetry', '2', '6'], 8),
        ('droxtal', ['--particle', '4', '44.565', '--symmetry', '2', '6'], 14),
        ('asymmetric', ['--particle-file', str(repo / 'tests/data/resident_asymmetric.particle'),
                        '--resize-dmax-um', '25', '--symmetry', '1', '1'], 8),
    ]
    common = ['--method', 'po', '--backend', 'cuda', '--wavelength-um', '.532',
              '--ri', '1.3116', '0', '--threads', '2', '--haar-alpha',
              '--theta-grid-file', str(grid), '--phi-points', '4', '--sobol-seed', '32', '23',
              '--orientation-pipeline', '--orientation-chunk', '16', '--trace-cache-mib', '0',
              '--close', '--allow-experimental-environment']
    records = []
    for name, shape, depth in shapes:
        for schedule in ['compact', 'facets4', 'facets8', 'facets16', 'facets']:
            for cutoff in ['off', 'balanced']:
                dest = root / (name + '_' + schedule + '_' + cutoff)
                dest.mkdir()
                env = os.environ.copy()
                env.update(MBS_GPU_TRACE_CONVEX_QUEUE='1', MBS_GPU_TRACE_CONVEX_LEVELS='1',
                           MBS_GPU_TRACE_CONVEX_DIRECT='0', MBS_GPU_TRACE_CONVEX_VERIFY='1',
                           MBS_GPU_TRACE_HYBRID='0', MBS_GPU_TRACE_CONVEX_SCHEDULE=schedule,
                           MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB='64',
                           MBS_HOST_MEM_BUDGET_MB='2048', OMP_WAIT_POLICY='PASSIVE')
                # --po-trace-backend resident would override DIRECT=0 and bypass
                # the raw CPU/GPU comparator. Select the queue through env here.
                command = [str(args.binary.resolve())] + common + shape + [
                    '--max-reflections', str(depth), '--cutoff-profile', cutoff,
                    '--output', str(dest / 'result')]
                (dest / 'command.args').write_text('\n'.join(command) + '\n')
                (dest / 'environment.json').write_text(json.dumps({
                    k: v for k, v in env.items() if k.startswith(('MBS_', 'OMP_', 'CUDA_'))}, indent=2))
                with (dest / 'stdout.log').open('w') as log:
                    result = subprocess.run(command, env=env, stdout=log, stderr=log, timeout=120)
                text = (dest / 'stdout.log').read_text()
                assert result.returncode == 0, (dest, result.returncode)
                accepted = sum(map(int, re.findall(r'GPU convex resident queue: orientations=\d+ accepted=(\d+)', text)))
                checks = re.findall(r'GPU convex verification: beams=(\d+) equivalent_vertex_changes=\d+ max_error=([\deE.+-]+)', text)
                assert accepted == 32 and not re.search(r'CPU_fallback=[1-9]', text), (dest, accepted)
                assert checks and sum(int(n) for n, _ in checks) > 0, dest
                error = max(float(e) for _, e in checks)
                assert error <= 1e-8, (dest, error)
                record = dict(case=dest.name, beams=sum(int(n) for n, _ in checks), max_error=error)
                records.append(record)
                print(json.dumps(record), flush=True)
    (root / 'result.json').write_text(json.dumps(dict(passed=True, records=records), indent=2))


if __name__ == '__main__':
    main()
