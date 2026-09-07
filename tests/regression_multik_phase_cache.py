#!/usr/bin/env python3
"""CUDA integration regression: cached/uncached multi-k vs per-size FP64.

Run on an idle GPU, for example:
  CUDA_VISIBLE_DEVICES=0 python3 tests/regression_multik_phase_cache.py \
    --binary gpu/bin/mbs_po_gpu_double --particle /path/to/particle.dat
All outputs are retained in a new temporary directory or --output.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


def matrix(path):
    with path.open() as stream:
        assert next(stream).split()[0] == 'ScAngle', path
        data = [[float(v) for v in line.split()] for line in stream if line.strip()]
    assert len(data) == 19, (path, len(data))
    assert all(len(r) == 18 and all(math.isfinite(v) for v in r) for r in data), path
    return data


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary', required=True, type=Path)
    ap.add_argument('--particle', required=True, type=Path)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    output = args.output or Path(tempfile.mkdtemp(prefix='mbs_phase_regression_'))
    output.mkdir(parents=True, exist_ok=True)
    sizes = output / 'sizes.dat'
    sizes.write_text('47.24199479082395\n413.3674544197096\n1181.049869770599\n')
    report = []
    for imaginary in ['0', '0.01']:
        datasets = {}
        for mode, switches in [
            ('legacy', {'MBS_GPU_MULTI_K_FULL': '0'}),
            ('uncached', {'MBS_GPU_MULTI_K_FULL': '1', 'MBS_GPU_MULTI_K_PHASE_CACHE': '0'}),
            ('vertex_array', {'MBS_GPU_MULTI_K_STREAM': '0'}),
            ('cached', {})
        ]:
            work = output / ('imag' + imaginary) / mode
            work.mkdir(parents=True, exist_ok=True)
            env = {k: v for k, v in os.environ.items() if not k.startswith('MBS_')}
            env.update(MBS_GPU_MULTI='0', MBS_GPU_MULTI_MAX='1',
                       MBS_SHARED_ORIENT_CHUNK='8', MBS_GPU_MEM_FRACTION='.7',
                       MBS_GPU_TRACE_INFLIGHT_BEAMS='32', MBS_GPU_TRACE_BATCH_BEAMS='2')
            env.update(switches)
            cmd = [str(args.binary.resolve()), '--method', 'po',
                   '--particle-file', str(args.particle.resolve()),
                   '--k-eq-list', str(sizes.resolve()),
                   '--refractive-index', '1.4', imaginary, '--wavelength-um', '.532',
                   '--max-reflections', '8', '--euler-grid', '4', '4',
                   '--scattering-grid', '0', '180', '72', '18',
                   '--backend', 'cuda', '--gpu-trace-prefilter', '--threads', '4',
                   '--beam-cutoff-jones', '.001', '--beam-cutoff-area', '.002',
                   '--trace-cutoff-importance', '.0001', '--trace-max-beams', '20000',
                   '--allow-experimental-environment', '--close',
                   '--output', str((work/'result').resolve())]
            with (work/'run.log').open('w') as log:
                subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=120)
            log_text = (work/'run.log').read_text()
            if mode != 'legacy':
                assert 'fused multi-k kernel enabled' in log_text
                assert 'falling back to per-size GPU path' not in log_text
            files = sorted((work/'result').glob('result_keq*.dat'))
            assert len(files) == 3
            datasets[mode] = {p.name: matrix(p) for p in files}
        for mode in ['uncached', 'vertex_array', 'cached']:
            worst = 0.0
            for name, ref in datasets['legacy'].items():
                test = datasets[mode][name]
                for a, b in zip(ref, test):
                    assert abs(a[0] - b[0]) < 1e-10
                    scale = max(abs(a[2]), 1e-300)
                    worst = max(worst, max(abs(x-y)/scale for x, y in zip(a[2:], b[2:])))
            assert worst < 1e-9, (imaginary, mode, worst)
            report.append(dict(imaginary=imaginary, mode=mode, max_mueller_error_over_M11=worst))
    (output/'summary.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print('PASS; artifacts:', output)


if __name__ == '__main__':
    main()
