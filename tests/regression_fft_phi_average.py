#!/usr/bin/env python3
"""Compare compact phi moments against dense FFT interpolation on an idle GPU.

Tests odd/even direct grids, two absorption values, multiple orientation chunks,
all Mueller components and both poles. Does NOT certify phi-grid convergence.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
from regression_multik_phase_cache import matrix


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--particle', type=Path, required=True)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    root = args.output or Path(tempfile.mkdtemp(prefix='mbs_fft_mean_regression_'))
    root.mkdir(parents=True, exist_ok=True)
    sizes = root / 'sizes.dat'
    sizes.write_text('47.24199479082395\n413.3674544197096\n1181.049869770599\n')
    report = []
    for imaginary in ['0', '0.01']:
        for nphi, factor in [(72, 2), (300, 4), (70, 1)]:
            datasets = {}
            for mode in ['fft', 'mean']:
                work = root / ('imag' + imaginary) / str(nphi) / mode
                work.mkdir(parents=True, exist_ok=True)
                env = {k: v for k, v in os.environ.items() if not k.startswith('MBS_')}
                env.update(MBS_GPU_MULTI='0', MBS_GPU_MULTI_MAX='1', MBS_GPU_MULTI_K_FULL='1',
                           MBS_SHARED_ORIENT_CHUNK='8', MBS_GPU_MEM_FRACTION='.7',
                           MBS_GPU_TRACE_INFLIGHT_BEAMS='32', MBS_GPU_TRACE_BATCH_BEAMS='2')
                if mode == 'mean':
                    env['MBS_FFT_PHI_AVERAGE_ONLY'] = '1'
                cmd = [str(args.binary.resolve()), '--method', 'po',
                       '--particle-file', str(args.particle.resolve()),
                       '--k-eq-list', str(sizes.resolve()), '--refractive-index', '1.4', imaginary,
                       '--wavelength-um', '.532', '--max-reflections', '8', '--euler-grid', '4', '4',
                       '--scattering-grid', '0', '180', str(nphi), '18', '--fft-factor', str(factor),
                       '--backend', 'cuda', '--gpu-trace-prefilter', '--threads', '4',
                       '--beam-cutoff-jones', '.001', '--beam-cutoff-area', '.002',
                       '--trace-cutoff-importance', '.0001', '--trace-max-beams', '20000',
                       '--allow-experimental-environment', '--close', '--output', str((work/'result').resolve())]
                with (work/'run.log').open('w') as log:
                    subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT,
                                   check=True, timeout=120)
                log_text = (work/'run.log').read_text()
                # The existing fused FFT dispatcher falls back at factor=1;
                # that direct per-size reference is valid for this edge case.
                if factor > 1 or mode == 'mean':
                    assert 'falling back to per-size GPU path' not in log_text
                if mode == 'mean':
                    assert 'FFT phi average-only: direct Nphi=' in log_text
                    assert 'GPU FFT multik angular interpolation:' not in log_text
                files = sorted((work/'result').glob('result_keq*.dat'))
                assert len(files) == 3
                datasets[mode] = {p.name: matrix(p) for p in files}
            worst = 0.0
            for name, reference in datasets['fft'].items():
                for a, b in zip(reference, datasets['mean'][name]):
                    assert a[:2] == b[:2]
                    worst = max(worst, max(abs(x-y)/max(abs(a[2]), 1e-300)
                                           for x, y in zip(a[2:], b[2:])))
            assert worst < 1e-9, (imaginary, nphi, factor, worst)
            report.append(dict(imaginary=imaginary, nphi=nphi, factor=factor,
                               max_mueller_error_over_M11=worst))
    # Safety switches must not silently bypass requested error checks/refinement.
    for switch in ['MBS_FFT_CHECK', 'MBS_FFT_ADAPTIVE_PHI', 'MBS_FFT_GLOBAL_REFINE', 'MBS_FFT_THETA_FACTOR']:
        guard_env = dict(env, **{switch: '2' if switch.endswith('FACTOR') else '1'})
        guard = subprocess.run(cmd, env=guard_env, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, universal_newlines=True, timeout=120)
        assert guard.returncode != 0 and 'FFT phi average-only' in guard.stdout, (switch, guard.stdout)
        (root / (switch + '.log')).write_text(guard.stdout)
    (root/'summary.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print('PASS; artifacts:', root)


if __name__ == '__main__':
    main()
