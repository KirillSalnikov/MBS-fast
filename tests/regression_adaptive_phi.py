#!/usr/bin/env python3
"""Idle-GPU test: adaptive phi against dense control, absorption and cropped theta.

Also verifies that a too-small refinement limit fails without final matrices.
"""
import argparse
import csv
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
    root = args.output or Path(tempfile.mkdtemp(prefix='mbs_adaptive_phi_regression_'))
    root.mkdir(parents=True, exist_ok=True)
    sizes = root/'sizes.dat'
    sizes.write_text('47.24199479082395\n413.3674544197096\n1181.049869770599\n')
    report = []
    for imaginary, first, last in [('0', '0', '180'), ('0.01', '20', '160')]:
        datasets = {}
        for mode in ['direct', 'adaptive', 'limit']:
            work = root/('imag'+imaginary)/mode
            work.mkdir(parents=True, exist_ok=True)
            env = {k:v for k,v in os.environ.items() if not k.startswith('MBS_')}
            env.update(MBS_GPU_MULTI='0', MBS_GPU_MULTI_MAX='1', MBS_GPU_MULTI_K_FULL='1',
                       MBS_SHARED_ORIENT_CHUNK='8', MBS_SHARED_BETA_GROUP='1',
                       MBS_GPU_MEM_FRACTION='.7', MBS_GPU_TRACE_INFLIGHT_BEAMS='32', MBS_GPU_TRACE_BATCH_BEAMS='2')
            if mode != 'direct':
                env['MBS_PHI_ADAPTIVE'] = '1'
            if mode == 'limit':
                env.update(MBS_PHI_ADAPT_MAX='300', MBS_PHI_ADAPT_M11_TOL='1e-14', MBS_PHI_ADAPT_POL_TOL='1e-14')
            cmd = [str(args.binary.resolve()), '--method', 'po', '--particle-file', str(args.particle.resolve()),
                   '--k-eq-list', str(sizes.resolve()), '--refractive-index', '1.4', imaginary,
                   '--wavelength-um', '.532', '--max-reflections', '8', '--euler-grid', '4', '4',
                   '--scattering-grid', first, last, '8192' if mode == 'direct' else '600', '18',
                   '--backend', 'cuda', '--gpu-trace-prefilter', '--threads', '4',
                   '--beam-cutoff-jones', '.001', '--beam-cutoff-area', '.002',
                   '--trace-cutoff-importance', '.0001', '--trace-max-beams', '20000',
                   '--allow-experimental-environment', '--close', '--output', str((work/'result').resolve())]
            with (work/'run.log').open('w') as log:
                proc = subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120)
            files = sorted((work/'result').glob('result_keq*.dat'))
            if mode == 'limit':
                assert proc.returncode != 0 and not files
                assert 'Adaptive phi did not converge at MAX' in (work/'run.log').read_text()
                continue
            assert proc.returncode == 0, (work/'run.log').read_text()
            assert len(files) == 3
            datasets[mode] = {p.name: matrix(p) for p in files}
            if mode == 'adaptive':
                reports = list(work.rglob('*phi_convergence.csv'))
                assert len(reports) == 1
                rows = list(csv.DictReader(reports[0].open()))
                assert len(rows) == 3*19*2  # 16 orientations in two chunks
                assert all(r['status'] == 'converged' and 300 <= int(r['nphi']) <= 9600 for r in rows)
        intensity = polarization = 0.0
        for name, ref in datasets['direct'].items():
            for a,b in zip(ref, datasets['adaptive'][name]):
                assert a[:2] == b[:2]
                intensity = max(intensity, abs(a[2]-b[2])/max(abs(a[2]), 1e-300))
                polarization = max(polarization, max(abs(x/a[2]-y/b[2]) for x,y in zip(a[2:],b[2:])))
        assert intensity < .01 and polarization < .005, (imaginary, intensity, polarization)
        report.append(dict(imaginary=imaginary, theta=[first,last], m11_max_relative=intensity,
                           normalized_mueller_max=polarization))
    (root/'summary.json').write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print('PASS; artifacts:', root)


if __name__ == '__main__':
    main()
