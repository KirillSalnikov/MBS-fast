#!/usr/bin/env python3
"""A collinear edge subdivision must not change scattering or extinction.

The same hexagonal prism has 6 or 48 cap vertices. The latter exercises the
prepared CPU fallback, also when a CUDA batch contains large apertures.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np


def particle(path, subdivisions):
    corners = [(math.cos(i * math.pi / 3), math.sin(i * math.pi / 3), -1.)
               for i in range(6)]
    bottom = [tuple((1 - t / subdivisions) * a[k] + t / subdivisions * b[k]
                    for k in range(3))
              for a, b in zip(corners, corners[1:] + corners[:1])
              for t in range(subdivisions)]
    top = [(x, y, 1.) for x, y, _ in bottom]
    faces = [list(reversed(bottom)), top]
    for i in range(6):
        faces.append([bottom[(i * subdivisions + j) % len(bottom)]
                      for j in range(subdivisions + 1)]
                     + [top[(i * subdivisions + j) % len(top)]
                        for j in range(subdivisions, -1, -1)])
    path.write_text('0\n0\n180 360\n\n' + '\n\n'.join(
        '\n'.join('%.17g %.17g %.17g' % v for v in face)
        for face in faces) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary', type=Path)
    parser.add_argument('--backend', choices=['cpu', 'cuda'], default='cpu')
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    temporary = tempfile.TemporaryDirectory(prefix='mbs-apertures-') if args.out is None else None
    root = Path(temporary.name) if temporary else args.out.resolve()
    root.mkdir(parents=True, exist_ok=True)
    for split in (1, 8):
        particle(root / f'hex{split}.particle', split)
    sizes = root / 'sizes.txt'
    sizes.write_text('20\n40\n')
    results = {}
    tolerance = 2e-6 if args.backend == 'cpu' else 8e-5

    def run(name, split, loss, size, extra_env=None):
        dest = root / name
        assert not dest.exists(), dest
        cmd = [str(args.binary.resolve()), '--method', 'po', '--backend', args.backend,
               '--particle-file', str(root / f'hex{split}.particle'),
               '--ri', '1.3116', str(loss), '--wavelength-um', '.532',
               '--max-reflections', '4', '--euler-grid', '3', '3',
               '--scattering-grid', '0', '180', '4', '36', '--cutoff-profile', 'off',
               '--threads', '2', '--output', str(dest), '--close']
        cmd += ['--k-eq', str(size)] if size else ['--k-eq-list', str(sizes)]
        if extra_env:
            cmd.append('--allow-experimental-environment')
        env = dict(os.environ, OMP_NUM_THREADS='2', MBS_GPU_MEM_FRACTION='.4')
        env.update(extra_env or {})
        (root / (name + '_command.json')).write_text(json.dumps(cmd, indent=2))
        with (root / (name + '.log')).open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=120, env=env)
        if args.backend == 'cuda' and split == 8:
            assert 'GPU aperture fallback:' in (root / (name + '.log')).read_text()
        outputs = {}
        for k in ([size] if size else [20, 40]):
            stem = name + ('' if size else f'_keq{k}')
            table = np.loadtxt(dest / (stem + '.dat'), skiprows=1)
            with (dest / (stem + '_integrals.tsv')).open() as stream:
                integral = next(csv.DictReader(stream, delimiter='\t'))
            outputs[k] = table, {key: float(integral[key]) for key in ('Cext', 'Cabs')}
        return outputs

    def check(name, actual, reference):
        a, integral = actual
        b, ref_integral = reference
        assert np.all(np.isfinite(a)), name
        np.testing.assert_allclose(a[:, :2], b[:, :2], rtol=0, atol=1e-12)
        error = float(np.max(abs(a[:, 2:] - b[:, 2:]) / np.maximum(b[:, 2, None], 1.)))
        integral_error = max(abs(integral[k] - ref_integral[k])
                             / max(abs(ref_integral[k]), 1.) for k in integral)
        assert error < tolerance, (name, 'Mueller', error)
        assert integral_error < tolerance, (name, 'integrals', integral_error)
        results[name] = dict(max_all16_over_max_M11_1=error, integral_error=integral_error)
        print('PASS', name, results[name], flush=True)

    for loss in (0., .02):
        reference = {}
        for k in (20, 40):
            name = f'loss{loss}_coarse{k}'
            reference[k] = run(name, 1, loss, k)[k]
            name = f'loss{loss}_fine{k}'
            check(name, run(name, 8, loss, k)[k], reference[k])
        variants = [('default', {})]
        if args.backend == 'cuda':
            variants += [('explicit', {'MBS_SHARED_GPU_EXPLICIT_SCALE': '1'}),
                         ('fused', {'MBS_GPU_MULTI_K_FULL': '1'})]
        for variant, env in variants:
            name = f'loss{loss}_sweep_{variant}'
            actual = run(name, 8, loss, None, env)
            for k in (20, 40):
                check(name + f'_{k}', actual[k], reference[k])
    (root / 'results.json').write_text(json.dumps(results, indent=2) + '\n')


if __name__ == '__main__':
    main()
