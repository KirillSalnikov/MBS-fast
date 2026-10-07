#!/usr/bin/env python3
"""Compare complete Mueller grids across serial/parallel orientation-file runs."""
import argparse, math, os, subprocess, tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def read_grid(path):
    with path.open() as stream:
        next(stream)
        return [list(map(float, line.split())) for line in stream if line.strip()]

def compare(reference, candidate):
    assert len(reference) == len(candidate), 'different row counts'
    scale = max(abs(v) for row in reference for v in row[2:])
    worst = 0.0
    for a, b in zip(reference, candidate):
        assert len(a) == len(b) and a[:2] == b[:2], 'different scattering grids'
        for x, y in zip(a[2:], b[2:]):
            assert math.isfinite(y), 'nonfinite Mueller element'
            worst = max(worst, abs(x-y))
    assert worst <= 1e-12 * max(scale, 1e-20), (worst, scale)
    return worst / max(scale, 1e-20)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, default=ROOT/'cpu/bin/mbs_po_mpi')
    args = parser.parse_args()
    env = {k: v for k, v in os.environ.items() if not k.startswith('MBS_')}
    with tempfile.TemporaryDirectory(prefix='mbs_orientation_file_') as directory:
        work = Path(directory)
        orientations = work/'orientations.txt'
        orientations.write_text('0 0\n27.5 33.25\n90 0\n90 0\n120 240\n180 90\n')
        cases = [
            ('convex', ['1','4','4'], ['1.31','0'], []),
            ('absorbing', ['1','4','4'], ['1.53','0.0018'], ['--abs-points','all']),
            ('nonconvex', ['10','4','2.8','29.570400970655'], ['1.31','0'], []),
            ('incoherent', ['1','4','4'], ['1.53','0.0018'], ['--incoherent']),
        ]
        for name, particle, ri, options in cases:
            reference = None
            for mode, flags, threads in [('serial', [], '4'),
                                         ('parallel_one', ['--parallel-trace'], '1'),
                                         ('parallel_four', ['--parallel-trace'], '4')]:
                target = work/(name+'_'+mode)/'result'
                command = [str(args.binary.resolve()), '--method','po','--backend','cpu',
                           '--particle',*particle,'--refractive-index',*ri,
                           '--wavelength-um','0.532','--max-reflections','4',
                           '--orientation-file',str(orientations),'--threads',threads,
                           '--scattering-grid','0','180','6','36','--cutoff-profile','off',
                           '--output',str(target),'--close',*options,*flags]
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.with_suffix('.log').open('w') as log:
                    subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                   check=True, timeout=120)
                grid = read_grid(target/'result.dat')
                if reference is None:
                    reference = grid
                else:
                    print('PASS:',name,mode,'normalized Mueller difference',compare(reference,grid))
    print('Orientation-file acceleration regression passed')

if __name__ == '__main__':
    main()
