#!/usr/bin/env python3
"""Verify dot/slash aliases preserve fullauto optics, replay and recovery."""
import argparse
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

from test_fullauto_memory_runtime import Suite, capture, compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    suite = Suite(SimpleNamespace(binary=args.binary, output=args.output, rounds=1,
                                  timeout=90, tolerance=1e-14))
    cube = (Path(__file__).resolve().parents[1] / 'examples/cube.particle').read_bytes()
    visibility = b'6\n1 1\n1 1\n1 1\n1 1\n1 1\n1 1\n'
    reference, reference_files = None, None
    try:
        for case in ['canonical', 'particle_dot', 'particle_double_slash', 'saved_nested_alias']:
            output = suite.root / case
            saved = output / 'geometry/deep/saved.particle'
            saved.parent.mkdir(parents=True)
            particle = output / 'particle.dat'
            particle.write_bytes(cube)
            Path(str(particle) + '.visibility').write_bytes(visibility)
            particle_argument = str(output) + ('/./particle.dat' if case == 'particle_dot'
                                               else '//particle.dat' if case == 'particle_double_slash'
                                               else '/particle.dat')
            saved_argument = str(output) + ('/geometry/./deep//saved.particle' if case == 'saved_nested_alias'
                                           else '/geometry/deep/saved.particle')
            command = suite.command('po_coherent', output, 'memory', checkpoint=True)
            start = command.index('--particle')
            command[start:start + 4] = ['--particle-file', particle_argument, '--resize-dmax-um', '5']
            command += ['--save-geometry', saved_argument]
            suite.run(case + '_fresh', command)
            actual = capture(output)
            files = [saved.read_bytes(), Path(str(saved) + '.visibility').read_bytes()]
            replay = json.loads((output / 'controller_command.json').read_text())['argv']
            assert replay[replay.index('--particle-file') + 1] == particle_argument, 'transport canonicalization changed recorded input CLI'
            assert replay[replay.index('--save-geometry') + 1] == saved_argument, 'transport canonicalization changed recorded save CLI'
            suite.run(case + '_resume', command)
            resume = compare(actual, capture(output), suite.args.tolerance)
            assert files == [saved.read_bytes(), Path(str(saved) + '.visibility').read_bytes()]
            if reference is None:
                reference, reference_files = actual, files
            comparison = compare(reference, actual, suite.args.tolerance)
            assert files == reference_files, 'aliases changed exported geometry/visibility bytes'
            comparison.update({'resume': resume, 'original_cli_preserved': True,
                               'geometry_sha256': hashlib.sha256(files[0]).hexdigest(),
                               'visibility_sha256': hashlib.sha256(files[1]).hexdigest()})
            suite.receipt['comparisons'][case] = comparison
            suite.save()
            print(case + ': optics, geometry bytes, original CLI and resume PASS', flush=True)
        suite.receipt['passed'] = True
    except Exception as error:
        suite.receipt['passed'] = False
        suite.receipt['failure'] = repr(error)
        raise
    finally:
        suite.save()
    print(json.dumps({'passed': True, 'receipt': str(suite.root / 'receipt.json')}, indent=2))


if __name__ == '__main__':
    main()
