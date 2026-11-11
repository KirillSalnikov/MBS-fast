#!/usr/bin/env python3
"""Check root/nested geometry exports against legacy fullauto calculations."""
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
    cube = Path(__file__).resolve().parents[1] / 'examples/cube.particle'
    try:
        for placement in ['root', 'nested']:
            results, files = {}, {}
            for mode in ['legacy', 'memory']:
                name = placement + '_' + mode
                output = suite.root / name
                saved = output / ('geometry/deep/saved.particle' if placement == 'nested' else 'saved.particle')
                saved.parent.mkdir(parents=True)
                command = suite.command('po_coherent', output, mode, checkpoint=True)
                start = command.index('--particle')
                command[start:start + 4] = ['--particle-file', str(cube), '--resize-dmax-um', '5']
                command += ['--save-geometry', str(saved)]
                suite.run(name, command)
                results[mode] = capture(output)
                sidecar = Path(str(saved) + '.visibility')
                assert saved.is_file() and sidecar.is_file(), 'geometry and sidecar must both be physically exported'
                files[mode] = [saved.read_bytes(), sidecar.read_bytes()]
                if mode == 'memory':
                    suite.run(name + '_resume', command)
                    compare(results[mode], capture(output), suite.args.tolerance)
                    assert files[mode] == [saved.read_bytes(), sidecar.read_bytes()], 'resume changed exported particle'
            assert files['memory'] == files['legacy'], 'memory transport changed saved geometry/visibility bytes'
            comparison = compare(results['memory'], results['legacy'], suite.args.tolerance)
            comparison['geometry_sha256'] = hashlib.sha256(files['memory'][0]).hexdigest()
            comparison['visibility_sha256'] = hashlib.sha256(files['memory'][1]).hexdigest()
            comparison['geometry_and_visibility_identical_to_legacy'] = True
            suite.receipt['comparisons'][placement] = comparison
            suite.save()
            print(placement + ': all16/CI and exported geometry/visibility match legacy PASS', flush=True)
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
