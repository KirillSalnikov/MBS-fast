#!/usr/bin/env python3
"""Compare real numeric products with/without duplicate reports; interrupt/resume fullauto."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command, log, ok=(0,)):
    with log.open('w') as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=240)
    assert result.returncode in ok, (result.returncode, log.read_text()[-3000:])
    return result.returncode


def reports(root):
    return sorted(str(p.relative_to(root)) for pattern in ('*_log.txt', '*_out.txt')
                  for p in root.rglob(pattern))


def compare_native(binary, root):
    base = [str(binary), '--backend', 'cpu', '--particle', '1', '.2', '.1',
            '--ri', '1.3116', '0', '--wavelength-um', '1', '--max-reflections', '4',
            '--sobol-seed', '64', '11', '--haar-alpha', '--symmetry', '1', '1',
            '--scattering-grid', '0', '180', '2', '8', '--threads', '1', '--close']
    variants = {
        'po': ['--method', 'po'],
        'po_incoherent_paired_prefix': ['--method', 'po', '--incoherent', '--paired-depth', '2',
                                      '--orientation-prefix-output', '32'],
        'po_noshadow': ['--method', 'po', '--no-shadow-output'],
        'go': ['--method', 'go'],
    }
    result = {}
    for name, options in variants.items():
        case = root/name
        case.mkdir()
        for mode in ('normal', 'compact'):
            (case/mode).mkdir()
            command = base+options+['--output', str(case/mode/'result')]
            if name == 'go':
                command.remove('--haar-alpha')
            if mode == 'compact':
                command += ['--no-text-reports']
            (case/(mode+'.args.json')).write_text(json.dumps(command, indent=2))
            run(command, case/(mode+'.stdout'))
        normal, compact = case/'normal/result', case/'compact/result'
        assert reports(normal), name
        assert not reports(compact), (name, reports(compact))
        products = [p for p in normal.rglob('*') if p.is_file()
                    and (p.suffix == '.dat' or p.name.endswith('_integrals.tsv')
                         or p.name.endswith('_integrals_revision.txt'))]
        assert products and any(p.name.endswith('_integrals.tsv') for p in products), name
        for path in products:
            other = compact/path.relative_to(normal)
            assert other.is_file() and other.read_bytes() == path.read_bytes(), (name, str(path))
        assert 'INTEGRAL' in (case/'compact.stdout').read_text().upper(), name
        removed = [normal/p for p in reports(normal)]
        result[name] = dict(equal_numeric_products=len(products),
                            omitted_files=len(removed),
                            omitted_bytes=sum(p.stat().st_size for p in removed))
    return result


def interrupt_persistent(binary, root):
    folder = root/'interrupted'
    command = [str(binary), '--method', 'po', '--backend', 'cpu', '--incoherent',
               '--particle', '1', '5', '2.5', '--ri', '1.3116', '0', '--wavelength-um', '.532',
               '--symmetry', '1', '1', '--max-reflections', '6', '--fullauto', '.04',
               '--threads', '1', '--fullauto-workers', 'persistent',
               '--fullauto-resources', 'off', '--fullauto-analytic-controls', 'off',
               '--fullauto-cutoffs', 'off', '--cutoff-profile', 'off',
               '--fullauto-mirror', 'off', '--fullauto-rotation-symmetry', 'off',
               '--fullauto-local-orientations', 'on', '--phi-points', '4',
               '--scattering-grid', '160', '170', '4', '8', '--fullauto-initial', '512',
               '--fullauto-pilot', '64', '--fullauto-kernel-count', '64',
               '--fullauto-min-correction', '64', '--fullauto-theta-start', '5',
               '--max-orientations', '2048', '--fullauto-max-rounds', '2',
               '--fullauto-final-audit', 'off', '--fullauto-final-checks', 'off',
               '--output', str(folder), '--close']
    (root/'interrupt.args.json').write_text(json.dumps(command, indent=2))
    with (root/'interrupt.stdout').open('w') as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
        started, interrupted = time.monotonic(), False
        while process.poll() is None:
            try:
                state = json.loads((folder/'fullauto_status.json').read_text())
                if state['round'] >= 1 and state['status'] == 'evaluating':
                    process.terminate()
                    interrupted = True
                    break
            except (FileNotFoundError, json.JSONDecodeError):
                pass
            if time.monotonic()-started > 200:
                process.terminate()
                process.wait(timeout=20)
                raise AssertionError('interruption test timed out')
            time.sleep(.05)
        process.wait(timeout=20)
    assert interrupted, 'test must interrupt a real unfinished round'
    completed = {str(p): digest(p) for p in folder.glob('jobs/*/complete.txt')}
    matrices = {str(p): digest(p) for p in folder.glob('round_snapshots/*/mueller_evaluated.dat')}
    assert completed and matrices and not reports(folder)
    run(command, root/'resume.stdout', ok=(0, 3, 4))
    assert all(digest(Path(p)) == value for p, value in completed.items())
    assert all(digest(Path(p)) == value for p, value in matrices.items())
    assert not reports(folder), reports(folder)
    status = json.loads((folder/'fullauto_status.json').read_text())
    assert status['angular_sampling_round'] >= 1
    for marker in folder.glob('jobs/j*/complete.txt'):
        recipe = marker.parent/'command.args'
        assert '--no-text-reports' in recipe.read_text().splitlines(), str(recipe)
        assert (marker.parent/'stdout.log').stat().st_size > 0
        assert (marker.parent/'integral_outputs.fnv').is_file()
    policy = json.loads((folder/'worker_output_policy.json').read_text())
    assert policy['checkpoint_graph_preserved'] and not policy['existing_files_deleted']
    return dict(retained_completed_blocks=len(completed), retained_snapshots=len(matrices),
                round=status['round'], duplicate_text_reports=0, persistent=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = args.binary.resolve()
    result = dict(binary_sha256=digest(binary))
    result['native'] = compare_native(binary, root)
    result['persistent_resume'] = interrupt_persistent(binary, root)
    result['passed'] = True
    (root/'verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
