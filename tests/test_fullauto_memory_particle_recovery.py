#!/usr/bin/env python3
"""Real fullauto checkpoint identity checks for particle visibility sidecars.

All fixtures are created beneath a new evidence directory. The supplied binary
is never rebuilt, production output is never opened, and no remote job runs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from test_fullauto_memory_runtime import capture, compare

TOLERANCE = 1e-14  # saved decimal intermediates may differ by final-bit roundoff


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            digest.update(block)
    return digest.hexdigest()


def run(binary, particle, output, evidence, label, receipt, expected):
    command = [str(binary), '--method', 'po', '--backend', 'cpu',
               '--particle-file', str(particle), '--resize-dmax-um', '5',
               '--ri', '1.3116', '0', '--wavelength-um', '.532',
               '--threads', '1', '--fullauto', '.04', '--max-reflections', '6',
               '--fullauto-initial', '64', '--max-orientations', '64',
               '--fullauto-max-rounds', '1', '--fullauto-resources', 'off',
               '--fullauto-cutoffs', 'off', '--cutoff-profile', 'off',
               '--fullauto-analytic-controls', 'off', '--fullauto-mirror', 'off',
               '--fullauto-rotation-symmetry', 'off', '--fullauto-local-orientations', 'off',
               '--fullauto-final-checks', 'off', '--symmetry', '1', '1',
               '--fullauto-workers', 'persistent', '--fullauto-checkpoint-seconds', '1',
               '--scattering-grid', '150', '170', '4', '6', '--phi-points', '4',
               '--fullauto-theta-start', '7', '--fullauto-theta-audit-count', '0',
               '--fullauto-pilot', '16', '--fullauto-min-correction', '64',
               '--output', str(output), '--close']
    (evidence / (label + '.args.json')).write_text(json.dumps(command, indent=2) + '\n')
    started = time.monotonic()
    with (evidence / (label + '.log')).open('w') as log:
        process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=90)
    record = {'label': label, 'returncode': process.returncode, 'seconds': time.monotonic() - started}
    receipt['runs'].append(record)
    log = (evidence / (label + '.log')).read_text()
    assert process.returncode in expected, (record, log[-5000:])
    return log


def execute(args, receipt):
    evidence = args.output.resolve()
    binary = args.binary.resolve()
    geometry = args.particle.read_bytes()
    visible = b'6\n1 1\n1 1\n1 1\n1 1\n1 1\n1 1\n'
    changed = b'6\n0 1\n1 1\n1 1\n1 1\n1 1\n1 1\n'
    for inside in [False, True]:
        for present in [False, True]:
            name = ('inside' if inside else 'outside') + ('_present' if present else '_absent')
            home = evidence / name
            output = home / 'output'
            output.mkdir(parents=True)
            particle = (output if inside else home) / 'particle.dat'
            sidecar = Path(str(particle) + '.visibility')
            particle.write_bytes(geometry)
            if present:
                sidecar.write_bytes(visible)
            run(binary, particle, output, evidence, name + '_fresh', receipt, [0, 3])
            checkpoint = output / 'fullauto_checkpoint.bin'
            assert checkpoint.is_file() and checkpoint.stat().st_size > 0
            before = capture(output)
            run(binary, particle, output, evidence, name + '_unchanged_resume', receipt, [0, 3])
            unchanged = compare(before, capture(output), TOLERANCE)

            def reject(label, mutate, restore, expected_error):
                previous_checkpoint = sha(checkpoint)
                history = output / 'fullauto_history.jsonl'
                history_size = history.stat().st_size
                mutate()
                try:
                    log = run(binary, particle, output, evidence, name + '_' + label, receipt, [2])
                    assert expected_error in log, (name, label, log)
                    assert sha(checkpoint) == previous_checkpoint, 'rejected identity overwrote valid recovery state'
                    with history.open() as stream:
                        stream.seek(history_size)
                        events = [json.loads(line) for line in stream if line.startswith('{')]
                    assert not any(event['event'] in ['worker_start', 'worker_complete'] for event in events)
                    assert 'Prepared trace reuse:' not in log, 'rejected identity reached native optics'
                    compare(before, capture(output), TOLERANCE)
                finally:
                    restore()

            rejected = []
            if present:
                reject('changed_visibility', lambda: sidecar.write_bytes(changed),
                       lambda: sidecar.write_bytes(visible), 'memory checkpoint inputs/binary changed')
                rejected.append('changed_visibility')
                reject('deleted_visibility', sidecar.unlink, lambda: sidecar.write_bytes(visible),
                       'memory checkpoint inputs/binary changed')
                rejected.append('deleted_visibility')
            else:
                reject('added_visibility', lambda: sidecar.write_bytes(visible), sidecar.unlink,
                       'memory checkpoint inputs/binary changed')
                rejected.append('added_visibility')
            # Fullauto constructs the particle before loading its snapshot.
            # An owned particle stored in the snapshot is still an input, not
            # a promise that the snapshot can replace the physical source.
            reject('missing_particle', particle.unlink, lambda: particle.write_bytes(geometry), 'particle')
            rejected.append('missing_particle')
            run(binary, particle, output, evidence, name + '_restored_resume', receipt, [0, 3])
            restored = compare(before, capture(output), TOLERANCE)
            receipt['cases'][name] = {'particle_inside_output': inside, 'visibility_present': present,
                                      'unchanged_resume': unchanged, 'restored_resume': restored,
                                      'rejected_before_native_work': rejected,
                                      'physical_particle_still_required': True}
            (evidence / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
            print(name + ': sidecar identity rejection and all16/CI recovery PASS', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--particle', type=Path, default=Path(__file__).resolve().parents[1] / 'examples/cube.particle')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists() or not any(args.output.iterdir()), 'Use a new empty evidence directory.'
    args.output.mkdir(parents=True, exist_ok=True)
    receipt = {'binary': str(args.binary.resolve()), 'binary_sha256': sha(args.binary),
               'fixture': str(args.particle.resolve()), 'runs': [], 'cases': {},
               'production_mutations': False, 'comparison_tolerance': TOLERANCE, 'passed': False}
    try:
        execute(args, receipt)
        receipt['passed'] = True
    except Exception as error:
        receipt['failure'] = repr(error)
        raise
    finally:
        (args.output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('PASS: ' + str(args.output / 'receipt.json'))


if __name__ == '__main__':
    main()
