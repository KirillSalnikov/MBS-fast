#!/usr/bin/env python3
"""Real CPU optical equivalence and recovery tests for fullauto working storage.

All calculations use one supplied binary. Fresh legacy/memory pairs perform
the same scientific work; timings are separate from cached restart timings.
This script never touches a production output or starts a remote job.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import signal
import subprocess
import time


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--timeout', type=float, default=600)
    parser.add_argument('--rounds', type=int, choices=[1, 2], default=2)
    parser.add_argument('--cases', default='po_coherent,po_incoherent,go')
    parser.add_argument('--tolerance', type=float, default=1e-11)
    parser.add_argument('--crash', action='store_true', help='also SIGKILL an isolated larger run after its first periodic checkpoint')
    return parser.parse_args()


def load(path):
    return json.loads(Path(path).read_text())


def table(path):
    lines = Path(path).read_text().splitlines()
    assert lines, str(path)
    delimiter = '\t' if '\t' in lines[0] else ','
    return list(csv.DictReader(lines, delimiter=delimiter))


def matrix(path):
    rows = [[float(value) for value in line.split()] for line in Path(path).read_text().splitlines()[1:] if line.strip()]
    assert rows and all(len(row) == 18 for row in rows), str(path)
    assert all(math.isfinite(value) for row in rows for value in row), str(path)
    return rows


def numeric_equal(left, right, tolerance):
    try:
        a, b = float(left), float(right)
    except (TypeError, ValueError):
        return left == right
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if not math.isfinite(a) or not math.isfinite(b):
        return a == b
    return abs(a - b) <= tolerance * max(1.0, abs(a), abs(b))


def fnv1a64(value):
    digest = 14695981039346656037
    for byte in value:
        digest = ((digest ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return digest


def verify_final_exports(output):
    manifest = output / 'final_outputs.fnv'
    rows = [line.split() for line in manifest.read_text().splitlines() if line.strip()]
    assert rows and all(len(row) == 2 for row in rows), str(manifest)
    for name, expected in rows:
        relative = Path(name)
        assert not relative.is_absolute() and '..' not in relative.parts, name
        path = output / relative
        assert path.is_file(), ('final manifest references a missing physical file', str(path))
        assert fnv1a64(path.read_bytes()) == int(expected, 16), ('final physical export hash mismatch', str(path))
    return {'authenticated_physical_files': len(rows), 'files': [row[0] for row in rows]}


def compare_table(left, right, label, tolerance):
    assert len(left) == len(right), (label, len(left), len(right))
    for index, (a, b) in enumerate(zip(left, right)):
        assert a.keys() == b.keys(), (label, index, a.keys(), b.keys())
        for key in a:
            # A repeated capped round has prior-history information that the
            # first round did not. It must not alter its sample budgets, matrix
            # or confidence intervals; these two scheduling flags may change.
            if label == 'parameters' and key in ['previous_available', 'half_prefix_stable']:
                continue
            assert numeric_equal(a[key], b[key], tolerance), (label, index, key, a[key], b[key])


def capture(output):
    current = load(output / 'intermediate_current.json')
    status = load(output / 'fullauto_status.json')
    folder = Path(current['snapshot_folder'])
    values = {'matrix': matrix(current['matrix']), 'parameters': table(current['accuracy_and_parameters']),
              'statistics': {key: value for key, value in current.items()
                             if key in ['statistical95', 'output_statistical95', 'interpolation_residual95',
                                        'combined_direct95', 'combined_output95', 'angular_sampling_noise95',
                                        'angular_sampling_failed_pairs', 'angular_shape_failed_stencils']},
              'target_depth': status['target_depth'], 'counts_per_seed': status['counts_per_seed'],
              'integrals': table(current['integrals']), 'integral_replicas': table(folder / 'integral_replicas.tsv'),
              'integral_confidence': table(folder / 'integral_confidence.tsv')}
    assert len(values['integral_replicas']) == 8, 'eight separate production seeds must be retained'
    assert len({row['seed'] for row in values['integral_replicas']}) == 8
    for name in ['orientation_budget.csv', 'interpolation_orientation_budget.csv']:
        path = folder / name
        if not path.exists():
            path = output / name
        if path.exists():
            values[name] = table(path)
    for name in ['mueller_forward.dat', 'mueller_backward.dat']:
        if (folder / name).exists():
            values[name] = matrix(folder / name)
    return values


def compare(a, b, tolerance, expected_depth=6):
    assert a.keys() == b.keys(), (a.keys(), b.keys())
    assert a['counts_per_seed'] == b['counts_per_seed'], ('different scientific work', a['counts_per_seed'], b['counts_per_seed'])
    assert a['target_depth'] == b['target_depth']
    if expected_depth is not None:
        assert a['target_depth'] == expected_depth
    assert a['statistics'].keys() == b['statistics'].keys()
    for key in a['statistics']:
        assert numeric_equal(a['statistics'][key], b['statistics'][key], tolerance), (key, a['statistics'][key], b['statistics'][key])
    worst = 0.0
    for key in ['matrix', 'mueller_forward.dat', 'mueller_backward.dat']:
        if key not in a:
            continue
        assert len(a[key]) == len(b[key])
        for x, y in zip(a[key], b[key]):
            assert all(numeric_equal(u, v, tolerance) for u, v in zip(x[:2], y[:2]))
            scale = max(abs(x[2]), abs(y[2]), 1e-300)
            error = max(abs(u - v) / scale for u, v in zip(x[2:], y[2:]))
            worst = max(worst, error)
            assert error <= tolerance, (key, x[0], error, tolerance)
    for key in a.keys() - {'matrix', 'mueller_forward.dat', 'mueller_backward.dat', 'statistics', 'target_depth', 'counts_per_seed'}:
        compare_table(a[key], b[key], key, tolerance)
    return {'max_all16_difference_scaled_M11': worst, 'counts_per_seed': a['counts_per_seed'],
            'target_depth': a['target_depth'], 'eight_integral_seeds_equal': True, 'confidence_tables_equal': True}


def proc_io(pid):
    try:
        return {key: int(value) for key, value in (line.split(':', 1) for line in Path(f'/proc/{pid}/io').read_text().splitlines())}
    except (OSError, ValueError):
        return {}


class Suite:
    def __init__(self, args):
        self.args = args
        self.root = args.output.resolve()
        assert not self.root.exists() or not any(self.root.iterdir()), 'Choose an empty evidence directory; existing data is preserved.'
        self.root.mkdir(parents=True, exist_ok=True)
        self.binary = args.binary.resolve()
        self.receipt = {'binary': str(self.binary), 'binary_sha256': hashlib.sha256(self.binary.read_bytes()).hexdigest(),
                        'runs': [], 'comparisons': {}, 'recovery': {}, 'production_mutations': False}
        self.save()

    def save(self):
        (self.root / 'receipt.json').write_text(json.dumps(self.receipt, indent=2) + '\n')

    def command(self, case, output, storage, checkpoint=False, large=False):
        go = case == 'go'
        command = [str(self.binary), '--method', 'go' if go else 'po', '--backend', 'cpu',
                   '--particle', '1', '5', '2.5', '--ri', '1.3116', '0', '--wavelength-um', '.532',
                   '--threads', '1', '--fullauto', '.04', '--max-reflections', '6',
                   '--fullauto-initial', '4096' if large else '64', '--max-orientations', '8192' if large else '64',
                   '--fullauto-max-rounds', '1' if large else str(self.args.rounds),
                   '--fullauto-resources', 'off', '--fullauto-cutoffs', 'off', '--cutoff-profile', 'off',
                   '--fullauto-analytic-controls', 'off', '--fullauto-mirror', 'off',
                   '--fullauto-rotation-symmetry', 'off', '--fullauto-local-orientations', 'off',
                   '--fullauto-final-checks', 'off', '--symmetry', '1', '1', '--fullauto-workers', 'persistent',
                   '--output', str(output), '--close']
        if go:
            command += ['--scattering-grid', '0', '180', '1', '12', '--phi-points', '1']
        else:
            command += ['--scattering-grid', '150', '170', '4', '6', '--phi-points', '4',
                        '--fullauto-theta-start', '7', '--fullauto-theta-audit-count', '0',
                        '--fullauto-pilot', '16', '--fullauto-min-correction', '16' if large else '64']
        if case == 'po_incoherent':
            command += ['--incoherent']
        if storage == 'legacy':
            command += ['--fullauto-storage', 'legacy', '--checkpoint']
        if checkpoint and storage != 'legacy':
            command += ['--fullauto-checkpoint-seconds', '1']
        return command

    def run(self, name, command, expected=(0, 3), kill_after_snapshot=None):
        log = self.root / (name + '.log')
        (self.root / (name + '.args.json')).write_text(json.dumps(command, indent=2) + '\n')
        before = resource.getrusage(resource.RUSAGE_CHILDREN)
        started = time.monotonic()
        maximum_io = {}
        killed = False
        with log.open('w') as stream:
            process = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                while process.poll() is None:
                    for key, value in proc_io(process.pid).items():
                        maximum_io[key] = max(value, maximum_io.get(key, 0))
                    if time.monotonic() - started > self.args.timeout:
                        raise TimeoutError(f'{name}: exceeded {self.args.timeout} seconds')
                    if kill_after_snapshot is not None and kill_after_snapshot.exists():
                        assert kill_after_snapshot.stat().st_size > 0
                        os.killpg(process.pid, signal.SIGKILL)
                        killed = True
                        process.wait(timeout=10)
                        break
                    time.sleep(.025)
                code = process.wait()
            finally:
                if process.poll() is None:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        record = {'name': name, 'returncode': code, 'wall_seconds': time.monotonic() - started,
                  'cpu_user_seconds': after.ru_utime - before.ru_utime, 'cpu_system_seconds': after.ru_stime - before.ru_stime,
                  'filesystem_input_blocks': after.ru_inblock - before.ru_inblock,
                  'filesystem_output_blocks': after.ru_oublock - before.ru_oublock,
                  'sampled_controller_io_lower_bound': maximum_io, 'killed_after_checkpoint': killed,
                  'io_scope': 'rusage child process tree; proc counters sampled on controller only and may miss its final writes'}
        self.receipt['runs'].append(record)
        self.save()
        if kill_after_snapshot is not None:
            assert killed and code == -signal.SIGKILL, (name, record, 'calculation finished before a live periodic checkpoint was observed')
        else:
            assert code in expected, (name, code, log.read_text()[-4000:])
        return record

    def check_recovery(self, case, mode, command, output, before):
        if mode == 'memory':
            assert (output / 'fullauto_checkpoint.bin').is_file()
        history = output / 'fullauto_history.jsonl'
        history_offset = history.stat().st_size
        self.run(f'{case}_{mode}_resume', command)
        after = capture(output)
        comparison = compare(before, after, self.args.tolerance)
        with history.open() as stream:
            stream.seek(history_offset)
            events = [json.loads(line) for line in stream if line.startswith('{')]
        restores = [event for event in events if event['event'] in ['cache_restored', 'cache_restore_reuse']]
        assert restores, (case, mode, 'resume did not report verified cache restoration')
        assert max(event['data'].get('indexed_jobs', 0) for event in restores) > 0
        comparison['cached_jobs_restored'] = max(event['data'].get('indexed_jobs', 0) for event in restores)
        self.receipt['recovery'][f'{case}_{mode}'] = comparison
        self.save()

    def check_completed_exports(self):
        captured = {}
        exports = {}
        for mode in ['legacy', 'memory']:
            output = self.root / ('completed_' + mode)
            command = self.command('po_coherent', output, mode, checkpoint=True)
            # Deliberately loose acceptance makes this a publication fixture;
            # it is not a claim that 128 orientations achieve production 4%.
            # Prefix stability requires a doubled initial prefix even at loose
            # epsilon; keep initial/minimum 64 and permit N128 in round two.
            command[command.index('--fullauto') + 1] = '.5'
            command[command.index('--max-orientations') + 1] = '128'
            command[command.index('--fullauto-max-rounds') + 1] = '2'
            self.run('completed_' + mode, command, expected=(0,))
            exports[mode] = verify_final_exports(output)
            captured[mode] = capture(output)
        self.receipt['comparisons']['completed_publication'] = compare(captured['legacy'], captured['memory'], self.args.tolerance)
        self.receipt['comparisons']['completed_publication']['scope'] = '50% acceptance used only to test final publication at identical N128; not an accuracy benchmark'
        self.receipt['final_exports'] = exports
        self.save()

    def check_automatic_depth(self):
        captures = {}
        for mode in ['legacy', 'memory']:
            output = self.root / ('automatic_depth_' + mode)
            command = self.command('po_coherent', output, mode, checkpoint=True)
            index = command.index('--max-reflections')
            del command[index:index + 2]
            for option, value in [('--fullauto-initial', '8'), ('--max-orientations', '128'),
                                  ('--fullauto-min-correction', '8'), ('--fullauto-pilot', '8')]:
                command[command.index(option) + 1] = value
            command += ['--fullauto-depth-start', '6', '--fullauto-depth-max', '8', '--fullauto-depth-pilot', '8']
            record = self.run('automatic_depth_' + mode, command, expected=(0, 3, 4))
            selection = load(output / 'depth_selection.json')
            for key in ['row_parameters', 'candidate_comparison', 'forecast_worker_seconds']:
                selection.pop(key, None)
            pilots = []
            for path in sorted(output.glob('depth_selection_candidate*.json')):
                pilot = load(path)
                pilot.pop('row_parameters', None)
                pilots.append(pilot)
            assert pilots, 'automatic depth fixture did not run paired pilots'
            candidates = table(output / 'depth_candidates.csv')
            for row in candidates:
                row.pop('forecast_worker_seconds', None)
            current = load(output / 'intermediate_current.json')
            captures[mode] = {'returncode': record['returncode'], 'selection': selection, 'pilots': pilots,
                              'candidates': candidates, 'matrix': matrix(current['matrix']),
                              'parameters': table(current['accuracy_and_parameters'])}
            (self.root / ('automatic_depth_' + mode + '_science.json')).write_text(json.dumps(captures[mode], indent=2) + '\n')
        left, right = captures['legacy'], captures['memory']
        assert left['returncode'] == right['returncode']
        for field in ['selection', 'pilots', 'candidates', 'parameters']:
            compare_table([left[field]] if field == 'selection' else left[field],
                          [right[field]] if field == 'selection' else right[field], field, self.args.tolerance)
        assert len(left['matrix']) == len(right['matrix'])
        worst = 0.0
        for a, b in zip(left['matrix'], right['matrix']):
            assert all(numeric_equal(x, y, self.args.tolerance) for x, y in zip(a[:2], b[:2]))
            scale = max(abs(a[2]), abs(b[2]), 1e-300)
            worst = max(worst, max(abs(x - y) / scale for x, y in zip(a[2:], b[2:])))
        assert worst <= self.args.tolerance, ('automatic depth latest matrix changed', worst)
        self.receipt['comparisons']['automatic_depth'] = {'returncode': left['returncode'], 'selection': left['selection'],
            'paired_pilot_receipts_equal': len(left['pilots']), 'max_all16_difference_scaled_M11': worst,
            'per_angle_decisions_and_confidence_equal': True}
        self.save()

    def execute(self):
        for case in self.args.cases.split(','):
            assert case in ['po_coherent', 'po_incoherent', 'go'], case
            captured = {}
            for mode in ['legacy', 'memory']:
                output = self.root / f'{case}_{mode}'
                command = self.command(case, output, mode, checkpoint=True)
                self.run(f'{case}_{mode}_fresh', command)
                captured[mode] = capture(output)
                (self.root / f'{case}_{mode}_science.json').write_text(json.dumps(captured[mode], indent=2) + '\n')
                if mode == 'memory':
                    assert load(output / 'fullauto_storage.json')['storage'] == 'memory'
                    assert not (output / 'jobs').exists() or not any((output / 'jobs').iterdir()), 'native jobs leaked to physical output'
                self.check_recovery(case, mode, command, output, captured[mode])
            self.receipt['comparisons'][case] = compare(captured['legacy'], captured['memory'], self.args.tolerance)
            runs = {record['name']: record for record in self.receipt['runs']}
            old, new = runs[f'{case}_legacy_fresh'], runs[f'{case}_memory_fresh']
            self.receipt['comparisons'][case]['fresh_wall_speedup'] = old['wall_seconds'] / new['wall_seconds']
            self.save()

        output = self.root / 'no_checkpoint'
        command = self.command('po_coherent', output, 'memory')
        self.run('no_checkpoint_fresh', command)
        assert (output / 'memory_storage_v1.txt').is_file() and not (output / 'fullauto_checkpoint.bin').exists()
        self.run('no_checkpoint_resume_rejected', command, expected=(1, 2))
        assert 'no recovery checkpoint' in (self.root / 'no_checkpoint_resume_rejected.log').read_text()
        self.receipt['recovery']['missing_snapshot_rejected'] = True

        first = self.args.cases.split(',')[0]
        old = self.root / f'{first}_legacy'
        command = self.command(first, old, 'legacy')
        command[command.index('--fullauto-storage') + 1] = 'memory'
        self.run('legacy_to_memory_rejected', command, expected=(1, 2))
        assert 'existing legacy checkpoints are preserved' in (self.root / 'legacy_to_memory_rejected.log').read_text()
        self.receipt['recovery']['unsafe_legacy_conversion_rejected'] = True
        self.save()

        self.check_completed_exports()
        self.check_automatic_depth()

        if self.args.crash:
            interrupted = self.root / 'crash_memory'
            command = self.command('po_coherent', interrupted, 'memory', checkpoint=True, large=True)
            self.run('crash_memory_killed', command, kill_after_snapshot=interrupted / 'fullauto_checkpoint.bin')
            self.run('crash_memory_resumed', command)
            reference = self.root / 'crash_reference'
            baseline = self.command('po_coherent', reference, 'memory', checkpoint=True, large=True)
            self.run('crash_memory_uninterrupted', baseline)
            self.receipt['recovery']['killed_vs_uninterrupted'] = compare(capture(interrupted), capture(reference), self.args.tolerance)
        self.receipt['passed'] = True
        self.save()
        print(json.dumps({'passed': True, 'receipt': str(self.root / 'receipt.json'), 'comparisons': self.receipt['comparisons']}, indent=2))


if __name__ == '__main__':
    suite = Suite(arguments())
    try:
        suite.execute()
    except Exception as error:
        suite.receipt['passed'] = False
        suite.receipt['failure'] = repr(error)
        suite.save()
        raise
