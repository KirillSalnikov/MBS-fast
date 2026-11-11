#!/usr/bin/env python3
"""Bounded, isolated 1/2-worker CUDA experiment; never restarts production.

Only execute when the experiment coordinator has reserved this measurement
window. The script also requires a quiet, stable production load before each
run and stops only its own process group if that condition changes.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--gpu', required=True, help='physical nvidia-smi index or UUID')
    p.add_argument('--cpus', required=True, help='exactly eight allowed distinct physical CPU cores, comma-separated')
    p.add_argument('--production-pids', default='', help='known production CUDA worker PIDs; all must be visible to NVML or the test refuses to start')
    p.add_argument('--timeout', type=float, default=120)
    p.add_argument('--quiet-seconds', type=float, default=8)
    p.add_argument('--maximum-production-utilization', type=float, default=10)
    p.add_argument('--count', type=int, default=8192)
    p.add_argument('--length', type=float, default=5)
    p.add_argument('--diameter', type=float, default=2.5)
    p.add_argument('--incoherent', action='store_true')
    p.add_argument('--tolerance', type=float, default=1e-6)
    return p.parse_args()


def fnv(value):
    digest = 14695981039346656037
    for byte in value:
        digest = ((digest ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return format(digest, 'x')


def dump(path, data):
    Path(path).write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')


def command_output(command):
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=8)
    if result.returncode:
        raise RuntimeError((command, result.returncode, result.stderr[-2000:]))
    return result.stdout


def number(value):
    try:
        return float(value.strip())
    except ValueError:
        return None


def gpu_snapshot(gpu):
    row = next(csv.reader(command_output(['nvidia-smi', '-i', gpu,
        '--query-gpu=index,uuid,utilization.gpu,memory.used,memory.free,power.draw',
        '--format=csv,noheader,nounits']).splitlines()))
    processes = {}
    text = command_output(['nvidia-smi', '-i', gpu,
        '--query-compute-apps=pid,used_gpu_memory', '--format=csv,noheader,nounits'])
    for values in csv.reader(text.splitlines()):
        if not values:
            continue
        if len(values) != 2 or not values[0].strip().isdigit():
            raise RuntimeError('GPU process enumeration is unavailable; cannot prove production isolation')
        processes[int(values[0])] = number(values[1])
    return {'index': int(row[0]), 'uuid': row[1].strip(), 'utilization': number(row[2]),
            'memory_used_mib': number(row[3]), 'memory_free_mib': number(row[4]),
            'power_watts': number(row[5]), 'processes': processes}


def process_utilization(gpu):
    text = command_output(['nvidia-smi', 'pmon', '-i', str(gpu), '-s', 'u', '-c', '1'])
    values = {}
    for line in text.splitlines():
        row = line.split()
        if len(row) >= 5 and row[0].isdigit() and row[1].isdigit():
            values[int(row[1])] = number(row[3])
    return values


def belongs(pid, group):
    try:
        return os.getpgid(pid) == group
    except ProcessLookupError:
        return False


def own_processes(pid):
    found = set()
    pending = [pid]
    while pending:
        child = pending.pop()
        if child in found:
            continue
        found.add(child)
        try:
            pending.extend(int(item) for item in Path(f'/proc/{child}/task/{child}/children').read_text().split())
        except (FileNotFoundError, ProcessLookupError):
            pass
    values = []
    for child in sorted(found):
        try:
            row = Path(f'/proc/{child}/stat').read_text().rsplit(')', 1)[1].split()
            io = {key: int(value) for key, value in (line.split(':', 1) for line in Path(f'/proc/{child}/io').read_text().splitlines())}
            values.append({'pid': child, 'cpu_ticks': int(row[11]) + int(row[12]),
                           'rss_bytes': int(row[21]) * os.sysconf('SC_PAGE_SIZE'), 'io': io})
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            pass
    return values


def stop_own(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=3)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait(timeout=5)
    finally:
        # The controller may exit before its children. Its private process
        # group remains ours until the last member has gone.
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def table(path):
    lines = Path(path).read_text().splitlines()
    return list(csv.DictReader(lines, delimiter='\t' if '\t' in lines[0] else ','))


def science(root):
    current = json.loads((root / 'intermediate_current.json').read_text())
    status = json.loads((root / 'fullauto_status.json').read_text())
    folder = Path(current['snapshot_folder'])
    matrix = [[float(x) for x in line.split()] for line in Path(current['matrix']).read_text().splitlines()[1:] if line.strip()]
    assert matrix and all(len(row) == 18 for row in matrix)
    assert all(math.isfinite(value) for row in matrix for value in row)
    return {'counts_per_seed': status['counts_per_seed'], 'target_depth': status['target_depth'], 'matrix': matrix,
            'parameters': table(current['accuracy_and_parameters']), 'integrals': table(current['integrals']),
            'integral_replicas': table(folder / 'integral_replicas.tsv'),
            'integral_confidence': table(folder / 'integral_confidence.tsv')}


def equal_value(a, b, tolerance):
    try:
        a, b = float(a), float(b)
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        if not math.isfinite(a) or not math.isfinite(b):
            return a == b
        return abs(a - b) <= tolerance * max(1., abs(a), abs(b))
    except (ValueError, TypeError):
        return a == b


def compare(left, right, tolerance):
    assert left['counts_per_seed'] == right['counts_per_seed']
    assert left['target_depth'] == right['target_depth'] == 6
    worst = 0.
    assert len(left['matrix']) == len(right['matrix'])
    for a, b in zip(left['matrix'], right['matrix']):
        assert all(equal_value(x, y, 1e-12) for x, y in zip(a[:2], b[:2]))
        scale = max(abs(a[2]), abs(b[2]), 1e-300)
        worst = max(worst, max(abs(x - y) / scale for x, y in zip(a[2:], b[2:])))
    assert worst <= tolerance, ('matrix difference / M11', worst)
    for key in ['parameters', 'integrals', 'integral_replicas', 'integral_confidence']:
        assert len(left[key]) == len(right[key])
        for a, b in zip(left[key], right[key]):
            assert a.keys() == b.keys()
            for column in a:
                assert equal_value(a[column], b[column], tolerance), (key, column, a[column], b[column])
    assert len(left['integral_replicas']) == 8
    return {'maximum_all16_difference_scaled_M11': worst, 'all_eight_integral_replicas_and_ci_equal': True}


class Experiment:
    def __init__(self, args):
        self.args = args
        self.known_production = {int(value) for value in args.production_pids.split(',') if value}
        if not 0 < args.timeout <= 120 or args.quiet_seconds < 5:
            raise ValueError('timeout must be <=120 seconds and quiet window >=5 seconds')
        cores = [int(value) for value in args.cpus.split(',')]
        allowed = os.sched_getaffinity(0)
        physical = set()
        for core in cores:
            base = Path(f'/sys/devices/system/cpu/cpu{core}/topology')
            physical.add(((base / 'physical_package_id').read_text().strip(), (base / 'core_id').read_text().strip()))
        if len(cores) != 8 or len(physical) != 8 or not set(cores) <= allowed:
            raise ValueError('choose exactly eight allowed distinct physical cores')
        if args.count < 16 or args.count & (args.count - 1):
            raise ValueError('count must be a power of two >=16')
        self.root = args.output.absolute()
        if self.root.exists():
            raise ValueError('output must be a new directory')
        self.root.mkdir(parents=True)
        self.active = self.root / 'active'
        self.baseline = self.root / 'baseline'
        (self.root / 'bin').mkdir()
        self.binary = self.root / 'bin' / 'mbs_gpu_benchmark'
        shutil.copy2(args.binary.resolve(), self.binary)
        payload = self.binary.read_bytes()
        identity = {'type': 'same immutable executable compared with itself',
                    'sha256': hashlib.sha256(payload).hexdigest(), 'fnv1a64': fnv(payload),
                    'claim': 'byte identity only; no claim about a different build or production compatibility'}
        proof = self.root / 'binary_self_identity.json'
        dump(proof, identity)
        receipt = f'FULLAUTO_NATIVE_COMPATIBILITY_V1\n{fnv(payload)} {fnv(payload)} {proof} {fnv(proof.read_bytes())}\n'
        Path(str(self.binary) + '.native-compatibility.tsv').write_text(fnv(receipt.encode()) + '\n' + receipt)
        self.env = dict(os.environ)
        for name in list(self.env):
            if name.startswith('MBS_'):
                del self.env[name]
        self.env.update(CUDA_VISIBLE_DEVICES=args.gpu, MBS_HOST_MEM_BUDGET_MB='1024', MBS_GPU_MEM_FRACTION='.20')
        self.receipt = {'binary': identity, 'production_mutations': False, 'cpu_affinity': cores,
                        'per_worker_host_budget_mib': 1024, 'per_worker_gpu_memory_fraction': .20,
                        'memory_comparison_scope': 'identical per-worker caps; aggregate cap increases with two workers; actual summed RSS recorded',
                        'runs': [], 'comparisons': []}
        self.save()

    def save(self):
        dump(self.root / 'receipt.json', self.receipt)

    def quiet(self, label):
        samples = []
        expected = None
        until = time.monotonic() + self.args.quiet_seconds
        while time.monotonic() < until:
            item = gpu_snapshot(self.args.gpu)
            samples.append(item)
            pids = set(item['processes'])
            if not self.known_production <= pids:
                raise RuntimeError('known production GPU PIDs are not visible; cannot prove production isolation')
            if expected is None:
                expected = pids
            if pids != expected or item['utilization'] is None or item['utilization'] > self.args.maximum_production_utilization:
                dump(self.root / (label + '_preflight.json'), samples)
                raise RuntimeError('production is not quiet/stable; benchmark was not started')
            time.sleep(.5)
        dump(self.root / (label + '_preflight.json'), samples)
        if samples[-1]['memory_free_mib'] < 2048:
            raise RuntimeError('less than 2 GiB free GPU memory; benchmark was not started')
        return samples[-1]

    def command(self):
        n = str(self.args.count)
        values = ['taskset', '--cpu-list', self.args.cpus, str(self.binary), '--method', 'po', '--backend', 'cuda',
            '--particle', '1', str(self.args.length), str(self.args.diameter), '--ri', '1.3116', '0', '--wavelength-um', '.532',
            '--threads', '8', '--fullauto', '.04', '--max-reflections', '6', '--phi-points', '16',
            '--scattering-grid', '150', '170', '16', '128', '--fullauto-theta-start', '129', '--fullauto-theta-audit-count', '0',
            '--fullauto-initial', n, '--fullauto-min-correction', n, '--max-orientations', n, '--fullauto-max-rounds', '1',
            '--fullauto-pilot', '16', '--fullauto-kernel-count', '16', '--fullauto-resources', 'off',
            '--fullauto-cutoffs', 'off', '--cutoff-profile', 'off', '--fullauto-analytic-controls', 'off',
            '--fullauto-mirror', 'off', '--fullauto-rotation-symmetry', 'off', '--fullauto-local-orientations', 'off',
            '--fullauto-final-checks', 'off', '--symmetry', '1', '1', '--fullauto-workers', 'persistent',
            '--fullauto-gpus', '0', '--fullauto-checkpoint-seconds', '1', '--fullauto-memory-mib', '512',
            '--output', str(self.active), '--close']
        if self.args.incoherent:
            values.append('--incoherent')
        return values

    def run(self, label, command, crash=False):
        before = self.quiet(label)
        dump(self.root / (label + '_args.json'), command)
        started = time.monotonic()
        record = {'name': label, 'preexisting_gpu_pids': sorted(before['processes']), 'samples': [], 'killed_for_checkpoint': False}
        with (self.root / (label + '.log')).open('w') as log:
            child = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, env=self.env, start_new_session=True)
            record['controller_pid'] = child.pid
            try:
                while child.poll() is None:
                    if time.monotonic() - started > self.args.timeout:
                        raise TimeoutError('case exceeded its bounded timeout')
                    item = gpu_snapshot(self.args.gpu)
                    other = {pid for pid in item['processes'] if not belongs(pid, child.pid)}
                    if other != set(before['processes']):
                        raise RuntimeError('production process set changed; comparison invalid')
                    utilization = process_utilization(item['index']) if other else {}
                    if any(utilization.get(pid) is None or utilization[pid] > self.args.maximum_production_utilization for pid in other):
                        raise RuntimeError('production per-process utilization changed/unavailable; comparison invalid')
                    for pid in other:
                        old, new = before['processes'][pid], item['processes'][pid]
                        if old is None or new is None or abs(new - old) > max(128., old * .1):
                            raise RuntimeError('production GPU memory changed; comparison invalid')
                    processes = own_processes(child.pid)
                    item.update(elapsed_seconds=time.monotonic() - started, own_processes=processes,
                                own_summed_rss_bytes=sum(value['rss_bytes'] for value in processes), other_sm_utilization=utilization)
                    record['samples'].append(item)
                    checkpoint = self.active / 'fullauto_checkpoint.bin'
                    if crash and checkpoint.exists() and checkpoint.stat().st_size:
                        os.killpg(child.pid, signal.SIGKILL)
                        child.wait(timeout=5)
                        record['killed_for_checkpoint'] = True
                        break
                    time.sleep(.2)
                code = child.wait()
                if crash:
                    assert code == -signal.SIGKILL and record['killed_for_checkpoint'], 'fixture finished before a partial periodic checkpoint'
                else:
                    assert code in (0, 3), (label, code)
                record['returncode'] = code
            except Exception as error:
                record['failure'] = repr(error)
                raise
            finally:
                stop_own(child)
                record['wall_seconds'] = time.monotonic() - started
                record['peak_summed_rss_bytes'] = max((sample['own_summed_rss_bytes'] for sample in record['samples']), default=0)
                dump(self.root / (label + '_metrics.json'), record)
                self.receipt['runs'].append({key: value for key, value in record.items() if key != 'samples'})
                self.save()
        return record

    def execute(self):
        self.run('seed_checkpoint', self.command(), crash=True)
        shutil.copytree(self.active, self.baseline)
        baseline_digest = hashlib.sha256((self.baseline / 'fullauto_checkpoint.bin').read_bytes()).hexdigest()
        reference_science = reference_work = None
        for label, slots in [('A1', 1), ('B1', 2), ('B2', 2), ('A2', 1)]:
            # These exact directories were created by this instance. The saved
            # root path must remain unchanged inside the authenticated snapshot.
            assert self.active.parent == self.root and not self.active.is_symlink()
            shutil.rmtree(self.active)
            shutil.copytree(self.baseline, self.active)
            assert hashlib.sha256((self.active / 'fullauto_checkpoint.bin').read_bytes()).hexdigest() == baseline_digest
            command = self.command()
            command[command.index('--fullauto-checkpoint-seconds') + 1] = '0'
            command += ['--fullauto-resume-compatible', str(self.binary), '--fullauto-dispatch-gpus', '0' if slots == 1 else '0,0',
                        '--fullauto-dispatch-threads', '8' if slots == 1 else '4']
            record = self.run(label, command)
            data = science(self.active)
            assert all(value == self.args.count for value in data['counts_per_seed'])
            history = [json.loads(line) for line in (self.active / 'fullauto_history.jsonl').read_text().splitlines() if line.startswith('{')]
            current = [event for event in history if event.get('controller_pid') == record['controller_pid']]
            starts = [event['data'] for event in current if event.get('event') == 'worker_start']
            fields = ['count', 'offset', 'seed', 'depth', 'paired_depth', 'phi', 'theta_rows']
            work = sorted((Path(item['folder']).name,) + tuple(item.get(key, 0) for key in fields) for item in starts)
            assert work, 'resume performed no new native work'
            assert any(event['event'] in ('cache_restored', 'cache_restore_reuse') for event in current), 'no checkpoint reuse receipt'
            data['remaining_native_work'] = work
            dump(self.root / (label + '_science.json'), data)
            progress = json.loads((self.active / 'execution_progress.json').read_text())
            dump(self.root / (label + '_execution_progress.json'), progress)
            if reference_science is None:
                reference_science, reference_work = data, work
            else:
                assert work == reference_work, 'different remaining native work; timing comparison invalid'
                result = compare(reference_science, data, self.args.tolerance)
                self.receipt['comparisons'].append({'reference': 'A1', 'candidate': label, **result})
            self.save()
        values = {run['name']: run['wall_seconds'] for run in self.receipt['runs']}
        self.receipt['same_remaining_work_wall_speedup'] = ((values['A1'] + values['A2']) / (values['B1'] + values['B2']))
        self.receipt['passed'] = True
        self.save()
        print(json.dumps({'passed': True, 'receipt': str(self.root / 'receipt.json'),
                          'same_remaining_work_wall_speedup': self.receipt['same_remaining_work_wall_speedup']}, indent=2))


if __name__ == '__main__':
    experiment = Experiment(arguments())
    try:
        experiment.execute()
    except Exception as error:
        experiment.receipt['passed'] = False
        experiment.receipt['failure'] = repr(error)
        experiment.save()
        raise
