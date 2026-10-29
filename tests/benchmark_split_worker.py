#!/usr/bin/env python3
"""Compare split tracing in persistent workers on an unchanged droxtal sample.

Only the test output directory is written. One request runs at a time. Warm-up
is excluded; five matched independent Sobol seeds are compared in alternating
order. Raw/prepared result caches are disabled to measure actual tracing.
"""
import argparse
import json
import os
import re
import socket
import statistics
import struct
import subprocess
import time
from pathlib import Path
from test_droxtal_resident import load, difference, integral_difference


class Worker:
    def __init__(self, binary, schedule, root):
        self.socket, child = socket.socketpair()
        self.socket.settimeout(300)
        env = os.environ.copy()
        env.update(MBS_GPU_TRACE_CONVEX_SCHEDULE=schedule,
                   MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB='0',
                   MBS_HOST_MEM_BUDGET_MB='2048', OMP_WAIT_POLICY='PASSIVE')
        self.binary = str(binary.resolve())
        (root / 'environment.json').write_text(json.dumps({
            k: v for k, v in env.items() if k.startswith(('MBS_', 'CUDA_', 'OMP_'))}, indent=2))
        with (root / 'service.log').open('w') as log:
            self.process = subprocess.Popen([self.binary, '--fullauto-worker-server', str(child.fileno())],
                pass_fds=[child.fileno()], stdout=log, stderr=log, env=env)
        child.close()

    def run(self, command, root):
        root.mkdir()
        log = root / 'stdout.log'
        args = command + ['--output', str(root / 'result')]
        (root / 'command.args').write_text('\n'.join([self.binary] + args) + '\n')
        words = [str(log), '0', 'auto'] + args
        start = time.monotonic()
        self.socket.sendall(struct.pack('=I', len(words)) + b''.join(
            struct.pack('=I', len(s.encode())) + s.encode() for s in words))
        reply = b''
        while len(reply) < 4:
            part = self.socket.recv(4 - len(reply))
            assert part, ('worker exited', log)
            reply += part
        elapsed = time.monotonic() - start
        assert struct.unpack('=i', reply)[0] == 0, log
        phases = {}
        for line in (root / 'result/result_phase_timings.tsv').read_text().splitlines()[1:]:
            name, value, _ = line.split()
            phases[name] = float(value)
        text = log.read_text()
        accepted = sum(map(int, re.findall(r'GPU convex resident direct: accepted=(\d+)', text)))
        assert not re.search(r'CPU_fallback=[1-9]', text), log
        assert 'GPU raw trace reuse: hit' not in text, log
        return dict(path=str(root), wall=elapsed, phase=phases, accepted=accepted)

    def close(self):
        try:
            if self.process.poll() is None:
                self.socket.sendall(struct.pack('=I', 0))
                self.process.wait(timeout=20)
        except (OSError, subprocess.TimeoutExpired):
            pass
        finally:
            self.socket.close()
            if self.process.poll() is None:
                self.process.terminate()
                self.process.wait(timeout=20)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--before', type=Path, required=True)
    parser.add_argument('--recipe', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--tolerance', type=float, default=1e-7)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True)
    recipe = json.loads(args.recipe.read_text())
    command = recipe['args'][1:]
    for flag in ['--output', '--po-trace-backend', '--orientation-chunk', '--trace-cache-mib']:
        while flag in command:
            at = command.index(flag)
            del command[at:at + 2]
    grid = root / 'theta.txt'
    grid.write_text(recipe['theta'])
    command[command.index('--theta-grid-file') + 1] = str(grid)
    command[command.index('--sobol-seed') + 1] = '1024'
    command[command.index('--threads') + 1] = '2'
    command += ['--orientation-chunk', '32', '--trace-cache-mib', '0',
                '--profile-phases', '--allow-experimental-environment']
    records = []
    for backend in ['resident', 'hybrid']:
        workers = {}
        try:
            for mode, binary, schedule in [('before', args.before, 'legacy'), ('split', args.binary, 'split8')]:
                dest = root / (backend + '_' + mode)
                dest.mkdir()
                workers[mode] = Worker(binary, schedule, dest)
            for repeat, seed in enumerate([23, 131, 107, 89, 71, 53]):
                cmd = list(command)
                cmd[cmd.index('--sobol-seed') + 2] = str(seed)
                cmd += ['--po-trace-backend', backend]
                runs = {}
                for mode in (['before', 'split'] if repeat % 2 == 0 else ['split', 'before']):
                    runs[mode] = workers[mode].run(cmd, root / (backend + '_' + mode) / ('seed' + str(seed)))
                    accepted = runs[mode]['accepted']
                    assert (accepted == 1024 if backend == 'resident' else accepted > 0), runs[mode]
                a, b = [Path(runs[mode]['path']) / 'result' for mode in ['before', 'split']]
                errors = {f.name: difference(load(f), load(b / f.name)) for f in a.glob('result*.dat')}
                errors.update({f.name: integral_difference(f, b / f.name) for f in a.glob('result*_integrals.tsv')})
                assert errors and max(errors.values()) <= args.tolerance, errors
                row = dict(backend=backend, seed=seed, warmup=repeat == 0, runs=runs, errors=errors)
                records.append(row)
                (root / 'comparison.json').write_text(json.dumps(records, indent=2))
                print(json.dumps(dict(backend=backend, seed=seed, wall={k: v['wall'] for k, v in runs.items()},
                                      processing={k: v['phase']['processing'] for k, v in runs.items()})), flush=True)
        finally:
            for worker in workers.values():
                worker.close()
    summary = {}
    for backend in ['resident', 'hybrid']:
        selected = [r for r in records if r['backend'] == backend and not r['warmup']]
        summary[backend] = {mode: dict(
            processing=statistics.median(r['runs'][mode]['phase']['processing'] for r in selected),
            wall=statistics.median(r['runs'][mode]['wall'] for r in selected)) for mode in ['before', 'split']}
    (root / 'summary.json').write_text(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
