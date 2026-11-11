#!/usr/bin/env python3
"""Native CPU optics equality: legacy, V1, streaming V2 RAM and disk stages."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import socket
import struct
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', required=True, type=Path)
p.add_argument('--output', required=True, type=Path)
a = p.parse_args()
root = a.output.resolve()
root.mkdir(parents=True, exist_ok=False)
binary = a.binary.resolve()


def digest(data):
    value = 14695981039346656037
    for byte in data:
        value = ((value ^ byte) * 1099511628211) & ((1 << 64) - 1)
    return value


def send_words(channel, words):
    channel.sendall(struct.pack('<I', len(words)))
    for word in words:
        data = word.encode()
        channel.sendall(struct.pack('<I', len(data)) + data)


def receive(channel, count):
    data = bytearray()
    while len(data) < count:
        piece = channel.recv(min(count - len(data), 65536))
        if not piece:
            raise RuntimeError('worker socket EOF')
        data.extend(piece)
    return bytes(data)


def send_files(channel, files, version):
    total = sum(len(name.encode()) + len(data) for name, data in files.items())
    channel.sendall(struct.pack('<IIIQ', 0x4d425342 if version == 1 else 0x4d425353, version, len(files), total))
    for name, data in files.items():
        path = name.encode()
        if version == 1:
            channel.sendall(struct.pack('<IQQ', len(path), len(data), digest(path + b'\0' + data)) + path)
        else:
            expected = digest(data)
            descriptor = digest(path + b'\0' + struct.pack('<QQ', len(data), expected))
            channel.sendall(struct.pack('<IQQQ', len(path), len(data), expected, descriptor) + path)
        for offset in range(0, len(data), 65536):
            channel.sendall(data[offset:offset + 65536])


def receive_files(channel, version):
    magic, actual, count, expected = struct.unpack('<IIIQ', receive(channel, 20))
    assert magic == (0x4d425342 if version == 1 else 0x4d425353) and actual == version and count <= 4096
    assert expected <= (64 * 1024 * 1024 if version == 1 else 8 * 1024**3)
    files = {}
    total = 0
    for _ in range(count):
        header = struct.unpack('<IQQ' if version == 1 else '<IQQQ', receive(channel, 20 if version == 1 else 28))
        length, size, expected_hash = header[:3]
        assert 0 < length <= 4096 and size + length <= expected - total
        path = receive(channel, length)
        if version == 2:
            assert digest(path + b'\0' + struct.pack('<QQ', size, expected_hash)) == header[3]
        data = receive(channel, size)
        assert digest(path + b'\0' + data if version == 1 else data) == expected_hash
        name = path.decode()
        assert not Path(name).is_absolute() and '..' not in Path(name).parts and name not in files
        files[name] = data
        total += length + size
    assert total == expected
    return files


def command(logical, index, incoherent, controls):
    args = ['--method', 'po', '--backend', 'cpu', '--particle', '1', '5', '2.5', '--wavelength-um', '.532',
            '--ri', '1.3116', '0', '--max-reflections', '6', '--sobol-seed', '64', '23' if controls else '11',
            '--haar-alpha', '--phi-points', '4', '--theta-grid-file', logical + '/grids/theta.txt',
            '--threads', '2', '--close', '--output', logical + f'/jobs/j{index}/result']
    if incoherent:
        args.append('--incoherent')
    if controls:
        args += ['--paired-depth', '4', '--cutoff-profile', 'balanced', '--analytic-facet-average',
                 '--analytic-shadow-control', 'circular', '--analytic-mean-cache', logical + '/means/shared.cache',
                 '--analytic-mean-reference-grid', logical + '/grids/theta.txt']
    return [logical + f'/jobs/j{index}/stdout.log', '64', 'auto'] + args


receipt = {'binary': str(binary), 'binary_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
           'runs': {}, 'matrix_comparisons': [], 'integral_comparisons': [], 'production_mutations': False}
cases = [(False, False), (True, False), (False, True), (True, True)]
for mode in ['legacy', 'v1', 'v2_tmpfs', 'v2_disk']:
    logical = root / (mode + '_logical')
    logical.mkdir()
    channel, remote = socket.socketpair()
    channel.settimeout(120)
    worker = subprocess.Popen([str(binary), '--fullauto-worker-server', str(remote.fileno())],
                              pass_fds=(remote.fileno(),), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True)
    remote.close()
    items = []
    cache = {}
    try:
        for index, (incoherent, controls) in enumerate(cases):
            inputs = {'grids/theta.txt': b'150\n155\n160\n165\n170\n', **cache}
            request = command(str(logical), index, incoherent, controls)
            started = time.monotonic()
            if mode == 'legacy':
                for name, data in inputs.items():
                    path = logical / name
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_bytes(data)
                (logical / f'jobs/j{index}').mkdir(parents=True, exist_ok=True)
                (logical / 'means').mkdir(exist_ok=True)
                send_words(channel, request)
            elif mode == 'v1':
                send_words(channel, ['MBS_WORKER_MEMORY_REQUEST_V1', str(logical)] + request)
                send_files(channel, inputs, 1)
            else:
                estimate = 32 * 1024 * 1024 if mode == 'v2_tmpfs' else 80 * 1024 * 1024
                send_words(channel, ['MBS_WORKER_MEMORY_REQUEST_V2', str(logical), str(estimate)] + request)
                send_files(channel, inputs, 2)
            status = struct.unpack('<i', receive(channel, 4))[0]
            outputs = ({str(path.relative_to(logical)): path.read_bytes() for path in (logical / f'jobs/j{index}').rglob('*') if path.is_file()}
                       if mode == 'legacy' else receive_files(channel, 1 if mode == 'v1' else 2))
            for name, data in outputs.items():
                saved = root / mode / f'case{index}' / name
                saved.parent.mkdir(parents=True, exist_ok=True)
                saved.write_bytes(data)
            log = outputs.get(f'jobs/j{index}/stdout.log', b'').decode(errors='replace')
            assert status == 0, (mode, index, status, log[-4000:])
            if mode == 'v2_disk':
                assert 'Worker memory staging: backend=disk' in log
            if mode.startswith('v2'):
                assert 'Worker memory staging: backend=' in log
            cache.update({name: data for name, data in outputs.items() if name.startswith('means/')})
            items.append({'index': index, 'incoherent': incoherent, 'controls': controls,
                          'seconds': time.monotonic() - started, 'files': list(outputs)})
        send_words(channel, [])
        channel.close()
        assert worker.wait(timeout=15) == 0
        if mode != 'legacy':
            assert not any(logical.iterdir()), ('stage files not cleaned', mode)
        assert not list(Path('/dev/shm').glob(f'mbs-worker-memory-{worker.pid}-*'))
    finally:
        if worker.poll() is None:
            os.killpg(worker.pid, signal.SIGKILL)
            worker.wait(timeout=10)
        channel.close()
    receipt['runs'][mode] = items

for mode in ['v1', 'v2_tmpfs', 'v2_disk']:
    for index, _ in enumerate(cases):
        old = root / 'legacy' / f'case{index}' / f'jobs/j{index}/result'
        new = root / mode / f'case{index}' / f'jobs/j{index}/result'
        for source in old.glob('*.dat'):
            target = new / source.name
            left = [[float(value) for value in row.split()] for row in source.read_text().splitlines()[1:] if row.strip()]
            right = [[float(value) for value in row.split()] for row in target.read_text().splitlines()[1:] if row.strip()]
            assert len(left) == len(right) and all(len(row) == 18 for row in left + right)
            maximum = 0.
            for x, y in zip(left, right):
                assert x[:2] == y[:2] and all(math.isfinite(value) for value in x + y)
                maximum = max(maximum, max(abs(a - b) / max(abs(x[2]), 1e-300) for a, b in zip(x[2:], y[2:])))
            assert maximum < 1e-13, (mode, index, source.name, maximum)
            receipt['matrix_comparisons'].append({'mode': mode, 'case': index, 'file': source.name, 'maximum_scaled_M11': maximum})
        for source in old.glob('*_integrals.tsv'):
            with source.open() as f:
                left = list(csv.DictReader(f, delimiter='\t'))
            with (new / source.name).open() as f:
                right = list(csv.DictReader(f, delimiter='\t'))
            assert len(left) == len(right)
            maximum = 0.
            for x, y in zip(left, right):
                assert x.keys() == y.keys()
                for key in x:
                    try:
                        a0, b0 = float(x[key]), float(y[key])
                    except ValueError:
                        assert x[key] == y[key]
                        continue
                    if math.isnan(a0):
                        assert math.isnan(b0)
                    elif not math.isfinite(a0):
                        assert a0 == b0
                    else:
                        maximum = max(maximum, abs(a0 - b0) / max(1, abs(a0)))
            assert maximum < 1e-12
            receipt['integral_comparisons'].append({'mode': mode, 'case': index, 'file': source.name, 'maximum_scaled_difference': maximum})
receipt['passed'] = True
(root / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'passed': True, 'receipt': str(root / 'receipt.json'),
                  'matrix_comparisons': len(receipt['matrix_comparisons']), 'integral_comparisons': len(receipt['integral_comparisons'])}))
