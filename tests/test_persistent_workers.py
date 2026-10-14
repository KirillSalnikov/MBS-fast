#!/usr/bin/env python3
"""Real-physics reentrancy/cache/paired-depth regression; keeps commands and data."""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import socket
import struct
import subprocess
import numpy as np


def compare(first, second, tolerance):
    a = np.loadtxt(str(first), skiprows=1); b = np.loadtxt(str(second), skiprows=1)
    assert a.shape == b.shape and np.array_equal(a[:, :2], b[:, :2])
    scale = np.maximum(abs(a[:, 2]), 1e-100)[:, None]
    error = float(np.max(abs(a[:, 2:]-b[:, 2:])/scale))
    assert error < tolerance, (str(first), str(second), error)
    return error


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--backend', choices=['cpu', 'cuda'], default='cpu')
    p.add_argument('--trace', choices=['cpu', 'resident', 'hybrid'], default='cpu')
    p.add_argument('--particle-file', type=Path)
    p.add_argument('--case-prefix', default='')
    opt = p.parse_args(); root = opt.output.resolve(); root.mkdir(parents=True)
    binary = str(opt.binary.resolve()); tolerance = 1e-7 if opt.backend == 'cuda' else 1e-12
    base = ['--method', 'po', '--backend', opt.backend, '--particle', '1', '2', '1',
            '--ri', '1.3116', '0', '--wavelength-um', '1', '--max-reflections', '8',
            '--sobol-seed', '64', '23', '--haar-alpha', '--symmetry', '1', '1',
            '--threads', '2', '--orientation-chunk', '32', '--cutoff-profile', 'off',
            '--po-trace-backend', opt.trace, '--no-gpu-trace-prefilter', '--close']
    if opt.backend == 'cpu': base.remove('--no-gpu-trace-prefilter')
    if opt.particle_file:
        index = base.index('--particle'); base[index:index+4] = ['--particle-file', str(opt.particle_file.resolve())]
    def setopt(args, key, values):
        args[args.index(key)+1:args.index(key)+1+len(values)] = values
    cases = [('base', 4, 8, []), ('angles', 8, 16, []),
             ('paired', 4, 8, ['--paired-depth', '4']),
             ('paired_angles', 8, 16, ['--paired-depth', '4']),
             ('cutoff', 4, 8, ['--paired-depth', '4', '--beam-cutoff', '.001']),
             ('cutoff_angles', 8, 16, ['--paired-depth', '4', '--beam-cutoff', '.001']),
             ('absorbing', 4, 8, ['--abs']), ('absorbing_angles', 8, 16, ['--abs']),
             ('incoherent', 4, 8, ['--incoherent']), ('incoherent_angles', 8, 16, ['--incoherent']),
             ('nonconvex', 4, 8, ['--geometry', 'nonconvex']),
             ('nonconvex_angles', 8, 16, ['--geometry', 'nonconvex']),
             ('controls', 4, 8, ['--analytic-facet-average', '--analytic-shadow-control', 'circular']),
             ('controls_angles', 8, 16, ['--analytic-facet-average', '--analytic-shadow-control', 'circular']),
             ('changed_seed', 4, 8, []), ('changed_shape', 4, 8, [])]
    connection, peer = socket.socketpair(); connection.settimeout(60)
    worker = subprocess.Popen([binary, '--fullauto-worker-server', str(peer.fileno())], pass_fds=(peer.fileno(),))
    peer.close(); results = []
    def run(args, destination, persistent):
        destination.mkdir(); log = destination/'stdout.log'
        command = args+['--output', str(destination/'result'), '--trace-cache-mib', '256' if persistent else '0']
        (destination/'command.json').write_text(json.dumps([binary]+command, indent=2)+'\n')
        if persistent:
            words = [str(log), '64', 'thread']+command
            connection.sendall(struct.pack('I', len(words))+b''.join(struct.pack('I', len(w.encode()))+w.encode() for w in words))
            response = b''
            while len(response) < 4:
                part = connection.recv(4-len(response)); assert part, str(log)
                response += part
            assert struct.unpack('i', response)[0] == 0, log.read_text()[-2000:]
        else:
            with log.open('w') as stream:
                subprocess.run([binary]+command, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=60)
        return destination/'result/result.dat'
    try:
        for name, phi, theta, extra in cases:
            if not name.startswith(opt.case_prefix): continue
            if opt.particle_file and name == 'changed_shape': continue
            args = base+['--scattering-grid', '0', '180', str(phi), str(theta)]+extra
            if name.startswith('cutoff'):
                index = args.index('--cutoff-profile'); del args[index:index+2]
                args += ['--trace-cutoff', '0']
            if name.startswith('absorbing'): setopt(args, '--ri', ['1.3116', '.01'])
            if name == 'changed_seed': setopt(args, '--sobol-seed', ['64', '37'])
            if name == 'changed_shape': setopt(args, '--particle', ['1', '3', '.7'])
            reference = list(args); setopt(reference, '--po-trace-backend', ['cpu'])
            fresh = run(reference, root/(name+'_fresh'), False)
            cached = run(args, root/(name+'_cached'), True)
            error = compare(fresh, cached, tolerance)
            log = (cached.parents[1]/'stdout.log').read_text()
            hits = re.search(r'Prepared trace reuse: hits=(\d+)', log)
            cpu_hits = int(hits.group(1)) if hits else 0
            gpu_hits = log.count('GPU raw trace reuse: hit')
            if name.endswith('angles'): assert cpu_hits or gpu_hits, (name, log[-2000:])
            if name in ['changed_seed', 'changed_shape']: assert cpu_hits == gpu_hits == 0
            # Energy accounting must survive reuse, not only matrix accumulation.
            with fresh.with_name('result_integrals.tsv').open() as f: integral_a = list(csv.DictReader(f, delimiter='\t'))
            with cached.with_name('result_integrals.tsv').open() as f: integral_b = list(csv.DictReader(f, delimiter='\t'))
            assert len(integral_a) == len(integral_b)
            for a, b in zip(integral_a, integral_b):
                for key in ['Cext', 'Csca', 'Cabs', 'Qext', 'Qsca', 'Qabs']:
                    av, bv = float(a[key]), float(b[key])
                    assert abs(av-bv) <= tolerance*max(abs(av), abs(float(a['Csca'])), 1e-100), (name, key, av, bv)
            paired_error = None
            if '--paired-depth' in args:
                standalone = list(reference); depth = standalone.pop(standalone.index('--paired-depth')+1); standalone.remove('--paired-depth')
                setopt(standalone, '--max-reflections', [depth])
                lower = run(standalone, root/(name+'_lower'), False)
                paired_error = compare(lower, cached.with_name('result_depth'+depth+'.dat'), tolerance)
            results.append(dict(case=name, all16_error=error, paired_error=paired_error, cpu_hits=cpu_hits, gpu_hits=gpu_hits))
            (root/'comparison.json').write_text(json.dumps(results, indent=2)+'\n')
            print(name, error, paired_error, cpu_hits, gpu_hits, flush=True)
    finally:
        connection.close()
        try: worker.wait(timeout=10)
        except subprocess.TimeoutExpired: worker.kill(); worker.wait()


if __name__ == '__main__': main()
