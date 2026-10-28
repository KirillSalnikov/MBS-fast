#!/usr/bin/env python3
"""Differential native benchmark: same poses, all Mueller entries and integrals.

Run old/new binaries in alternating order. Preserve every command, DAT and
phase log; report tracing separately from complete processing. No speed gate
is imposed on shared servers. Matrix agreement is mandatory.
"""
import argparse
import json
import os
import statistics
import subprocess
import time
from pathlib import Path
from test_droxtal_resident import load, difference, integral_difference


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--before', type=Path, required=True)
    p.add_argument('--after', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--backend', choices=['cpu', 'cuda'], default='cpu')
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--tolerance', type=float, default=1e-12)
    a = p.parse_args()
    root = a.output.resolve()
    root.mkdir(parents=True)
    repo = Path(__file__).resolve().parents[1]
    grid = root / 'angles.txt'
    grid.write_text('0\n0.1\n79.9976\n150\n163.75909223159732\n170\n179.9\n180\n')
    small_grid = root / 'timing_angles.txt'
    small_grid.write_text('79.9976\n163.75909223159732\n')
    common = ['--backend', a.backend, '--wavelength-um', '.532', '--ri', '1.3116', '0',
              '--threads', str(a.threads), '--haar-alpha', '--theta-grid-file', str(grid),
              '--phi-points', '8', '--sobol-offset', '256', '--orientation-chunk', '64',
              '--trace-cache-mib', '0', '--profile-phases', '--close']
    if a.backend == 'cuda':
        common += ['--po-trace-backend', 'cpu', '--orientation-pipeline']
    env = os.environ.copy()
    env['OMP_WAIT_POLICY'] = 'PASSIVE'
    env['MBS_HOST_MEM_BUDGET_MB'] = '2048'
    binaries = {'before': str(a.before.resolve()), 'after': str(a.after.resolve())}
    droxtal = ['--particle', '4', '44.565', '--symmetry', '2', '6']
    column = ['--particle', '1', '398.1', '138.9', '--symmetry', '2', '6']
    asym = ['--particle-file', str(repo/'tests/data/resident_asymmetric.particle'),
            '--resize-dmax-um', '25', '--symmetry', '1', '1']
    concave = ['--particle-file', str(repo/'examples/particles/concave_hexagonal.particle'),
               '--resize-dmax-um', '5']
    records = []

    def run(label, version, args):
        dest = root/(label+'_'+version)
        dest.mkdir()
        cmd = [binaries[version]]+args+['--output', str(dest/'result')]
        (dest/'command.args').write_text('\n'.join(cmd)+'\n')
        start = time.monotonic()
        with (dest/'stdout.log').open('w') as stream:
            result = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT,
                                    env=env, timeout=300)
        assert result.returncode == 0, (label, version, dest/'stdout.log')
        phase = {}
        for line in (dest/'result/result_phase_timings.tsv').read_text().splitlines()[1:]:
            k, v, _ = line.split()
            phase[k] = float(v)
        return dict(directory=str(dest), seconds=time.monotonic()-start, phase=phase)

    def compare(label, args, order):
        results = {version: run(label, version, args) for version in order}
        before, after = [Path(results[k]['directory'])/'result' for k in ['before', 'after']]
        errors = {}
        products = sorted(before.glob('result*.dat'))
        assert products
        assert [f.name for f in products] == sorted(f.name for f in after.glob('result*.dat'))
        for f in products:
            errors[f.name] = difference(load(f), load(after/f.name))
        for f in sorted(before.glob('result*_integrals.tsv')):
            errors[f.name] = integral_difference(f, after/f.name)
        passed = max(errors.values()) <= a.tolerance
        record = dict(case=label, passed=passed, errors=errors, runs=results)
        records.append(record)
        (root/'comparison.json').write_text(json.dumps(records, indent=2)+'\n')
        print(json.dumps(dict(case=label, passed=passed, max_error=max(errors.values()),
                              phases={k:v['phase'] for k,v in results.items()})), flush=True)
        assert passed, record

    for name, shape in [('droxtal89', droxtal), ('column398', column),
                        ('asymmetric', asym), ('concave', concave)]:
        for incoherent in [False, True]:
            args = common+['--method', 'po']+shape+['--sobol-seed', '64', '23',
                '--max-reflections', '12' if name != 'concave' else '4',
                '--cutoff-profile', 'balanced' if incoherent else 'off']
            if name != 'concave':
                args += ['--paired-depth', '10', '--orientation-prefix-output', '32']
            if incoherent:
                args += ['--incoherent']
            compare(name+('_incoherent' if incoherent else '_coherent'), args,
                    ['after', 'before'] if incoherent else ['before', 'after'])
    for cutoff in ['off', 'balanced']:
        args = common+['--method', 'po']+droxtal+['--sobol-seed', '64', '7013',
            '--max-reflections', '10', '--cutoff-profile', cutoff, '--absorption']
        args[args.index('--ri')+2] = '.00042407'
        compare('absorbing_'+cutoff, args, ['before', 'after'])

    for name, shape, count in [('droxtal89', droxtal, 256), ('column398', column, 1024)]:
        for repeat in range(a.repeats):
            args = common+['--method', 'po']+shape+['--sobol-seed', str(count), '23',
                '--max-reflections', '12', '--cutoff-profile', 'off', '--incoherent']
            args[args.index('--theta-grid-file')+1] = str(small_grid)
            compare('timing_'+name+'_'+str(repeat), args,
                    ['before', 'after'] if repeat%2 == 0 else ['after', 'before'])
    summary = {}
    for name in ['droxtal89', 'column398']:
        group = [r for r in records if r['case'].startswith('timing_'+name+'_')]
        if not group:
            continue
        summary[name] = {}
        for phase in ['trace', 'trace_prepare', 'processing']:
            medians = {v: statistics.median(r['runs'][v]['phase'][phase] for r in group)
                       for v in binaries}
            summary[name][phase] = dict(medians, speedup=medians['before']/medians['after'])
    (root/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
