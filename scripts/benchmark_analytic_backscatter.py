#!/usr/bin/env python3
"""Compare raw and analytic-control MBS means on independent Owen scrambles.

Example:
  python3 scripts/benchmark_analytic_backscatter.py --orientations 512 2048 \
    --output results/analytic_check -- --particle 1 101.02 101.02 \
    --refractive-index 1.31 0 --wavelength-um .532 --max-reflections 20 \
    --scattering-grid 179 180 4 4 --cutoff-profile off

The native calculation supplies both estimators from the same coherent fields.
Reported scatter is an empirical error estimate, not an accuracy certificate
or a measured runtime-to-tolerance speedup. Uses only the Python standard library.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def positive(text):
    value = int(text)
    if value <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return value


def nonnegative(text):
    value = int(text)
    if value < 0 or value > 2147483647:
        raise argparse.ArgumentTypeError('must be in [0, 2147483647]')
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--binary', type=Path, default=ROOT/'cpu/bin/mbs_po_mpi')
    parser.add_argument('--orientations', nargs='+', type=positive, default=[512,2048])
    parser.add_argument('--seeds', nargs='+', type=nonnegative, default=[11,23,37,53,71,89,107,131])
    parser.add_argument('--threads', type=positive, default=4)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('solver_arguments', nargs=argparse.REMAINDER,
                        help='particle, physics and fixed-grid options after --')
    args = parser.parse_args()
    physics = args.solver_arguments
    if physics and physics[0] == '--':
        physics = physics[1:]
    if not physics:
        parser.error('provide particle, refractive index, wavelength and reflection/grid options after --')
    forbidden = {'--method','--po','--go','--backend','--gpu','--cpu',
                 '--sobol','--sobol-seed','--sobol_seed','--hammersley','--lattice',
                 '--fixed','--fixed-orientation','--orientation-file','--orientfile',
                 '--analytic-backscatter','--output','--ofile','--threads','--close'}
    if any(option in forbidden for option in physics):
        parser.error('method/backend, orientation rule, analytic flag, output and threads are set by this script')
    if len(args.seeds) < 2 or len(set(args.seeds)) != len(args.seeds):
        parser.error('provide at least two distinct independent seeds')
    if len(set(args.orientations)) != len(args.orientations):
        parser.error('orientation counts must be distinct')
    binary = args.binary.resolve()
    if not binary.is_file():
        parser.error(f'binary does not exist: {binary}')
    output = args.output.resolve()
    conflicts = [output/f'N{count}_seed{seed}' for count in args.orientations for seed in args.seeds
                 if (output/f'N{count}_seed{seed}').exists()]
    if conflicts:
        parser.error(f'output already contains a run: {conflicts[0]}; choose a new --output directory')
    output.mkdir(parents=True,exist_ok=True)
    env = {key:value for key,value in os.environ.items() if not key.startswith('MBS_')}
    manifest = dict(binary=str(binary), binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    solver_arguments=physics, orientations=args.orientations, seeds=args.seeds,
                    threads=args.threads, scope='M11 at exactly 180 degrees; same full coherent fields for both estimators')
    (output/'configuration.json').write_text(json.dumps(manifest,indent=2)+'\n')
    runs, summaries = [], []
    for count in args.orientations:
        for seed in args.seeds:
            name=f'N{count}_seed{seed}'
            target=output/name
            command=[str(binary),'--method','po','--backend','cpu',*physics,
                     '--sobol-seed',str(count),str(seed),'--analytic-backscatter',
                     '--threads',str(args.threads),'--output',str(target),'--close']
            started=time.monotonic()
            with (output/(name+'.log')).open('w') as log:
                completed=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            if completed.returncode:
                raise SystemExit(f'MBS failed for {name}; see {output/(name+".log")}')
            with (target/(name+'_analytic_backscatter.tsv')).open() as stream:
                row=next(csv.DictReader(stream,delimiter='\t'))
            raw=float(row['raw_M11']); hybrid=float(row['hybrid_M11'])
            known=float(row['analytic_control_mean']); residual=float(row['residual_mean'])
            if not all(math.isfinite(value) for value in [raw,hybrid,known,residual]):
                raise SystemExit(f'nonfinite result for {name}')
            if abs(hybrid-known-residual)>1e-9*max(abs(raw),abs(known),1e-100):
                raise SystemExit(f'full residual identity failed for {name}')
            runs.append(dict(N=count,seed=seed,raw=raw,hybrid=hybrid,
                             seconds=time.monotonic()-started,
                             setup_profile_seconds=float(row['setup_profile_seconds']),
                             known_control=known, moment_refinement=float(row['moment_refinement']),
                             command=command))
            (output/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
            print(f'{name}: raw={raw:.8g}, hybrid={hybrid:.8g}',flush=True)
        rows=[run for run in runs if run['N']==count]
        raw=[run['raw'] for run in rows]; hybrid=[run['hybrid'] for run in rows]
        raw_mean,hybrid_mean=statistics.mean(raw),statistics.mean(hybrid)
        raw_se=statistics.stdev(raw)/math.sqrt(len(rows))
        hybrid_se=statistics.stdev(hybrid)/math.sqrt(len(rows))
        raw_relative=100*raw_se/abs(raw_mean) if raw_mean else None
        hybrid_relative=100*hybrid_se/abs(hybrid_mean) if hybrid_mean else None
        summary=dict(N=count,scrambles=len(rows),raw_mean=raw_mean,hybrid_mean=hybrid_mean,
                     raw_standard_error=raw_se,hybrid_standard_error=hybrid_se,
                     raw_relative_SE_percent=raw_relative,hybrid_relative_SE_percent=hybrid_relative,
                     variance_ratio=raw_se**2/hybrid_se**2 if hybrid_se else None,
                     mean_wall_seconds=statistics.mean(run['seconds'] for run in rows),
                     mean_setup_profile_seconds=statistics.mean(run['setup_profile_seconds'] for run in rows))
        summaries.append(summary)
        (output/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
        with (output/'summary.tsv').open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(summary),delimiter='\t')
            writer.writeheader();writer.writerows(summaries)
        print(f'N={count}: raw SE={raw_relative}%, hybrid SE={hybrid_relative}%; '
              f'variance ratio={summary["variance_ratio"]}',flush=True)
    print('Variance ratios are not measured solver speedups or certified accuracy.',flush=True)


if __name__ == '__main__':
    main()
