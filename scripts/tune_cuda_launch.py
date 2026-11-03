#!/usr/bin/env python3
"""Explore CUDA launch settings with a fixed, fully validated optical workload.

This never changes defaults or produces a hardware-ranking score. Chunk size
and environment overrides are recorded separately for each candidate.
"""
import argparse
import copy
import json
from pathlib import Path

import benchmark_hardware as bench


VARIANTS = [
    ('default', {}, 64),
    ('block128-3d', {'MBS_GPU_BLOCK': '128', 'MBS_GPU_THREAD_GRID_3D': '1'}, 64),
    ('block256-3d', {'MBS_GPU_BLOCK': '256', 'MBS_GPU_THREAD_GRID_3D': '1'}, 64),
    ('block64-flat', {'MBS_GPU_THREAD_GRID_3D': '0'}, 64),
    ('block128-flat', {'MBS_GPU_BLOCK': '128'}, 64),
    ('warp64', {'MBS_GPU_WARP_BEAMS': '1'}, 64),
    ('stage-mueller', {'MBS_GPU_STAGE_MUELLER': '1'}, 64),
    ('no-split4', {'MBS_GPU_COMPACT_BEAM4_SPLIT': '0'}, 64),
    ('chunk128', {}, 128),
    ('chunk256', {}, 256),
    ('default-end', {}, 64),
]


def main():
    suite = json.loads((bench.SUITE / 'suite.json').read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prefix', type=Path)
    parser.add_argument('--profile', choices=[p for p in suite['profiles'] if p.startswith('gpu')], default='gpu-fp32')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--device', type=int, default=0)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--timeout', type=int, default=180)
    parser.add_argument('--variants', help='Comma-separated subset of variant names')
    args = parser.parse_args()
    if args.repeats < 3 or args.threads < 1 or args.device < 0 or args.timeout < 1:
        parser.error('At least three repeats, positive threads/timeout and nonnegative device required')
    args.output = args.output.expanduser().absolute()
    if args.output.exists():
        parser.error('Use a fresh output directory')
    selected = set(args.variants.split(',')) if args.variants else {v[0] for v in VARIANTS}
    if selected - {v[0] for v in VARIANTS}:
        parser.error('Unknown variant name')
    args.binary = args.python = None
    args.label = 'CUDA launch tuning'
    build = bench.build_info(args, suite)
    hardware = bench.hardware(args, build)
    if hardware['before'] and (hardware['before'].get('utilization.gpu') or 0) > 10:
        raise RuntimeError('GPU busy; retry when idle')
    case = suite['cases'][0]
    reference = json.loads((bench.SUITE / 'reference' / (case['id'] + '__' + args.profile + '.json')).read_text())
    signature = bench.digest(dict(suite=suite['id'], case=case, profile=args.profile, threads=args.threads))
    if signature != reference['workload_sha256'] or bench.digest(reference['matrix']) != reference['matrix_sha256']:
        raise ValueError('Reference/workload mismatch; keep the frozen eight-thread workload')
    limits = suite['profiles'][args.profile]['validation']
    env = bench.environment(args, build)
    report = dict(schema=1, purpose='Exploratory tuning; not a hardware-ranking score',
                  build=build, hardware=hardware, runner_sha256=bench.sha(Path(__file__)),
                  case=case, profile=args.profile, threads=args.threads,
                  reference_matrix_sha256=reference['matrix_sha256'], variants=[])
    args.output.mkdir(parents=True)
    for name, overrides, chunk in VARIANTS:
        if name not in selected:
            continue
        candidate = copy.deepcopy(case)
        candidate['arguments'][candidate['arguments'].index('--orientation-chunk') + 1] = str(chunk)
        record = dict(name=name, environment=overrides, chunk=chunk,
                      arguments=candidate['arguments'], runs=[], passed=False)
        report['variants'].append(record)
        try:
            for index in range(args.repeats + 1):
                directory = args.output / name / str(index)
                command = bench.command_for(candidate, build, args.profile, args.threads, directory)
                run = bench.timed_process(command, dict(env, **overrides), directory, args.device, args.timeout)
                matrix, metrics, work = bench.read_result(candidate, directory, args.profile)
                error = bench.matrix_error(matrix, reference['matrix'])
                if any(error[key] > value for key, value in limits.items()):
                    raise ValueError('Mueller reference gate failed: ' + str(error))
                if bench.sha(build['native_binary']) != build['binary_sha256']:
                    raise RuntimeError('Native binary changed during tuning')
                run.update(warmup=index == 0, command=command, metrics=metrics,
                           work=work, validation=error, matrix_sha256=bench.digest(matrix))
                record['runs'].append(run)
                bench.write(args.output / 'results.json', report)
            measured = record['runs'][1:]
            record['statistics'] = {key: bench.statistics_of([r[key] for r in measured])
                                    for key in ['wall_s', 'cpu_user_s', 'cpu_system_s', 'peak_child_rss_MiB']}
            record['kernel_statistics'] = bench.statistics_of([r['metrics']['gpu_kernels_s'] for r in measured])
            record['timing_stable'] = record['statistics']['wall_s']['relative_MAD'] <= .1
            record['passed'] = True
            print(name, record['statistics']['wall_s'], flush=True)
        except (RuntimeError, ValueError, OSError) as error:
            record['error'] = str(error)
            print(name, 'FAILED', str(error), flush=True)
        bench.write(args.output / 'results.json', report)
    if any(not v['passed'] for v in report['variants']):
        raise SystemExit('One or more candidates failed; inspect results.json. Failed candidates must not be ranked.')


if __name__ == '__main__':
    main()
