#!/usr/bin/env python3
"""Validate worker resource rejection before any CUDA native job is dispatched."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--output', type=Path)
args = parser.parse_args()
binary = str(args.binary.resolve())
base = [binary, '--method', 'po', '--backend', 'cuda', '--particle', '1', '.2', '.1',
        '--ri', '1.31', '0', '--wavelength-um', '1', '--fullauto', '.04',
        '--fullauto-resources', 'off', '--fullauto-analytic-controls', 'off']
cases = [
    ('count', ['--fullauto-gpu-workers', '3'], {}, 'must be 1, 2 or 4'),
    ('threads_missing', ['--fullauto-gpu-workers', '4'], {}, 'require explicit dispatch threads'),
    ('automatic_devices', ['--fullauto-gpu-workers', '4', '--fullauto-dispatch-threads', '1',
                          '--fullauto-gpus', 'auto'], {}, 'require fixed GPU selection'),
    ('cpu_oversubscription', ['--fullauto-gpu-workers', '4', '--fullauto-dispatch-threads', '1000000'],
                           {}, 'exceed allowed physical CPU cores'),
    ('gpu_budget', ['--fullauto-gpu-workers', '4', '--fullauto-dispatch-threads', '1'],
                   {'MBS_GPU_MEM_FRACTION': '0.5'}, 'exceeds shared GPU budget'),
    ('nonfinite_gpu_budget', ['--fullauto-gpu-workers', '4', '--fullauto-dispatch-threads', '1'],
                   {'MBS_GPU_MEM_FRACTION': 'nan'}, 'exceeds shared GPU budget'),
    ('invalid_device', ['--fullauto-gpu-workers', '4', '--fullauto-dispatch-threads', '1',
                       '--fullauto-gpus', '999999'], {}, 'invalid/duplicate visible GPU index'),
    ('compatible_override', ['--fullauto-gpu-workers', '1', '--fullauto-resume-compatible', binary],
                           {}, 'conflict with compatible dispatch overrides'),
]
records = []
with tempfile.TemporaryDirectory(prefix='mbs_gpu_workers_') as temporary:
    for name, flags, changes, expected in cases:
        output = Path(temporary) / name
        env = dict(os.environ)
        env.pop('MBS_GPU_MEM_FRACTION', None)
        env.pop('MBS_HOST_MEM_BUDGET_MB', None)
        env.update(changes)
        run = subprocess.run(base + flags + ['--output', str(output)], env=env,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True, timeout=30)
        assert run.returncode != 0 and expected in run.stdout, (name, run.returncode, run.stdout)
        history = output / 'fullauto_history.jsonl'
        assert not history.exists() or not any(
            json.loads(line).get("event") == "worker_start"
            for line in history.read_text().splitlines() if line.strip()), name
        records.append({'case': name, 'returncode': run.returncode, 'rejected_before_dispatch': True})
if args.output:
    args.output.write_text(json.dumps(records, indent=2))
print(json.dumps(records, indent=2))
