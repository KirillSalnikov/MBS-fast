#!/usr/bin/env python3
"""Targeted scheduling and fresh holdouts after a real failed native audit.

Build fullauto_scheduler_probe and fullauto_queue_probe with the CPU objects
(except main.o and FullAuto.o). --worker is the unmodified native executable.
All commands, worker matrices, checkpoints and logs remain under --output.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scheduler-probe', type=Path, required=True)
    parser.add_argument('--queue-probe', type=Path, required=True)
    parser.add_argument('--worker', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    command = [str(args.scheduler_probe.resolve()), str(out/'decisions'),
               str(ROOT/'tests/data/fullauto_priority_backscatter.txt')]
    with (out/'decisions.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    command = [str(args.queue_probe.resolve()), '--method', 'po', '--backend', 'cpu',
               '--particle', '1', '.2', '.1', '--ri', '1.3116', '0', '--wavelength-um', '1',
               '--max-reflections', '4', '--threads', '1', '--fullauto', '.02',
               '--fullauto-analytic-controls', 'off', '--fullauto-mirror', 'off',
               '--fullauto-resources', 'off', '--cutoff-profile', 'off', '--phi-points', '4',
               '--scattering-grid', '0', '180', '4', '8', '--fullauto-initial', '256',
               '--fullauto-orientation-block', '128', '--fullauto-depth-sharing', 'off',
               '--fullauto-work-queue', 'on', '--output', str(out/'failed_audit')]
    (out/'failed_audit.args').write_text('\n'.join(command)+'\n')
    started = time.monotonic()
    with (out/'failed_audit.log').open('w') as log:
        subprocess.run(command, env=dict(os.environ, TASK_MBS_WORKER=str(args.worker.resolve()),
                       TASK_MBS_AUDIT_FAILURE='1'), stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=600)
    audits = [json.loads(p.read_text()) for p in
              (out/'failed_audit').glob('independent_final_validation_serial*.json')]
    seeds = [s for audit in audits for s in audit['seeds']]
    assert len(audits) >= 2 and all(not audit['passed'] for audit in audits)
    assert len(seeds) == len(set(seeds)), 'Exposed holdout reused after recovery'
    report = dict(passed=True, seconds=time.monotonic()-started,
                  independent_failed_audits=len(audits), seeds=seeds,
                  recovery_preserved_fresh_holdout=True)
    (out/'result.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
