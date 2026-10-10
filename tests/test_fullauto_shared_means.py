#!/usr/bin/env python3
"""Verify one authoritative analytic mean cache beyond 4097 theta rows.

Uses completed worker outputs and an independent direct-grid audit. Optional
previous binary checks the old per-grid path on the same physical problem.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from audit_fullauto import run as audit

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT/'bin/mbs_po')
    parser.add_argument('--reference-binary', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    options = parser.parse_args()
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    records = []
    variants = [('shared', options.binary)]
    if options.reference_binary:
        variants.append(('per_grid', options.reference_binary))
    for name, binary in variants:
        result = output/name
        command = [str(binary.resolve()), '--method', 'po', '--backend', 'cpu',
                   '--particle', '1', '.2', '.1', '--ri', '1.31', '0',
                   '--wavelength-um', '1', '--max-reflections', '1',
                   '--fullauto', '.01', '--fullauto-mirror', 'off',
                   '--scattering-grid', '0', '25', '4', '4098', '--phi-points', '4',
                   '--fullauto-pilot', '32', '--fullauto-initial', '256',
                   '--fullauto-theta-start', '5', '--fullauto-max-rounds', '16',
                   '--fullauto-min-correction', '32', '--max-orientations', '65536',
                   '--threads', '2', '--output', str(result)]
        (output/(name+'_command.args')).write_text('\n'.join(command)+'\n')
        started = time.monotonic()
        with (output/(name+'.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=1800)
        seconds = time.monotonic()-started
        state = json.loads((result/'fullauto_status.json').read_text())
        assert state['requested_theta'] == 4099
        assert state['shared_analytic_reference'] == (name == 'shared')
        commands = [(job/'command.args').read_text().splitlines()
                    for job in (result/'jobs').iterdir() if (job/'complete.txt').exists()]
        if name == 'shared':
            assert all('--analytic-mean-reference-grid' in c for c in commands)
            assert len(list(result.glob('*.cache'))) == 1
        report = audit(result, ROOT/'scripts/gpu_campaign.py', dense=True, dense_workers=2)
        records.append(dict(variant=name, controller_seconds=seconds, state=state, audit=report))
    comparison = dict(cases=records)
    if len(records) == 2:
        a, b = [np.loadtxt(output/name/'mueller_fullauto.dat', skiprows=1) for name in ('shared','per_grid')]
        difference = float(np.max(abs(a[:,2:]-b[:,2:])/b[:,2,None]))
        comparison['max_matrix_difference_scaled_M11'] = difference
        # A mean cache may change the adaptation history through timings;
        # use the measured independent pointwise intervals in that case.
        same_counts = records[0]['state']['counts_per_seed'] == records[1]['state']['counts_per_seed']
        comparison['same_counts'] = same_counts
        if same_counts:
            assert difference < 1e-10, comparison
        else:
            assert difference <= sum(r['state']['max_pointwise95_scaled_M11'] for r in records), comparison
    (output/'comparison.json').write_text(json.dumps(comparison, indent=2)+'\n')
    print(json.dumps(comparison, indent=2))


if __name__ == '__main__':
    main()
