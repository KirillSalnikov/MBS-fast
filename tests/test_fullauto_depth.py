#!/usr/bin/env python3
"""Check geometric initial budgets, fixed override, and automatic physical-depth audits."""
import argparse
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT/'bin/mbs_po')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve();output.mkdir(parents=True, exist_ok=True)
    records = []
    # Each one-round run computes real matrices and deliberately stops before
    # stability certification. This measures the controller's actual start.
    for name, height, fixed in [('small', '.2', None), ('larger', '1', None),
                                ('explicit', '1', '256')]:
        result = output/name
        command = [str(args.binary.resolve()), '--method', 'po', '--backend', 'cpu',
            '--particle', '1', height, str(float(height)/2), '--ri', '1.31', '0',
            '--wavelength-um', '1', '--max-reflections', '1', '--fullauto', '.1',
            '--scattering-grid', '0', '25', '1', '2', '--fullauto-theta-start', '3',
            '--fullauto-max-rounds', '1', '--max-orientations', '262144',
            '--fullauto-analytic-controls', 'off', '--fullauto-resources', 'off',
            '--fullauto-cutoffs', 'off', '--fullauto-mirror', 'off',
            '--fullauto-rotation-symmetry', 'off', '--fullauto-final-audit', 'off',
            '--fullauto-depth-sharing', 'off', '--threads', '2', '--output', str(result), '--close']
        if fixed:command += ['--fullauto-initial', fixed]
        (output/(name+'.args')).write_text('\n'.join(command)+'\n')
        with (output/(name+'.log')).open('w') as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=300)
        assert completed.returncode == 3, (completed.returncode, (output/(name+'.log')).read_text()[-2000:])
        parameters = json.loads((result/'fullauto_parameters.json').read_text())
        edge = parameters['maximum_edge_um'];step = .69/edge/16
        beta, gamma = max(3, math.ceil(math.pi/step)), max(3, math.ceil(2*math.pi/step))
        desired = max(512, beta*gamma)
        expected = int(fixed) if fixed else 2**math.ceil(math.log2(min(desired, 262144//4)))
        assert parameters['initial_count'] == expected
        assert parameters['oldauto_N_beta_estimate'] == beta
        assert parameters['oldauto_N_gamma_estimate'] == gamma
        assert parameters['physical_target_depth'] == 1 and not parameters['automatic_depth']
        assert (result/'intermediate_current.json').exists()
        records.append(dict(case=name, maximum_edge_um=edge, beta=beta, gamma=gamma,
                            estimated=beta*gamma, selected=expected, explicit_initial=fixed))
    (output/'initial_budget_report.json').write_text(json.dumps(records, indent=2)+'\n')
    print(json.dumps(records, indent=2))


if __name__ == '__main__':main()
