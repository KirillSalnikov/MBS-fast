#!/usr/bin/env python3
"""Independent all-Mueller checks of adaptive scheduling on five real shapes.

Explicit phi/pilot isolates allocation; the staging probe separately exercises
late phi screening, manual overrides and interrupted promotion recovery.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from audit_fullauto_v2 import run as audit
from test_fullauto_arbitrary import asymmetric_tetrahedron, rotated_box

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    tetra, box = out/'tetra.particle', out/'box.particle'
    asymmetric_tetrahedron(tetra)
    rotated_box(box)
    cases = [('prism', ['--particle', '1', '.2', '.1'], '0', []),
             ('asymmetric', ['--particle-file', str(tetra)], '0', []),
             ('absorbing', ['--particle-file', str(box)], '.02', []),
             ('concave', ['--particle-file', str(ROOT/'examples/particles/concave_hexagonal.particle'),
                          '--resize-dmax-um', '.3'], '0', []),
             ('incoherent', ['--particle', '1', '.2', '.1'], '0', ['--incoherent'])]
    binaries = [('new', args.binary.resolve())]
    if args.reference:
        binaries.insert(0, ('reference', args.reference.resolve()))
    results = []
    for name, shape, imag, extra in cases:
        pair = []
        for label, binary in binaries:
            dest = out/(name+'_'+label)
            command = [str(binary), '--method', 'po', '--backend', 'cpu']+shape+extra+[
                '--ri', '1.3116', imag, '--wavelength-um', '1', '--max-reflections', '4',
                '--threads', '1', '--fullauto', '.02', '--fullauto-resources', 'off',
                '--fullauto-mirror', 'off', '--fullauto-cutoffs', 'off', '--cutoff-profile', 'off',
                '--symmetry', '1', '1', '--phi-points', '4', '--scattering-grid', '0', '25', '4', '16',
                '--fullauto-workers', 'persistent', '--fullauto-cpu-workers', '2', '--fullauto-pilot', '32',
                '--fullauto-initial', '256', '--fullauto-min-correction', '64', '--fullauto-theta-start', '5',
                '--fullauto-max-rounds', '32', '--fullauto-orientation-block', '256',
                '--max-orientations', '65536', '--output', str(dest), '--close']
            (out/(name+'_'+label+'.args')).write_text('\n'.join(command)+'\n')
            started = time.monotonic()
            with (out/(name+'_'+label+'.log')).open('w') as stream:
                subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=600)
            report = audit(dest, ROOT/'scripts/gpu_campaign.py', False)
            report.update(case=name, label=label, seconds=time.monotonic()-started)
            markers = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in (dest/'jobs').glob('*/complete.txt')}
            with (out/(name+'_'+label+'_resume.log')).open('w') as stream:
                subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True, timeout=60)
            assert all(hashlib.sha256(p.read_bytes()).hexdigest() == value for p, value in markers.items())
            report['completed_blocks_preserved'] = len(markers)
            pair.append(report)
            results.append(report)
            print(name, label, 'PASS', flush=True)
        if len(pair) == 2:
            matrices = [np.loadtxt(out/(name+'_'+label)/'mueller_fullauto.dat', skiprows=1) for label, _ in binaries]
            difference = float(np.max(abs(matrices[0][:, 2:]-matrices[1][:, 2:])/matrices[0][:, 2, None]))
            bound = sum(r['statistical95_scaled_M11']+r['guard_residual95_scaled_M11'] for r in pair)
            assert difference <= bound, (name, difference, bound)
            for report in pair:
                report['all16_difference_scaled_M11'] = difference
        (out/'results.json').write_text(json.dumps(results, indent=2)+'\n')


if __name__ == '__main__':
    main()
