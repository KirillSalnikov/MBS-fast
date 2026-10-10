#!/usr/bin/env python3
"""End-to-end accuracy gate for native fullauto on file-based particles.

Runs automatic theta selection and an independent direct-grid NumPy/SciPy
audit of all 16 Mueller entries. Keeps commands, matrices and logs for review.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

from audit_fullauto import run as audit

ROOT = Path(__file__).resolve().parents[1]


def rotated_box(path):
    # Unequal sides and a rotation remove the body xz reflection assumed by
    # the prism shortcut. Face order follows outward normals.
    vertices = [(-.5, -.4, -.3), (.5, -.4, -.3),
                (.5, .4, -.3), (-.5, .4, -.3),
                (-.5, -.4, .3), (.5, -.4, .3),
                (.5, .4, .3), (-.5, .4, .3)]
    faces = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    c, s = math.cos(.613), math.sin(.613)
    rotated = [(c*x-s*y, s*x+c*y, z) for x, y, z in vertices]
    text = '0\n0\n180 360\n\n'
    for face in faces:
        text += ''.join(' '.join(format(v, '.17g') for v in rotated[i])+'\n'
                        for i in face)+'\n'
    path.write_text(text)


def asymmetric_tetrahedron(path):
    # Unequal edges and an offset apex give no body xz reflection. Orient
    # each triangular face outwards rather than relying on vertex order.
    vertices = [(-.47, -.31, -.23), (.61, -.27, -.19),
                (-.16, .53, -.14), (.12, .07, .68)]
    center = [sum(p[k] for p in vertices)/4 for k in range(3)]
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
    text = '0\n0\n180 360\n\n'
    for face in faces:
        a, b, c = [vertices[i] for i in face]
        u, v = [[p[k]-a[k] for k in range(3)] for p in (b, c)]
        normal = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2],
                  u[0]*v[1]-u[1]*v[0]]
        if sum(normal[k]*(center[k]-a[k]) for k in range(3)) > 0:
            face = tuple(reversed(face))
        text += ''.join(' '.join(format(v, '.17g') for v in vertices[i])+'\n'
                        for i in face)+'\n'
    path.write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT/'bin/mbs_po')
    parser.add_argument('--reference-binary', type=Path,
                        help='Also audit and compare the previous implementation at the selected phi.')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--case', choices=['all', 'convex', 'concave', 'asymmetric',
                                         'absorbing'], default='all')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--phi-points', type=int, default=4,
                        help='Fixed phi count; use 0 to test automatic selection.')
    settings = parser.parse_args()
    assert settings.threads > 0
    assert settings.phi_points >= 0
    output = settings.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    shape = output/'rotated_box.particle'
    rotated_box(shape)
    tetrahedron = output/'asymmetric_tetrahedron.particle'
    asymmetric_tetrahedron(tetrahedron)
    cases = [('convex', shape, 2, '1.3116', '0'),
             ('concave', ROOT/'examples/particles/concave_hexagonal.particle', 4,
              '1.3116', '0'),
             ('asymmetric', tetrahedron, 3, '1.3116', '0'),
             ('absorbing', shape, 3, '1.45', '.02')]
    records = []
    for name, particle, depth, real, imaginary in cases:
        if settings.case != 'all' and settings.case != name:
            continue
        result = output/name
        command = [str(settings.binary.resolve()), '--method', 'po', '--backend', 'cpu',
                   '--particle-file', str(particle), '--symmetry', '1', '1',
                   '--ri', real, imaginary, '--wavelength-um', '1',
                   '--max-reflections', str(depth), '--fullauto', '.01',
                   '--fullauto-analytic-controls', 'auto',
                   '--fullauto-pilot', '64',
                   '--fullauto-initial', '256', '--fullauto-min-correction', '64',
                   '--fullauto-theta-start', '9', '--fullauto-max-rounds', '24',
                   '--max-orientations', '131072', '--threads', str(settings.threads),
                   '--output', str(result)]
        if settings.phi_points:
            command += ['--phi-points', str(settings.phi_points)]
        (output/(name+'_command.args')).write_text('\n'.join(command)+'\n')
        started = time.monotonic()
        with (output/(name+'.log')).open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=1800)
        seconds = time.monotonic()-started
        state = json.loads((result/'fullauto_status.json').read_text())
        assert state['interpolation_budget'] <= .0025
        assert state['pass_streak'] >= 2
        assert (result/'dense_theta_pilot.json').exists()
        if name != 'concave':
            assert not state['mirror_gamma'], 'Unverified mirror reduction'
        if not settings.phi_points:
            calibration = json.loads((result/'calibration.json').read_text())
            assert len(calibration['phi_candidates']) > 1
        report = audit(result, ROOT/'scripts/gpu_campaign.py', dense=True)
        record = dict(case=name, controller_seconds=seconds, state=state, audit=report,
                      automatic_phi=settings.phi_points == 0,
                      refractive_index=[float(real), float(imaginary)])
        if settings.reference_binary:
            reference_result = output/(name+'_reference')
            reference_command = command.copy()
            reference_command[0] = str(settings.reference_binary.resolve())
            reference_command[reference_command.index('--output')+1] = str(reference_result)
            if '--phi-points' in reference_command:
                reference_command[reference_command.index('--phi-points')+1] = str(state['phi'])
            else:
                reference_command += ['--phi-points', str(state['phi'])]
            (output/(name+'_reference_command.args')).write_text('\n'.join(reference_command)+'\n')
            started = time.monotonic()
            with (output/(name+'_reference.log')).open('w') as log:
                subprocess.run(reference_command, stdout=log, stderr=subprocess.STDOUT,
                               check=True, timeout=1800)
            reference_seconds = time.monotonic()-started
            reference_report = audit(reference_result, ROOT/'scripts/gpu_campaign.py', dense=True)
            previous = np.loadtxt(reference_result/'mueller_fullauto.dat', skiprows=1)
            current = np.loadtxt(result/'mueller_fullauto.dat', skiprows=1)
            assert np.allclose(current[:, 0], previous[:, 0], rtol=0, atol=1e-10)
            difference = float(np.max(abs(current[:, 2:]-previous[:, 2:])/previous[:, 2, None]))
            bound = (report['statistical95_scaled_M11']+reference_report['statistical95_scaled_M11']+
                     report['independent_dense']['residual_plus95_scaled_M11']+
                     reference_report['independent_dense']['residual_plus95_scaled_M11'])
            assert difference <= bound+1e-12, (difference, bound)
            record['previous_implementation'] = dict(controller_seconds=reference_seconds,
                audit=reference_report, selected_phi=state['phi'],
                max_matrix_difference_scaled_M11=difference, comparison_bound_scaled_M11=bound,
                same_counts=state['counts_per_seed'] == reference_report['final_counts_per_seed'])
        records.append(record)
        (output/'comparison.json').write_text(json.dumps(records, indent=2)+'\n')


if __name__ == '__main__':
    main()
