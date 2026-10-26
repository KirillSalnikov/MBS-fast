#!/usr/bin/env python3
"""Relative output cutoffs must not change with absorption path length.

The two geometrically similar sizes use exactly the same rays and poses.
With a reference binary, also verify that absorption itself is unchanged
when output cutoffs are off, and that nonabsorbing results are unchanged.
"""
import argparse
import json
import math
import re
import subprocess
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--reference', type=Path)
    p.add_argument('--backend', choices=['cpu', 'cuda'], default='cpu')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    root = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    records = []

    def run(binary, label, radius, absorption, cutoffs):
        tag = '%s_r%g_k%g_%s' % (label, radius, absorption, cutoffs)
        folder = output/tag
        folder.mkdir(exist_ok=True)
        levels = ['.001', '.002', '.0001'] if cutoffs == 'production' else ['0']*3
        command = [str(binary.resolve()), '--method', 'po', '--backend', args.backend,
                   '--particle-file', str(root/'tests/data/cutoff_absorbing_20.dat'),
                   '--geometry', 'auto', '--k-eq', format(2*math.pi*radius/.34, '.17g'),
                   '--refractive-index', '1.465', str(absorption), '--wavelength-um', '.34',
                   '--max-reflections', '8', '--symmetry', '1', '1',
                   '--so3-rule', str(root/'tests/data/cutoff_absorbing_so3.txt'),
                   '--scattering-grid', '0', '180', '4', '180', '--threads', '2',
                   '--beam-cutoff-jones', levels[0], '--beam-cutoff-area', levels[1],
                   '--beam-cutoff-importance', levels[2], '--trace-max-beams', '20000',
                   '--trace-limit-retries', '0', '--trace-cache-mib', '0',
                   '--output', str(folder/'result'), '--close']
        (folder/'command.json').write_text(json.dumps(command, indent=2)+'\n')
        with (folder/'console.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=300)
        log = (folder/'console.log').read_text()
        def counter(name):
            match = re.search(re.escape(name)+r':\s*(\d+)', log)
            assert match, (name, folder)
            return int(match[1])
        assert counter('Hard tree-limit hits') == 0
        assert counter('Orientations still incomplete') == 0
        data = np.loadtxt(folder/'result/result.dat', skiprows=1)
        assert data.shape == (181, 18) and np.isfinite(data).all()
        matrix = data[:, -16:]
        assert (matrix[:, 0] > 0).all()
        kept = counter('Kept output beams')
        candidates = counter('Candidate output beams')
        record = dict(case=tag, candidates=candidates, kept=kept,
                      M11_180=float(matrix[-1, 0]), M22_180=float(matrix[-1, 5]))
        records.append(record)
        return matrix, (candidates, kept)

    for absorption in [0., .009653489]:
        selected = []
        for radius in [5., 10.]:
            for mode in ['production', 'off']:
                matrix, counts = run(args.binary, 'candidate', radius, absorption, mode)
                if mode == 'production':
                    assert 0 < counts[1] < counts[0], counts
                    selected.append(counts)
                else:
                    assert counts[1] == counts[0], counts
                if args.reference and (absorption == 0 or mode == 'off'):
                    reference, _ = run(args.reference, 'reference', radius, absorption, mode)
                    error = float(np.max(np.abs(matrix-reference)/reference[:, :1]))
                    assert error < 1e-11, (absorption, radius, mode, error)
        # Absorption in the retained field changes with radius. The relative
        # *selection* reference must not: this assertion fails before the fix.
        assert selected[0] == selected[1], (absorption, selected)
        print('PASS: size-invariant cutoff selection, k=', absorption, selected, flush=True)
    (output/'report.json').write_text(json.dumps(dict(passed=True, records=records), indent=2)+'\n')


if __name__ == '__main__':
    main()
