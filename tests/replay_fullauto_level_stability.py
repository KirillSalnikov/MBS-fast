"""Build the small native regression fixture from saved optical replicas.

No server access or optical calculation. The input directory is the archived
fullauto_back_validation_20261021/all16 directory containing saved_replicas.json,
saved_controls.json, and all16_replay.json. All 16 elements are corrected exactly
as Engine::Data, including the saved analytic polarization coefficients.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    raw = json.loads((args.source / 'saved_replicas.json').read_text())
    controls = json.loads((args.source / 'saved_controls.json').read_text())
    reference = json.loads((args.source / 'all16_replay.json').read_text())['results']
    weights = np.loadtxt(io.StringIO(next(
        value['text'] for name, value in controls['meta'].items()
        if name.endswith('/polarization_weights.tsv'))), skiprows=1)
    seeds = [11, 23, 37, 53, 71, 89, 107, 131]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w') as output:
        output.write(f'{len(reference)} 3 8 16\n')
        for expected in reference:
            theta = expected['theta']
            output.write(f'{theta:.17g}\n')
            for kind in ['low', 'high']:
                counts = expected[kind + '_counts']
                output.write(' '.join(map(str, counts)) + '\n')
                values = np.zeros((3, 8, 16))
                seen = np.zeros((3, 8), dtype=int)
                polar = np.array([
                    np.interp(theta, weights[:, 0], weights[:, j + 1])
                    for j in range(32)]).reshape(16, 2)
                for record in raw['records']:
                    level, n = record['level'], record['N']
                    if record['kind'] != kind or n != counts[level]:
                        continue
                    path = record['path']
                    rows = [row for row in raw['jobs'][path]['rows']
                            if abs(row['theta'] - theta) < 1e-10]
                    if not rows:
                        continue
                    assert len(rows) == 1
                    control = [row['controls'] for row in controls['jobs'][path]['rows']
                               if abs(row['theta'] - theta) < 1e-10]
                    assert len(control) == 1
                    native = np.array(rows[0]['matrix'], copy=True)
                    c = control[0]
                    native[1:] += (polar[1:, 0] * (c[1] - c[3])
                                   + polar[1:, 1] * (c[2] - c[4]))
                    seed = seeds.index(record['seed'])
                    values[level, seed] += record['sign'] * native
                    seen[level, seed] += 1
                assert np.all(seen[0] == 1) and np.all(seen[1:] == 2)
                combined = values.sum(axis=0)
                np.testing.assert_allclose(combined.mean(axis=0),
                                           expected[kind]['mean'], rtol=0, atol=1e-10)
                np.savetxt(output, values.reshape(-1, 16), fmt='%.17g')
                for field in ['mean', 'stat95_scaled_M11', 'ratio95']:
                    np.savetxt(output, [expected[kind][field]], fmt='%.17g')
    provenance = {
        'description': 'Saved prism L446.7 D147.1, backward PO, wavelength 0.532, m=1.3116.',
        'source_directory': str(args.source.resolve()),
        'source_sha256': {name: sha(args.source / name) for name in
                          ['saved_replicas.json', 'saved_controls.json', 'all16_replay.json']},
        'fixture_sha256': sha(args.output),
        'seeds': seeds,
        'level_depths': [8, 10, 12],
        'level_definition': ['M(n=8)', 'M(n=10)-M(n=8)', 'M(n=12)-M(n=10)'],
        'analytic_controls': 'Saved per-angle polarization weights applied before level subtraction.',
        'precision': '17 significant digits; eight complete replicas, all 16 Mueller entries.',
        'scope': ('Low/high saved estimators use different prefix lengths and mixed cutoff history. '
                  'They test acceptance of a changed correction, not cutoff-isolated bias or an '
                  'exact half-prefix optical identity. No new optical runs.'),
        'theta': raw['target'],
        'format': ('Header: angle_count, levels, seeds, elements. Per angle: theta; '
                   'then low and high: three counts, level-major seed-major matrix rows; '
                   'expected complete mean, scaled Student-t halfwidths, Fieller halfwidths.'),
    }
    args.output.with_suffix('.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(f'Wrote {args.output} ({args.output.stat().st_size} bytes), four real saved angles.')


if __name__ == '__main__':
    main()
