#!/usr/bin/env python3
"""Full controller convergence, integral receipts and computation-free backfill."""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    dest = root/'calculation'
    command = [str(args.binary.resolve()), '--method', 'po', '--backend', 'cpu',
        '--particle', '1', '.2', '.1', '--ri', '1.3116', '0', '--wavelength-um', '1',
        '--max-reflections', '2', '--threads', '1', '--fullauto', '.03',
        '--fullauto-resources', 'off', '--fullauto-cutoffs', 'off', '--cutoff-profile', 'off',
        '--symmetry', '1', '1', '--fullauto-mirror', 'off', '--phi-points', '4',
        '--fullauto-workers', 'persistent', '--fullauto-cpu-workers', '2',
        '--fullauto-pilot', '32', '--fullauto-initial', '256', '--fullauto-min-correction', '64',
        '--fullauto-orientation-block', '128', '--fullauto-theta-start', '5',
        '--fullauto-max-rounds', '16', '--max-orientations', '16384',
        '--scattering-grid', '0', '25', '4', '16', '--output', str(dest), '--close']
    (root/'command.args').write_text('\n'.join(command)+'\n')

    def run(label, expected=0):
        with (root/(label+'.log')).open('w') as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, timeout=300)
        assert result.returncode == expected, (label, (root/(label+'.log')).read_text()[-2000:])

    def starts():
        return sum(json.loads(line)['event'] == 'worker_start'
                   for line in (dest/'fullauto_history.jsonl').read_text().splitlines())

    run('initial')
    assert json.loads((dest/'fullauto_status.json').read_text())['status'] == 'converged'
    summary = dest/'mueller_fullauto_integrals.tsv'
    with summary.open() as stream:
        row = next(csv.DictReader(stream, delimiter='\t'))
    assert row['status'] == 'partial_angular_integral' and row['Cext'] == row['Csca']
    assert float(row['Cabs']) == 0 and row['full_angular_coverage'] == '0'
    marker = {p: p.read_bytes() for p in (dest/'jobs').glob('*/complete.txt')}
    saved = summary.read_bytes()
    before = starts()
    run('resume')
    assert starts() == before and summary.read_bytes() == saved
    snapshots = list((dest/'round_snapshots').glob('*/mueller_fullauto_integrals.tsv'))
    assert snapshots
    products = {line.split()[0] for line in (dest/'final_outputs.fnv').read_text().splitlines()
                if line.startswith(('mueller_fullauto_integrals.tsv ', 'fullauto_integrals.log ',
                                    'integral_confidence.tsv ', 'integral_replicas.tsv ',
                                    'integral_sources.fnv ', 'integral_estimator_sources.tsv ',
                                    'integral_summary.json '))}
    assert len(products) == 7
    summary.write_bytes(saved+b' ')
    run('corrupt_summary', 2)
    assert starts() == before
    summary.write_bytes(saved)
    # Simulate an older completed result and an interrupted first backfill:
    # one uncommitted stale TSV remains, but the final receipt has no integrals.
    receipt = dest/'final_outputs.fnv'
    receipt.write_text(''.join(line+'\n' for line in receipt.read_text().splitlines()
                               if line.split()[0] not in products))
    for name in products:
        (dest/name).unlink()
    summary.write_text('interrupted stale output\n')
    run('backfill')
    assert summary.read_bytes() == saved and starts() == before
    assert all(p.read_bytes() == value for p, value in marker.items())
    run('backfill_resume')
    assert starts() == before
    # The new source receipt protects scalar input files as well as matrices.
    paths = [line.rsplit('\t', 1)[0] for line in (dest/'integral_sources.fnv').read_text().splitlines()]
    source = next(Path(path) for path in paths if path.endswith('_integrals.tsv'))
    original = source.read_bytes()
    source.write_bytes(original+b' ')
    run('corrupt_source', 2)
    source.write_bytes(original)
    run('restored_source')
    assert starts() == before
    result = dict(converged=True, snapshots=len(snapshots), native_workers=before,
                  completed_blocks_preserved=len(marker), no_optical_work_on_resume=True,
                  interrupted_backfill_recovered=True, source_corruption_rejected=True,
                  summary_corruption_rejected=True)
    # GO must persist full-sphere finite-bin mass plus both integrated poles.
    dest = root/'go_calculation'
    command = [str(args.binary.resolve()), '--method', 'go', '--backend', 'cpu',
        '--particle', '1', '.2', '.1', '--ri', '1.3116', '.02', '--wavelength-um', '1',
        '--max-reflections', '2', '--threads', '1', '--fullauto', '.1',
        '--cutoff-profile', 'off', '--symmetry', '1', '1', '--fullauto-initial', '256',
        '--fullauto-workers', 'persistent', '--fullauto-cpu-workers', '2',
        '--fullauto-max-rounds', '16', '--max-orientations', '32768',
        '--scattering-grid', '0', '180', '1', '6', '--output', str(dest), '--close']
    (root/'go_command.args').write_text('\n'.join(command)+'\n')
    run('go_initial')
    with (dest/'mueller_fullauto_integrals.tsv').open() as stream:
        row = next(csv.DictReader(stream, delimiter='\t'))
    angular = sum(float(v[1])*float(v[2]) for v in
                  (line.split() for line in (dest/'mueller_fullauto.dat').read_text().splitlines()[1:]))
    for pole in ['forward', 'backward']:
        angular += float((dest/('mueller_'+pole+'.dat')).read_text().splitlines()[1].split()[2])
    assert row['full_angular_coverage'] == '1' and math.isclose(float(row['Csca']), angular, rel_tol=1e-12)
    assert math.isclose(float(row['Cext'])-float(row['Csca']), float(row['Cabs']), rel_tol=1e-12)
    before = starts()
    run('go_resume')
    assert starts() == before
    result['go_full_controller_poles_and_resume'] = True
    (root/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
