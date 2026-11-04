"""Exercise an actual fullauto output migration and unchanged resume identity."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=Path)
    parser.add_argument('--tool', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    result = root / 'calculation'
    command = [str(args.binary.resolve()), '--method', 'po', '--backend', 'cpu',
               '--particle', '1', '.2', '.1', '--ri', '1.3116', '0',
               '--wavelength-um', '.532', '--incoherent', '--max-reflections', '6',
               '--fullauto', '.5', '--threads', '1', '--fullauto-workers', 'persistent',
               '--fullauto-resources', 'off', '--fullauto-analytic-controls', 'off',
               '--fullauto-cutoffs', 'off', '--cutoff-profile', 'off',
               '--fullauto-mirror', 'off', '--fullauto-rotation-symmetry', 'off',
               '--fullauto-local-orientations', 'on', '--phi-points', '4',
               '--scattering-grid', '160', '170', '4', '8', '--fullauto-initial', '128',
               '--fullauto-pilot', '64', '--fullauto-kernel-count', '64',
               '--fullauto-min-correction', '64', '--fullauto-theta-start', '5',
               '--max-orientations', '512', '--fullauto-max-rounds', '1',
               '--fullauto-final-audit', 'off', '--fullauto-final-checks', 'off',
               '--output', str(result), '--close']
    environment = dict(os.environ)
    environment.pop('MBS_FULLAUTO_PACK_JOBS', None)
    with (root/'before.log').open('w') as log:
        before = subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT, timeout=180)
    assert before.returncode in (0, 3, 4), (root/'before.log').read_text()[-2000:]
    identity = (result/'configuration.txt').read_bytes()
    jobs = list((result/'jobs').glob('*/complete.txt'))
    assert jobs
    tool = [str(args.tool.resolve()), str(result)]
    with (result/'controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        rejected = subprocess.run(tool, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        assert rejected.returncode != 0 and 'controller is active' in rejected.stderr
        assert not (result/'job_pack_format.txt').exists(), 'active migration wrote policy'
    with (root/'migration.log').open('w') as log:
        subprocess.run(tool, stdout=log, stderr=subprocess.STDOUT, timeout=180, check=True)
    assert not list((result/'jobs').glob('*/complete.txt'))
    snapshots = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                 for path in (result/'round_snapshots').glob('*/mueller_evaluated.dat')}
    assert snapshots
    environment['MBS_FULLAUTO_PACK_JOBS'] = '1'
    with (root/'after.log').open('w') as log:
        after = subprocess.run(command, env=environment, stdout=log, stderr=subprocess.STDOUT, timeout=180)
    assert after.returncode in (0, 3, 4), (root/'after.log').read_text()[-2000:]
    assert (result/'configuration.txt').read_bytes() == identity
    assert all(hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest for path, digest in snapshots.items())
    assert 'resume inputs/binary/environment changed' not in (root/'after.log').read_text()
    summary = dict(passed=True, migrated_jobs=len(jobs), preserved_snapshots=len(snapshots),
                   old_identity_preserved=True, active_controller_migration_rejected=True,
                   return_codes=[before.returncode, after.returncode])
    (root/'verification.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
