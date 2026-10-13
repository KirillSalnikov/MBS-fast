#!/usr/bin/env python3
"""Run one durable, serial GPU lane from an immutable campaign manifest.

Each job uses its own output directory. Restarting this runner replays an
interrupted job with the same executable and arguments; terminal failures are
recorded and the remaining independent jobs still run. Python 3.6 compatible.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(str(temporary), str(path))


def main():
    manifest_path = Path(sys.argv[1]).resolve()
    manifest = json.loads(manifest_path.read_text())
    root = manifest_path.parent
    binary = manifest['binary']
    digest = hashlib.sha256(Path(binary).read_bytes()).hexdigest()
    if digest != manifest['binary_sha256']:
        raise RuntimeError('Frozen executable changed; checkpoint replay refused')
    os.chdir(str(root))
    # The campaign does not inherit experimental overrides or GPU remapping.
    environment = dict(os.environ)
    for key in list(environment):
        if key.startswith('MBS_') or key == 'CUDA_VISIBLE_DEVICES':
            del environment[key]
    environment.update(manifest['environment'])
    statuses = {}
    state_file = root / 'campaign_status.json'
    if state_file.exists():
        statuses = json.loads(state_file.read_text()).get('jobs', {})
    started = time.time()

    def save(state, current=None):
        write_json(state_file, {'schema': 1, 'status': state,
                   'lane': manifest['lane'], 'gpus': manifest['gpus'],
                   'binary_sha256': digest, 'controller_runner_pid': os.getpid(),
                   'runner_started_unix': started, 'updated_unix': time.time(),
                   'current_job': current, 'jobs': statuses})

    save('running')
    for job in manifest['jobs']:
        identifier = job['id']
        previous = statuses.get(identifier, {})
        if previous.get('status') in ('converged', 'failed'):
            continue
        command = [binary] + job['arguments']
        output = Path(job['output'])
        output.parent.mkdir(parents=True, exist_ok=True)
        prior_seconds = previous.get('elapsed_seconds', 0)
        attempt = previous.get('attempt', 0) + 1
        begin = time.time()
        statuses[identifier] = {'status': 'running', 'attempt': attempt,
                                'started_unix': begin,
                                'elapsed_seconds': prior_seconds,
                                'output': str(output), 'command': command}
        save('running', identifier)
        print('START {} attempt={} gpus={}'.format(identifier, attempt,
                                                   manifest['gpus']), flush=True)
        with (root / (identifier + '.console.log')).open('a') as log:
            log.write('\nCAMPAIGN ATTEMPT {} {}\n'.format(attempt, begin))
            log.flush()
            process = subprocess.Popen(command, env=environment, stdout=log,
                                       stderr=subprocess.STDOUT)
            statuses[identifier]['pid'] = process.pid
            save('running', identifier)
            while process.poll() is None:
                time.sleep(10)
                statuses[identifier]['elapsed_seconds'] = prior_seconds + time.time() - begin
                save('running', identifier)
        state = {}
        adaptive_state = output / 'fullauto_status.json'
        if adaptive_state.exists():
            state = json.loads(adaptive_state.read_text())
        passed = process.returncode == 0 and state.get('status') == 'converged'
        statuses[identifier].update(status='converged' if passed else 'failed',
                                    returncode=process.returncode,
                                    finished_unix=time.time(),
                                    elapsed_seconds=prior_seconds + time.time() - begin,
                                    adaptive_state=state)
        save('running')
        print('END {} returncode={} status={}'.format(identifier,
              process.returncode, statuses[identifier]['status']), flush=True)
    failures = sum(x.get('status') == 'failed' for x in statuses.values())
    save('completed_with_failures' if failures else 'completed')


if __name__ == '__main__':
    main()
