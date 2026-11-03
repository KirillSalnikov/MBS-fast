#!/usr/bin/env python3
"""Durable serial queue with per-job manifest reload and verified process adoption.

Python 3.6 compatible. An adoption receipt proves the PID, boot, start time,
command and binary of the already running calculation. It never restarts it.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from fullauto_storage import check_storage, record_runtime, StorageBlocked


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(str(temporary), str(path))


def process_info(pid):
    try:
        root = Path('/proc') / str(pid)
        fields = (root / 'stat').read_text().rsplit(') ', 1)[1].split()
        return {'state': fields[0], 'start_ticks': fields[19],
                'argv': (root / 'cmdline').read_bytes().decode().split('\0')[:-1]}
    except (FileNotFoundError, ProcessLookupError):
        return None


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    manifest_path = Path(sys.argv[1]).resolve()
    root = manifest_path.parent
    os.chdir(str(root))
    lock = (root / 'campaign_runner.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    manifest = json.loads(manifest_path.read_text())
    infrastructure = {k: manifest[k] for k in ['lane', 'gpus', 'environment']}
    state_file = root / 'campaign_status.json'
    statuses = json.loads(state_file.read_text()).get('jobs', {}) if state_file.exists() else {}
    receipt_path = root / 'campaign_handoff.json'
    handoff = json.loads(receipt_path.read_text()) if receipt_path.exists() else None
    started = time.time()
    poll_seconds = float(manifest.get('poll_seconds', 5))
    if not 0.05 <= poll_seconds <= 30:
        raise RuntimeError('Invalid queue poll interval')

    last_storage_check = [0.0]

    def save(state, current=None):
        if time.time() - last_storage_check[0] >= 30:
            health = check_storage(root)
            last_storage_check[0] = time.time()
            if not health['ready']:
                record_runtime(root, dict(status='storage_blocked', current_job=current, storage=health))
                raise StorageBlocked(health['reason'])
        try:
            return save_on_disk(state, current)
        except OSError as error:
            record_runtime(root, dict(status='storage_blocked', current_job=current,
                                     errno=error.errno, reason=str(error)))
            raise

    def save_on_disk(state, current=None):
        write_json(state_file, {'schema': 2, 'status': state,
                   'lane': manifest['lane'], 'gpus': manifest['gpus'],
                   'binary_sha256': manifest['binary_sha256'],
                   'controller_runner_pid': os.getpid(),
                   'runner_started_unix': started, 'updated_unix': time.time(),
                   'current_job': current, 'jobs': statuses})

    while True:
        # An edit changes only work which has not started. This process keeps
        # the current job's command and environment fixed until it exits.
        manifest = json.loads(manifest_path.read_text())
        if any(manifest[k] != v for k, v in infrastructure.items()):
            raise RuntimeError('Live queue infrastructure changed')
        # A held failure remains held after the runner itself is restarted.
        job = next((j for j in manifest['jobs'] if
                    statuses.get(j['id'], {}).get('status') != 'converged' and
                    (statuses.get(j['id'], {}).get('status') != 'failed' or
                     j.get('hold_queue_on_failure', False) or
                     statuses.get(j['id'], {}).get('handoff_recovery'))), None)
        if job is None:
            break
        if statuses.get(job['id'], {}).get('status') == 'failed':
            save('needs_attention', job['id'])
            print('HOLD {}: failed calculation awaits an explicit resume'.format(job['id']), flush=True)
            return
        health = check_storage(root)
        if not health['ready']:
            record_runtime(root, dict(status='storage_blocked', current_job=job['id'], storage=health))
            print('WAIT storage: ' + health['reason'], flush=True)
            time.sleep(30)
            continue
        record_runtime(root, dict(status='storage_ready', current_job=job['id'], storage=health))
        identifier = job['id']
        previous = statuses.get(identifier, {})
        binary = job.get('binary', manifest['binary'])
        expected = job.get('binary_sha256', manifest['binary_sha256'])
        if digest(binary) != expected:
            raise RuntimeError('Job executable changed; checkpoint replay refused')
        command = [binary] + job['arguments']
        if previous.get('command') and previous['command'] != command:
            raise RuntimeError('Interrupted job command changed; checkpoint replay refused')
        output = Path(job['output'])
        prior_seconds = previous.get('elapsed_seconds', 0)
        begin = time.time()
        token = handoff if handoff and handoff['job'] == identifier else None
        pid = previous.get('pid')
        info = process_info(pid) if pid else None
        # stale-handoff-compatible-resume-v1
        adopting = False
        if token is not None:
            alive = bool(info and info['state'] != 'Z')
            identity_ok = (isinstance(pid, int) and pid > 0 and token.get('pid') == pid
                           and token.get('argv') == command and token.get('binary_sha256') == expected)
            if alive:
                identity_ok = (identity_ok
                               and token.get('boot_id') == Path('/proc/sys/kernel/random/boot_id').read_text().strip()
                               and info['start_ticks'] == token.get('start_ticks')
                               and info['argv'] == command
                               and digest('/proc/{}/exe'.format(pid)) == expected)
            if not identity_ok:
                statuses[identifier] = dict(previous, recovery_error='Process adoption identity mismatch')
                save('needs_attention', identifier)
                print('HOLD {}: adoption identity mismatch; no process started'.format(identifier), flush=True)
                return
            if alive:
                adopting = True
            else:
                # A receipt describes a former process, not a completed job.
                # Keep the exact manifest command and replay its own checkpoint.
                archive = root / 'handoff_archive'
                archive.mkdir(exist_ok=True)
                destination = archive / ('campaign_handoff_{}_{}_{}.json'.format(
                    identifier, int(time.time() * 1000000), os.getpid()))
                os.replace(str(receipt_path), str(destination))
                previous = dict(previous, status='pending_resume', pid=None,
                                handoff_recovery={'policy': 'stale-handoff-compatible-resume-v1',
                                                  'archive': str(destination), 'former_pid': pid,
                                                  'former_boot_id': token.get('boot_id'),
                                                  'recovered_unix': time.time()})
                statuses[identifier] = previous
                save('pending_resume', identifier)
                handoff = None
                token = None
                print('RESUME {}: archived dead-process receipt; same command and checkpoint'.format(identifier), flush=True)
        elif info and info['state'] != 'Z':
            raise RuntimeError('A live calculation needs an explicit adoption receipt')
        environment = dict(os.environ)
        for key in list(environment):
            if key.startswith('MBS_') or key == 'CUDA_VISIBLE_DEVICES':
                del environment[key]
        environment.update(manifest['environment'])
        if adopting:
            statuses[identifier] = dict(previous, status='running', adopted_running_process=True)
            save('running', identifier)
            write_json(root / 'campaign_handoff_verified.json',
                       dict(token, runner_pid=os.getpid(), verified_unix=time.time(),
                            process_alive=bool(info and info['state'] != 'Z')))
            print('ADOPT {} pid={} (no native restart)'.format(identifier, pid), flush=True)
            while info and info['state'] != 'Z':
                if info['start_ticks'] != token['start_ticks'] or info['argv'] != command:
                    raise RuntimeError('Adopted process identity changed')
                time.sleep(poll_seconds)
                statuses[identifier]['elapsed_seconds'] = prior_seconds + time.time() - begin
                save('running', identifier)
                info = process_info(pid)
            returncode = None  # We are not its parent; do not invent an exit code.
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            attempt_history = list(previous.get('attempt_history', []))
            if previous:
                attempt_history.append({k: v for k, v in previous.items() if k != 'attempt_history'})
            statuses[identifier] = dict(previous, status='running', attempt=previous.get('attempt', 0) + 1,
                                        started_unix=begin, elapsed_seconds=prior_seconds,
                                        output=str(output), command=command, binary_sha256=expected,
                                        attempt_history=attempt_history)
            for obsolete in ('pid', 'returncode', 'finished_unix', 'adaptive_state',
                             'adopted_running_process', 'exit_confirmation', 'recovery_error'):
                statuses[identifier].pop(obsolete, None)
            save('running', identifier)
            print('START {} binary={}'.format(identifier, binary), flush=True)
            with (root / (identifier + '.console.log')).open('a') as log:
                log.write('\nCAMPAIGN ATTEMPT {} {}\n'.format(statuses[identifier]['attempt'], begin));log.flush()
                process = subprocess.Popen(command, env=environment, stdout=log, stderr=subprocess.STDOUT)
                statuses[identifier]['pid'] = process.pid
                save('running', identifier)
                while process.poll() is None:
                    time.sleep(poll_seconds)
                    statuses[identifier]['elapsed_seconds'] = prior_seconds + time.time() - begin
                    save('running', identifier)
            returncode = process.returncode
        state = {}
        state_path = output / 'fullauto_status.json'
        if state_path.exists():
            state = json.loads(state_path.read_text())
        passed = state.get('status') == 'converged' and (returncode == 0 or adopting)
        if adopting:
            # Do not accept a pre-existing convergence record as a new result.
            passed = passed and state_path.stat().st_mtime >= token['prepared_unix']
            statuses[identifier]['exit_confirmation'] = 'solver_status_after_adoption'
        statuses[identifier].update(status='converged' if passed else 'failed', returncode=returncode,
                                    finished_unix=time.time(), elapsed_seconds=prior_seconds + time.time() - begin,
                                    adaptive_state=state)
        save('running')
        if adopting:
            os.replace(str(receipt_path), str(root / 'campaign_handoff_completed.json'))
            handoff = None
        print('END {} status={}'.format(identifier, statuses[identifier]['status']), flush=True)
        # fullauto_target_failure_hold_v1: only the explicitly marked job.
        # A normal return keeps systemd Restart=on-failure from retrying or
        # advancing after a failed deployment. The operator repairs its state.
        if not passed and (job.get('hold_queue_on_failure', False) or statuses[identifier].get('handoff_recovery')):
            save('needs_attention', identifier)
            print('HOLD {}: failed scoped refinement; following jobs untouched'.format(identifier), flush=True)
            return

    failures = sum(x.get('status') == 'failed' and statuses.get(x.get('retry_job_id'), {}).get('status') != 'converged'
                   for x in statuses.values())
    save('completed_with_failures' if failures else 'completed')


if __name__ == '__main__':
    main()
