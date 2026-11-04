"""Storage admission and out-of-filesystem diagnostics for serial queues."""
import errno
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time


def runtime_path(root):
    base = Path(os.environ.get('XDG_RUNTIME_DIR', '/dev/shm')) / ('mbs-fullauto-' + str(os.getuid()))
    base.mkdir(mode=0o700, parents=True, exist_ok=True)
    if base.stat().st_uid != os.getuid():
        raise RuntimeError('Runtime diagnostic directory has another owner')
    key = hashlib.sha256(str(Path(root).resolve()).encode()).hexdigest()[:20]
    return base / (key + '.json')


def record_runtime(root, state):
    state = dict(state, updated_unix=time.time(), runner_pid=os.getpid(), lane=str(root))
    path = runtime_path(root)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(state))
    os.replace(str(temp), str(path))


def check_storage(root, minimum_bytes=2*1024**3, minimum_inodes=250000):
    """Btrfs reports zero total inodes: that is not inode exhaustion."""
    root = Path(root)
    st = os.statvfs(str(root))
    result = dict(ready=False, free_bytes=st.f_bavail*st.f_frsize,
                  free_inodes=st.f_favail if st.f_files else None)
    if result['free_bytes'] < minimum_bytes:
        return dict(result, reason='byte reserve')
    if st.f_files and st.f_favail < minimum_inodes:
        return dict(result, reason='inode reserve')
    old = new = None
    try:
        fd, old = tempfile.mkstemp(prefix='.fullauto-storage-', dir=str(root))
        with os.fdopen(fd, 'wb') as stream:
            stream.write(b'fullauto storage admission\n')
            stream.flush()
            os.fsync(stream.fileno())
        new = old + '.renamed'
        os.replace(old, new)
        return dict(result, ready=True, reason='write/fsync/rename verified')
    except OSError as error:
        return dict(result, reason=str(error), errno=error.errno)
    finally:
        for name in (old, new):
            if name:
                try:
                    os.unlink(name)
                except FileNotFoundError:
                    pass


class StorageBlocked(RuntimeError):
    pass


def check_job_storage(lane, output=None):
    """Probe both write locations; the result may be on a different mount.

    For a new result use its nearest existing parent.  Probe distinct paths
    even on one device: a read-only bind mount can share the same st_dev.
    """
    locations = [('queue', Path(lane))]
    if output is not None:
        target = Path(output)
        if not target.is_absolute():
            target = Path(lane) / target
        while not target.exists() and target.parent != target:
            target = target.parent
        locations.append(('result', target))
    checks = []
    seen = set()
    for role, path in locations:
        key = str(path.absolute())
        if key in seen:
            continue
        seen.add(key)
        try:
            health = check_storage(path)
        except OSError as error:
            health = dict(ready=False, reason=str(error), errno=error.errno)
        health = dict(health, role=role, path=str(path))
        checks.append(health)
        if not health['ready']:
            return dict(health, checks=checks)
    return dict(ready=True, reason='queue and result storage verified', checks=checks)
