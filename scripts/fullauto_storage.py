"""Storage admission and out-of-filesystem diagnostics for serial queues."""
import errno
import fcntl
import hashlib
import json
import os
import re
import struct
import subprocess
from pathlib import Path
import tempfile
import time
import uuid


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


# Admission margins, not guarantees of future allocation. Global reserve is
# already included in metadata/bytes_may_use; never subtract it a second time.
BTRFS_METADATA_MARGIN = 64 * 1024**2
BTRFS_ALLOCATION_CHUNK = 1024**3
BTRFS_USAGE_TIMEOUT = 3
BTRFS_PRESSURE_ATTEMPTS = 4
BTRFS_PRESSURE_RECHECK_SECONDS = 2
BTRFS_SYSFS = Path('/sys/fs/btrfs')
# Linux UAPI: _IOR(0x94, 31, struct btrfs_ioctl_fs_info_args[1024 bytes]).
BTRFS_IOC_FS_INFO = 0x8400941f


def filesystem_type(root):
    """Choose the deepest mount in this namespace, including bind mounts."""
    path = str(Path(root).resolve())
    best = ('', 'unknown')
    try:
        rows = Path('/proc/self/mountinfo').read_text().splitlines()
        for row in rows:
            left, right = row.split(' - ', 1)
            mount = re.sub(r'\\([0-7]{3})', lambda m: chr(int(m.group(1), 8)), left.split()[4])
            if (path == mount or path.startswith(mount.rstrip('/') + '/')) and len(mount) >= len(best[0]):
                best = (mount, right.split()[0])
    except (OSError, ValueError, IndexError):
        return 'unknown'
    return best[1]


def _btrfs_identity(root):
    fd = os.open(str(root), os.O_RDONLY | os.O_DIRECTORY)
    try:
        info = bytearray(1024)
        fcntl.ioctl(fd, BTRFS_IOC_FS_INFO, info, True)
        return str(uuid.UUID(bytes=bytes(info[16:32]))), struct.unpack_from('=Q', info, 8)[0]
    finally:
        os.close(fd)


def _btrfs_usage(root):
    """Read-only, non-root diagnostic; no balance or device scan is requested."""
    executable = next((p for p in ('/usr/sbin/btrfs', '/sbin/btrfs', '/usr/bin/btrfs', '/bin/btrfs')
                       if os.path.isfile(p) and os.access(p, os.X_OK)), None)
    if executable is None:
        raise RuntimeError('btrfs-progs executable unavailable')
    env = dict(os.environ, LC_ALL='C', LANG='C')
    result = subprocess.run([executable, 'filesystem', 'usage', '-b', str(root)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, timeout=BTRFS_USAGE_TIMEOUT, env=env)
    if result.returncode:
        raise RuntimeError('btrfs usage failed: ' + result.stderr.strip()[:500])
    match = re.search(r'^\s*Device unallocated:\s*(\d+)\s*$', result.stdout, re.M)
    if not match:
        raise RuntimeError('btrfs usage did not report unallocated bytes')
    return dict(device_unallocated_bytes=int(match.group(1)),
                usage_warning=result.stderr.strip()[:500])


def btrfs_headroom(root):
    """Read allocator reservations, not df/data free space.

    Known Btrfs with unavailable allocator information is blocked explicitly.
    Healthy non-Btrfs never enters this check. Low existing metadata slack may
    pass only with physical unallocated room for any reservation deficit plus
    a conventional new chunk on a known single-device profile, followed by the
    usual write/fsync/rename probe. Aggregate free
    space alone cannot establish a multi-device placement. These are admission
    margins; fragmentation, new reservations and device errors can still fail.
    """
    health = dict(ready=False, source='btrfs sysfs allocation',
                  metadata_margin_bytes=BTRFS_METADATA_MARGIN)
    try:
        fsid, devices = _btrfs_identity(root)
        health.update(fsid=fsid, device_count=devices)
        allocation = BTRFS_SYSFS / fsid / 'allocation'
        def number(name, optional=False):
            try:
                value = int((allocation / name).read_text().strip())
            except FileNotFoundError:
                if optional:
                    return 0
                raise
            if value < 0:
                raise ValueError('negative Btrfs counter ' + name)
            return value
        names = ('total_bytes', 'bytes_used', 'bytes_may_use', 'bytes_reserved',
                 'bytes_pinned', 'bytes_readonly', 'disk_total', 'flags')
        meta = {name: number('metadata/' + name) for name in names}
        meta['bytes_zone_unusable'] = number('metadata/bytes_zone_unusable', True)
        health['metadata'] = meta
        health['global_rsv_reserved'] = number('global_rsv_reserved', True)
        health['global_rsv_size'] = number('global_rsv_size', True)
        if not meta['total_bytes'] or not meta['flags'] & 4:
            raise ValueError('invalid metadata space-info counters')
        slack = meta['total_bytes'] - sum(meta[name] for name in
                ('bytes_used', 'bytes_may_use', 'bytes_reserved', 'bytes_pinned',
                 'bytes_readonly', 'bytes_zone_unusable'))
        health['metadata_headroom_bytes'] = slack
        if slack >= BTRFS_METADATA_MARGIN:
            return dict(health, ready=True, reason='Btrfs existing metadata headroom')
        health.update(_btrfs_usage(root))
        # disk_total/total_bytes is a physical/logical ratio, not data-free space.
        ratio, remainder = divmod(meta['disk_total'], meta['total_bytes'])
        if devices != 1 or remainder or ratio not in (1, 2) or meta['flags'] != 4:
            return dict(health, failure_kind='metadata_pressure',
                        reason='Btrfs low metadata headroom; allocator placement unknown')
        # Btrfs permits metadata reservations against unallocated space, so
        # negative existing-chunk slack alone is not ENOSPC. Cover that whole
        # logical deficit plus one extra GiB, then account for physical copies.
        # Keep the measured negative slack visible; never substitute df space.
        deficit = max(0, -slack)
        required = (deficit + BTRFS_ALLOCATION_CHUNK) * ratio
        health['metadata_reservation_deficit_bytes'] = deficit
        health['physical_metadata_ratio'] = ratio
        health['allocation_margin_bytes'] = required
        if health['device_unallocated_bytes'] < required:
            return dict(health, failure_kind='metadata_pressure',
                        reason='Btrfs metadata reserve: insufficient metadata and unallocated headroom')
        return dict(health, ready=True,
                    reason='Btrfs low metadata headroom; single-device allocation margin available',
                    warning='Admission margin only; chunk placement and future allocation are not guaranteed')
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        return dict(health, reason='Btrfs allocation information unavailable: ' + str(error)[:500],
                    diagnostic='Read /sys/fs/btrfs/<fsid>/allocation and run btrfs filesystem usage -b on this path')


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
    result['filesystem'] = filesystem_type(root)
    if result['filesystem'] == 'btrfs':
        result['btrfs'] = btrfs_headroom(root)
        if not result['btrfs']['ready']:
            return dict(result, reason=result['btrfs']['reason'])
    elif result['filesystem'] == 'unknown':
        result['filesystem_warning'] = 'Mount type unavailable; existing byte/inode/write checks remain in effect'
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


def check_job_storage_confirmed(lane, output=None):
    """Confirm Btrfs reservation pressure; never retry a real I/O failure.

    Sysfs exposes allocator counters through separate reads. A negative sum is
    retained as evidence, not clamped or assumed to prove a filesystem error.
    Four consecutive low samples wait at most six seconds. If committed space
    is still pinned by a transaction, attempt one bounded filesystem sync and
    recheck the same margins. Real I/O errors and failed syncs fail closed.
    """
    samples = []
    waited = 0
    for attempt in range(BTRFS_PRESSURE_ATTEMPTS):
        health = check_job_storage(lane, output)
        samples.append(dict(health))
        allocator_only = (not health['ready'] and 'errno' not in health
                          and health.get('btrfs', {}).get('failure_kind') == 'metadata_pressure')
        if health['ready'] or not allocator_only or attempt + 1 == BTRFS_PRESSURE_ATTEMPTS:
            drain = None
            if (allocator_only and attempt + 1 == BTRFS_PRESSURE_ATTEMPTS
                    and health.get('btrfs', {}).get('metadata', {}).get('bytes_pinned', 0) > 0):
                drain_path = health['path']
                drain = drain_btrfs_transaction(drain_path)
                if drain['success']:
                    health = check_job_storage(lane, output)
                    samples.append(dict(health))
                    allocator_only = (not health['ready'] and 'errno' not in health
                                      and health.get('btrfs', {}).get('failure_kind') == 'metadata_pressure')
                else:
                    health = dict(health, reason='Btrfs transaction sync failed: ' + drain['error'])
                post = health.get('btrfs', {})
                for check in health.get('checks', []):
                    if check.get('path') == drain_path:
                        post = check.get('btrfs', post)
                diagnostic = dict(policy='bounded-btrfs-transaction-sync-v1', path=drain_path,
                                  success=drain['success'], seconds=drain.get('seconds'),
                                  child_pending=drain.get('child_pending', False),
                                  post_ready=health['ready'],
                                  post_metadata_headroom_bytes=post.get('metadata_headroom_bytes'))
                try:
                    os.write(2, ('STORAGE transaction_sync ' + json.dumps(diagnostic, sort_keys=True) + '\n').encode())
                except OSError:
                    pass  # Diagnostics cannot replace the original storage decision.
            return dict(health, storage_confirmation=dict(
                policy='bounded-btrfs-pressure-confirmation-v2',
                attempts=len(samples), wait_seconds=waited, samples=samples,
                transaction_sync=drain,
                persistent_allocator_pressure=bool(allocator_only)))
        time.sleep(BTRFS_PRESSURE_RECHECK_SECONDS)
        waited += BTRFS_PRESSURE_RECHECK_SECONDS


def drain_btrfs_transaction(path):
    """Flush this filesystem once; never balance, remount, or lower a margin.

    The runner remains able to stop its owned calculation even if syncfs hangs:
    a killed uninterruptible sync child is not waited for without a timeout.
    """
    started = time.monotonic()
    process = None
    try:
        executable = next(p for p in ('/usr/bin/sync', '/bin/sync')
                          if os.path.isfile(p) and os.access(p, os.X_OK))
        process = subprocess.Popen([executable, '-f', str(path)],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                   close_fds=True)
        code = process.wait(timeout=10)
        if code != 0:
            raise RuntimeError('sync -f exit ' + str(code))
        return dict(success=True, path=str(path), seconds=time.monotonic()-started)
    except (StopIteration, OSError, RuntimeError, subprocess.TimeoutExpired) as error:
        pending = False
        if process is not None and process.poll() is None:
            try:
                process.kill()
                process.wait(timeout=1)
            except subprocess.TimeoutExpired:
                pending = True
            except OSError:
                pending = process.poll() is None
        return dict(success=False, path=str(path), seconds=time.monotonic()-started,
                    error=str(error) or 'sync executable unavailable', child_pending=pending)


def read_runtime_record(root):
    """A queue-file write may fail; preserve its hold across service restarts."""
    path = runtime_path(root)
    try:
        if path.stat().st_size > 1024 * 1024:
            raise RuntimeError('Runtime storage diagnostic exceeds its size bound')
        value = json.loads(path.read_text())
    except FileNotFoundError:
        return None
    if value.get('lane') != str(root):
        raise RuntimeError('Runtime storage diagnostic belongs to another lane')
    return value
