#!/usr/bin/env python3
"""Archive obsolete compact-layout publications after a stopped/completed run.
Never reads/deletes work blocks or recovery checkpoints. Dry run by default.
"""
import argparse
import contextlib
import fcntl
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import tarfile
import time


def digest(stream):
    h = hashlib.sha256()
    for block in iter(lambda: stream.read(1024 * 1024), b''):
        h.update(block)
    return h.hexdigest()


def inventory(directories, root):
    result = {}
    for directory in directories:
        if directory.is_symlink() or not directory.is_dir():
            raise ValueError('Publication must be a real directory: ' + str(directory))
        for p in directory.iterdir():
            if p.is_symlink() or not p.is_file():
                raise ValueError('Unexpected publication member: ' + str(p))
            with p.open('rb') as f:
                result[str(p.relative_to(root))] = (p.stat().st_size, digest(f))
    return result


@contextlib.contextmanager
def reader(path, compression):
    if compression == 'gzip':
        with tarfile.open(str(path), 'r:gz') as archive:
            yield archive
    else:
        process = subprocess.Popen(['zstd', '-q', '-d', '-c', str(path)], stdout=subprocess.PIPE)
        try:
            with tarfile.open(fileobj=process.stdout, mode='r|') as archive:
                yield archive
            process.stdout.close()
            if process.wait() != 0:
                raise RuntimeError('zstd decompression failed')
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()


def archive_results(root, keep=2, apply=False, compression='zstd'):
    root = Path(root).resolve()
    if (root / '.mbs-result-layout-v1').read_text() != 'mbs_compact_result_v1\n':
        raise ValueError('Not a compact fullauto result')
    # The controller holds this same lock; no archive/cleanup of a live run.
    with (root / 'work/controller.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        latest = (root / 'latest').resolve(strict=True)
        parent = root / 'results'
        if latest.parent != parent or not re.fullmatch(r'(snapshot_[0-9a-f]+|round_(?:[0-9]+|unknown)_[0-9a-f]+)', latest.name):
            raise ValueError('latest points outside owned publications')
        candidates = []
        for p in parent.iterdir():
            if re.fullmatch(r'(snapshot_[0-9a-f]+|round_(?:[0-9]+|unknown)_[0-9a-f]+)', p.name):
                if p.is_symlink() or not (p / 'manifest.fnv').is_file():
                    raise ValueError('Unrecognized publication: ' + str(p))
                candidates.append(p)
        candidates.sort(key=lambda p: p.stat().st_mtime_ns, reverse=True)
        retain = {latest}
        for p in candidates:
            if len(retain) >= max(2, keep):
                break
            retain.add(p)
        selected = [p for p in candidates if p not in retain]
        report = dict(selected=[p.name for p in selected], kept=sorted(p.name for p in retain), applied=False)
        if not apply or not selected:
            return report
        if compression == 'zstd' and not shutil.which('zstd'):
            raise RuntimeError('Install zstd or select --compression gzip; originals unchanged')
        before = inventory(selected, root)
        archive_dir = root / 'archive'
        archive_dir.mkdir(exist_ok=True)
        suffix = '.tar.zst' if compression == 'zstd' else '.tar.gz'
        target = archive_dir / ('snapshots-' + str(int(time.time() * 1000000)) + suffix)
        temporary = Path(str(target) + '.partial')
        if compression == 'gzip':
            with tarfile.open(str(temporary), 'x:gz', compresslevel=1) as tar:
                for p in selected:
                    tar.add(str(p), arcname=str(p.relative_to(root)))
        else:
            with temporary.open('xb') as f:
                process = subprocess.Popen(['zstd', '-q', '-1', '-T1', '-c'], stdin=subprocess.PIPE, stdout=f)
                try:
                    with tarfile.open(fileobj=process.stdin, mode='w|') as tar:
                        for p in selected:
                            tar.add(str(p), arcname=str(p.relative_to(root)))
                    process.stdin.close()
                    if process.wait() != 0:
                        raise RuntimeError('zstd compression failed')
                finally:
                    if process.poll() is None:
                        process.kill()
                        process.wait()
        found = {}
        with reader(temporary, compression) as tar:
            for member in tar:
                if member.isdir():
                    continue
                if not member.isfile() or member.name not in before or member.name in found:
                    raise ValueError('Unexpected archive member')
                with tar.extractfile(member) as f:
                    found[member.name] = (member.size, digest(f))
        if found != before or inventory(selected, root) != before:
            raise ValueError('Checksum mismatch: originals retained')
        with temporary.open('rb') as f:
            os.fsync(f.fileno())
            checksum = digest(f)
        os.replace(str(temporary), str(target))
        receipt = dict(sha256=checksum, files=before, restore='Extract this archive into ' + str(root))
        receipt_path = Path(str(target) + '.json')
        with receipt_path.open('x') as f:
            json.dump(receipt, f, indent=2)
            f.flush()
            os.fsync(f.fileno())
        fd = os.open(str(archive_dir), os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
        for p in selected:
            shutil.rmtree(str(p))
        fd = os.open(str(parent), os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
        report.update(applied=True, archive=str(target), files=len(before))
        return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result')
    parser.add_argument('--keep', type=int, default=2)
    parser.add_argument('--apply', action='store_true')
    parser.add_argument('--compression', choices=('zstd', 'gzip'), default='zstd')
    args = parser.parse_args()
    print(json.dumps(archive_results(args.result, args.keep, args.apply, args.compression), indent=2))
