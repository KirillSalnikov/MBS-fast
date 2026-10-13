#!/usr/bin/env python3
"""Resident, restartable scattering benchmark; stdlib only (Python >=3.6).

Every job retains its command, environment, log, matrix and checked completion
record. Completed jobs survive disconnection/restart; a partial job is rerun.
"""
import argparse
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def save(path,data):
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2)+'\n');os.replace(str(tmp),str(path))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--backend',choices=['cpu','gpu'],default='gpu')
    p.add_argument('--gpu',type=int,default=0)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--cases',nargs='+')
    p.add_argument('--limit',type=int,default=0)
    p.add_argument('--shard-index',type=int,default=0)
    p.add_argument('--shard-count',type=int,default=1)
    p.add_argument('--wait-lock',action='store_true',help='Queue behind another benchmark using this device/backend.')
    a=p.parse_args();root=a.root.resolve();binary=a.binary.resolve()
    manifest=json.loads((root/'manifest.json').read_text());binary_sha=sha(binary)
    env=os.environ.copy();env['OMP_NUM_THREADS']=str(a.threads)
    if a.backend=='gpu':env['CUDA_VISIBLE_DEVICES']=str(a.gpu)
    lease=Path('/tmp')/('mbs_so3_'+str(os.getuid())+'_'+a.backend+'_'+str(a.gpu)+'.lock')
    with lease.open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX if a.wait_lock else fcntl.LOCK_EX|fcntl.LOCK_NB)
        status=root/('status_'+a.backend+'_'+str(a.gpu)+'.json')
        progress=dict(state='running',pid=os.getpid(),binary_sha256=binary_sha,
                      backend=a.backend,gpu=a.gpu,completed=0,reused=0,total=0)
        jobs=[j for j in manifest['jobs'] if not a.cases or j['case'] in a.cases]
        if a.shard_count<1 or not 0<=a.shard_index<a.shard_count:raise ValueError('invalid shard index/count')
        jobs=[job for index,job in enumerate(jobs) if index%a.shard_count==a.shard_index]
        if a.limit:jobs=jobs[:a.limit]
        progress['total']=len(jobs)
        child=[None]
        def stop(signum,frame):
            if child[0] is not None and child[0].poll() is None:
                child[0].terminate()
                try:child[0].wait(timeout=30)
                except subprocess.TimeoutExpired:child[0].kill();child[0].wait()
            progress.update(state='stopped',updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
            save(status,progress);raise SystemExit(128+signum)
        signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)
        for job in jobs:
            folder=root/'results'/a.backend/job['id'];folder.mkdir(parents=True,exist_ok=True)
            dest=folder/'matrix';matrix=dest/'matrix.dat'
            args=[str(binary),'--method','po','--backend',a.backend]+job['args']
            if a.backend=='cpu' and '--orientation-pipeline' in args:
                args.remove('--orientation-pipeline')
            args += ['--threads',str(a.threads),'--output',str(dest),'--close']
            if job['method']=='sobol':
                args+=['--haar-alpha','--sobol-seed',str(job['orientations']),str(job['seed'])]
            else:args+=['--so3-rule',str(root/job['rule'])]
            # Particle and theta paths are relative to the portable experiment.
            for flag in ['--particle-file','--theta-grid-file','--analytic-mean-cache']:
                if flag in args:
                    i=args.index(flag)+1;args[i]=str(root/args[i])
            signature=hashlib.sha256(json.dumps([args,binary_sha],sort_keys=True).encode()).hexdigest()
            complete=folder/'complete.json'
            progress.update(job=job['id'],updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
            save(status,progress)
            if complete.exists():
                previous=json.loads(complete.read_text())
                if previous['signature']!=signature:raise RuntimeError('Refusing incompatible cached job '+job['id'])
                if matrix.exists() and sha(matrix)==previous['matrix_sha256']:
                    progress['reused']+=1;progress['completed']+=1;continue
                raise RuntimeError('Corrupt completed matrix '+job['id'])
            # The native executable adds '(1)' to an existing output directory.
            # Preserve interrupted attempts, then give it a fresh target so the
            # expected output path remains stable across queue restarts.
            attempt='attempt_'+str(int(time.time()*1000000))+'_'
            for path in list(folder.iterdir()):
                if path.name.startswith('matrix') or path.name=='stdout.log':
                    path.rename(folder/(attempt+path.name))
            (folder/'command.args').write_text('\n'.join(args)+'\n')
            save(folder/'environment.json',{k:env.get(k) for k in ['CUDA_VISIBLE_DEVICES','OMP_NUM_THREADS','LD_LIBRARY_PATH']})
            started=time.monotonic()
            with (folder/'stdout.log').open('w') as log:
                child[0]=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT,env=env)
                code=child[0].wait();child[0]=None
            seconds=time.monotonic()-started
            if code or not matrix.exists():
                message='Job {} failed: {}'.format(job['id'],code) if code else 'No matrix emitted by '+job['id']
                progress.update(state='failed',error=message,updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
                save(status,progress);raise RuntimeError(message)
            rows=0
            with matrix.open() as stream:
                next(stream)
                for text in stream:
                    if not text.strip():continue
                    row=list(map(float,text.split()))
                    if len(row)!=18 or not all(math.isfinite(v) for v in row):
                        raise RuntimeError('Invalid Mueller matrix '+job['id'])
                    rows+=1
            if rows!=job['theta_rows']:raise RuntimeError('Wrong theta grid '+job['id'])
            record=dict(job=job,seconds=seconds,signature=signature,binary_sha256=binary_sha,
                        matrix_sha256=sha(matrix),theta_rows=rows)
            save(complete,record);progress['completed']+=1
            print('{} {:.6f}s {}/{}'.format(job['id'],seconds,progress['completed'],len(jobs)),flush=True)
        progress.update(state='complete',updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
        save(status,progress)


if __name__=='__main__':main()
