#!/usr/bin/env python3
"""Stdlib-only finite native test queue; no numerical parameter selection."""
import argparse,fcntl,json,os,subprocess,time
from pathlib import Path
def save(path,data):
    temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(data,indent=2)+'\n');os.replace(str(temporary),str(path))
def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',default='manifest.json');a=p.parse_args();root=a.root.resolve()
    jobs=json.loads((root/a.manifest).read_text())['jobs'];status=root/'acceptance_status.json'
    with open('/tmp/mbs_native_acceptance_'+str(os.getuid())+'.lock','a') as lease:
        fcntl.flock(lease,fcntl.LOCK_EX)
        for index,job in enumerate(jobs):
            record=root/(job['id']+'_complete.json')
            if record.exists():continue
            args=job['args'];(root/(job['id']+'.args')).write_text('\n'.join(args)+'\n')
            save(status,dict(status='running',completed=index,total=len(jobs),job=job['id'],pid=os.getpid()))
            started=time.monotonic()
            with (root/(job['id']+'.log')).open('w') as log:r=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT)
            if r.returncode:
                save(status,dict(status='failed',job=job['id'],returncode=r.returncode));raise RuntimeError(job['id'])
            save(record,dict(job=job['id'],seconds=time.monotonic()-started,args=args))
        save(status,dict(status='complete',completed=len(jobs),total=len(jobs)))
if __name__=='__main__':main()
