#!/usr/bin/env python3
"""Changing analytic weight/cache files must retain exact prepared beams."""
import argparse
import json
import os
from pathlib import Path
import re
import socket
import struct
import subprocess
import time

p=argparse.ArgumentParser()
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--allow-misses',action='store_true')
a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True)
binary=str(a.binary.resolve());grid=root/'theta.csv'
grid.write_text('\n'.join(str(150+10*i/8)for i in range(9))+'\n')
base=['--method','po','--backend','cuda','--particle','1','794.3','196.2',
      '--ri','1.3116','0','--wavelength-um','.532','--symmetry','2','6',
      '--haar-alpha','--sobol-seed','64','23','--max-reflections','10',
      '--theta-grid-file',str(grid),'--phi-points','32','--threads','2',
      '--analytic-mean-threads','2','--orientation-chunk','32',
      '--cutoff-profile','balanced','--incoherent','--analytic-facet-average',
      '--analytic-shadow-control','circular','--no-text-reports','--close']
connection,peer=socket.socketpair();connection.settimeout(180)
worker=subprocess.Popen([binary,'--fullauto-worker-server',str(peer.fileno())],pass_fds=(peer.fileno(),))
peer.close();report=[]
def run(name,seed,weight,persistent):
    dest=root/name;dest.mkdir();weights=dest/'weights.tsv'
    weights.write_text('theta_deg reflection_weight shadow_weight\n'+
        ''.join('{} {} {}\n'.format(150+10*i/8,weight,weight/2)for i in range(9)))
    cmd=base[:];cmd[cmd.index('--sobol-seed')+2]=str(seed)
    cmd+=['--analytic-control-weights',str(weights),'--analytic-mean-cache',str(dest/'means.cache'),
          '--trace-cache-mib','256'if persistent else'0','--output',str(dest/'result')]
    (dest/'command.json').write_text(json.dumps([binary]+cmd,indent=2));log=dest/'stdout.log';start=time.time()
    if persistent:
        words=[str(log),'64','thread']+cmd
        connection.sendall(struct.pack('I',len(words))+b''.join(struct.pack('I',len(w.encode()))+w.encode()for w in words))
        response=b''
        while len(response)<4:
            part=connection.recv(4-len(response));assert part;response+=part
        assert struct.unpack('i',response)[0]==0,log.read_text()[-2000:]
    else:
        with log.open('w')as f:subprocess.run([binary]+cmd,stdout=f,stderr=f,check=True,timeout=180)
    seconds=time.time()-start;text=log.read_text();m=re.search(r'Prepared trace reuse: hits=(\d+) misses=(\d+)',text)
    rows=[]
    for line in(dest/'result/result.dat').read_text().splitlines():
        try:row=list(map(float,line.split()))
        except ValueError:continue
        if row:rows.append(row)
    return rows,int(m.group(1))if m else 0,seconds
try:
    run('warm',23,1,True)
    for name,seed,weight,want_hits in [('changed_weights',23,.25,True),('changed_seed',37,.5,False)]:
        reference,_,fresh_seconds=run(name+'_fresh',seed,weight,False)
        result,hits,seconds=run(name+'_cached',seed,weight,True)
        assert len(reference)==len(result)==9
        error=max(abs(x-y)/max(abs(u[2]),1e-100)for u,v in zip(reference,result)for x,y in zip(u[2:],v[2:]))
        assert error<1e-6,(name,error)
        if want_hits and not a.allow_misses:assert hits==64,(name,hits)
        if not want_hits:assert hits==0,(name,hits)
        report.append(dict(case=name,hits=hits,max_all16_scaled_M11=error,fresh_seconds=fresh_seconds,cached_seconds=seconds))
        (root/'comparison.json').write_text(json.dumps(report,indent=2));print(report[-1],flush=True)
finally:
    connection.close()
    try:worker.wait(timeout=15)
    except subprocess.TimeoutExpired:worker.kill();worker.wait()
