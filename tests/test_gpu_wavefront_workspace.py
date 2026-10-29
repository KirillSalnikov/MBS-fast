#!/usr/bin/env python3
"""Validate reused CUDA workspace across requests, optics modes and raw-cache hits."""
import argparse, json, os, re, socket, struct, subprocess
from pathlib import Path
from test_droxtal_resident import load, difference, integral_difference

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--tolerance',type=float,default=1e-7);a=p.parse_args()
    root=a.output.resolve();root.mkdir(parents=True);grid=root/'theta.txt';grid.write_text('0\n0.1\n170\n176\n180\n')
    env=os.environ.copy();env.update(MBS_GPU_TRACE_CONVEX_SCHEDULE='facets8',MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB='64',
        MBS_HOST_MEM_BUDGET_MB='2048',OMP_WAIT_POLICY='PASSIVE')
    common=['--method','po','--backend','cuda','--particle','4','5','--symmetry','2','6',
        '--ri','1.3116','0','--wavelength-um','.532','--haar-alpha','--threads','2','--orientation-pipeline',
        '--orientation-chunk','16','--theta-grid-file',str(grid),'--phi-points','8','--close',
        '--profile-phases','--allow-experimental-environment']
    variants=[('coherent',32,8,'off',[]),('coherent_repeat',32,8,'off',[]),
        ('incoherent',64,8,'balanced',['--incoherent']),
        ('absorbing',16,6,'balanced',['--absorption']),('concave',16,6,'off',[]),
        ('resized',16,6,'off',[]),('coherent_return',32,8,'off',[])]
    left,right=socket.socketpair();left.settimeout(180)
    with(root/'service.log').open('w')as log:
        worker=subprocess.Popen([str(a.binary),'--fullauto-worker-server',str(right.fileno())],pass_fds=[right.fileno()],
            env=env,stdout=log,stderr=log)
    right.close();records=[];hits=0;reuse=0
    try:
        for label,count,depth,cutoff,extra in variants:
            cmd=common+['--sobol-seed',str(count),'23','--max-reflections',str(depth),'--paired-depth',str(depth-2),
                '--orientation-prefix-output',str(count//2),'--cutoff-profile',cutoff]+extra
            if label=='resized':cmd[cmd.index('--particle')+2]='7'
            if label=='absorbing':cmd[cmd.index('--ri')+2]='.00042407'
            if label=='concave':
                at=cmd.index('--particle');del cmd[at:at+3]
                at=cmd.index('--symmetry');cmd[at+1:at+3]=['1','1']
                # The existing paired-depth feature supports convex geometry.
                at=cmd.index('--paired-depth');del cmd[at:at+2]
                cmd+=['--particle-file',str(Path(__file__).resolve().parents[1]/'examples/particles/concave_hexagonal.particle'),
                    '--resize-dmax-um','5']
            dirs={}
            for mode in ['cpu','resident']:
                dest=root/(label+'_'+mode);dest.mkdir();dirs[mode]=dest/'result';log=dest/'stdout.log'
                args=cmd+['--trace-cache-mib','64'if mode=='resident'else'0','--po-trace-backend',mode,'--output',str(dest/'result')]
                (dest/'command.args').write_text('\n'.join([str(a.binary)]+args)+'\n')
                if mode=='cpu':
                    with log.open('w')as stream:r=subprocess.run([str(a.binary)]+args,env=env,stdout=stream,stderr=stream,timeout=180)
                    assert r.returncode==0,(label,mode,log)
                else:
                    words=[str(log),'64','auto']+args
                    left.sendall(struct.pack('=I',len(words))+b''.join(struct.pack('=I',len(s.encode()))+s.encode()for s in words))
                    reply=b''
                    while len(reply)<4:
                        part=left.recv(4-len(reply));assert part,'worker exited';reply+=part
                    assert struct.unpack('=i',reply)[0]==0,(label,log)
                    text=log.read_text();accepted=sum(map(int,re.findall(r'GPU convex resident direct: accepted=(\d+)',text)))
                    if label in ['absorbing','concave']:
                        assert accepted==0,(label,accepted,'expected CPU fallback')
                    else:
                        assert accepted==count and 'CPU_fallback=0'in text,(label,accepted,count)
                    hits+=text.count('GPU raw trace reuse: hit')
                    reuse+=sum(map(int,re.findall(r'GPU convex workspace: reused_buffers=(\d+)',text)))
                    for size in re.findall(r'GPU convex workspace:.*retained_MiB=([0-9.]+)',text):assert float(size)<=64
            errors={}
            for f in dirs['cpu'].glob('result*.dat'):errors[f.name]=difference(load(f),load(dirs['resident']/f.name))
            for f in dirs['cpu'].glob('result*_integrals.tsv'):errors[f.name]=integral_difference(f,dirs['resident']/f.name)
            assert errors and max(errors.values())<=a.tolerance,(label,errors)
            records.append({'case':label,'errors':errors});print(records[-1],flush=True)
        assert hits>0 and reuse>0,(hits,reuse)
        left.sendall(struct.pack('=I',0));worker.wait(timeout=20);assert worker.returncode==0
        (root/'result.json').write_text(json.dumps({'passed':True,'records':records,'raw_cache_hits':hits,'reused_buffers':reuse},indent=2))
    finally:
        left.close()
        if worker.poll()is None:worker.terminate();worker.wait(timeout=20)

if __name__=='__main__':main()
