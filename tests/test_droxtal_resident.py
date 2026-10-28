#!/usr/bin/env python3
"""All-16/ratio and prefix/lower-depth gates for the convex resident path.

Run on one CUDA device. CPU preparation and resident preparation use identical
orientations, cutoffs and diffraction precision. Failed native runs are never
accepted. The report preserves commands, matrices and measured phase timings.
"""
import argparse,csv,json,math,os,re,socket,struct,subprocess,time
from pathlib import Path

def load(path):
    return [[float(x)for x in line.split()]for line in path.read_text().splitlines()[1:]if line.strip()]

def difference(a,b):
    assert len(a)==len(b)and a
    worst=0.
    for x,y in zip(a,b):
        assert x[0]==y[0]and len(x)==len(y)==18
        assert all(math.isfinite(v)for v in x+y)and x[2]>0 and y[2]>0
        for i in range(2,18):
            worst=max(worst,abs(x[i]-y[i])/x[2])
            if i>2:worst=max(worst,abs(x[i]/x[2]-y[i]/y[2])/max(1.,abs(x[i]/x[2]),abs(y[i]/y[2])))
    return worst

def integral_difference(a,b):
    with a.open()as stream:x=list(csv.DictReader(stream,delimiter='\t'))
    with b.open()as stream:y=list(csv.DictReader(stream,delimiter='\t'))
    assert len(x)==len(y)and x
    worst=0.
    for u,v in zip(x,y):
        assert u['method']==v['method']and u['label']==v['label']
        scale=max(abs(float(u['Cext'])),abs(float(u['Csca'])),1e-30)
        for key in ['Cext','Csca','Cabs','projected_area','angular_scattering_integral']:
            a,b=float(u[key]),float(v[key])
            assert math.isfinite(a)and math.isfinite(b)
            worst=max(worst,abs(a-b)/max(scale,abs(a)))
    return worst

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--tolerance',type=float,default=1e-7);p.add_argument('--threads',type=int,default=2)
    p.add_argument('--backend',choices=['cpu','cuda'],default='cuda');p.add_argument('--warm-only',action='store_true')
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=True);binary=str(a.binary.resolve())
    grid=root/'angles.txt';grid.write_text('\n'.join(map(str,[0,.1,1,25,73,120,170,175,179,179.9,180]))+'\n')
    common=['--method','po','--backend',a.backend,'--wavelength-um','.532','--ri','1.3116','0',
        '--threads',str(a.threads),'--haar-alpha','--theta-grid-file',str(grid),'--phi-points','8',
        '--orientation-pipeline','--profile-phases','--trace-cache-mib','0','--close','--full-only']
    if a.backend=='cpu':common.remove('--orientation-pipeline')
    def args(shape,n,depth,seed,profile='balanced'):
        return common+shape+['--max-reflections',str(depth),'--sobol-seed',str(n),str(seed),'--cutoff-profile',profile]
    def run(label,cmd,service=None):
        out=root/label/'result';out.parent.mkdir(parents=True,exist_ok=True);log=out.parent/'stdout.log'
        cmd=cmd+['--output',str(out)];(out.parent/'command.args').write_text('\n'.join([binary]+cmd)+'\n')
        start=time.monotonic()
        if service is None:
            with log.open('w')as stream:r=subprocess.run([binary]+cmd,stdout=stream,stderr=subprocess.STDOUT,timeout=240)
            assert r.returncode==0,(label,r.returncode,log)
        else:
            words=[str(log),'64','auto']+cmd
            service.sendall(struct.pack('=I',len(words))+b''.join(struct.pack('=I',len(s.encode()))+s.encode()for s in words))
            reply=b''
            while len(reply)<4:
                b=service.recv(4-len(reply));assert b,(label,'worker exited');reply+=b
            assert struct.unpack('=i',reply)[0]==0,(label,log)
        phase={}
        for line in (out/'result_phase_timings.tsv').read_text().splitlines()[1:]:
            k,v,_=line.split();phase[k]=float(v)
        return dict(label=label,output=str(out),seconds=time.monotonic()-start,phase=phase,log=log.read_text())
    records=[]
    if not a.warm_only:
        repo=Path(__file__).resolve().parents[1]
        cases=[('droxtal5',['--particle','4','2.5','--symmetry','2','6'],128,6,[],True),
            ('droxtal79',['--particle','4','39.715','--symmetry','2','6'],256,14,[],True),
            ('droxtal100',['--particle','4','50','--symmetry','2','6'],128,18,[],True),
            ('droxtal_paired',['--particle','4','35.395','--symmetry','2','6'],256,14,['--paired-depth','12','--orientation-prefix-output','128'],True),
            ('prism',['--particle','1','398.1','138.9','--symmetry','2','6'],256,8,[],True),
            ('asymmetric',['--particle-file',str(repo/'tests/data/resident_asymmetric.particle'),'--resize-dmax-um','25','--symmetry','1','1'],128,8,[],True),
            ('concave',['--particle-file',str(repo/'examples/particles/concave_hexagonal.particle'),'--resize-dmax-um','5'],64,6,[],False),
            ('incoherent',['--particle','1','398.1','138.9','--symmetry','2','6'],128,8,['--incoherent'],True),
            ('incoherent_droxtal',['--particle','4','39.715','--symmetry','2','6'],128,14,['--incoherent'],True),
            ('absorbing',['--particle','4','5','--symmetry','2','6'],64,6,['--absorption'],False)]
        for name,shape,n,depth,extra,expected in cases:
            command=args(shape,n,depth,11)+extra
            if name=='incoherent_droxtal':command.remove('--full-only')
            if name=='absorbing':command[command.index('--ri')+2]='.00042407'
            pair=[run(name+'_'+mode,command+['--po-trace-backend',mode])for mode in ['cpu','resident']]
            errors={}
            products=sorted(Path(pair[0]['output']).glob('result*.dat'))
            for f in products:
                other=Path(pair[1]['output'])/f.name
                assert other.exists(),other
                errors[f.name]=difference(load(f),load(other))
            errors['integrals']=integral_difference(Path(pair[0]['output'])/'result_integrals.tsv',Path(pair[1]['output'])/'result_integrals.tsv')
            accepted=sum(map(int,re.findall(r'GPU convex resident direct: accepted=(\d+)',pair[1]['log'])))
            assert accepted==(n if expected else 0),(name,'resident coverage',accepted,n,expected)
            assert max(errors.values())<=a.tolerance,(name,errors,a.tolerance)
            result=dict(case=name,errors=errors,resident_orientations=accepted,
                processing_speedup=pair[0]['phase']['processing']/pair[1]['phase']['processing'],
                cpu=pair[0]['phase'],resident=pair[1]['phase'])
            records.append(result);print(json.dumps(result),flush=True)
    # Repeated requests in one worker: exercise learned memory, a lower-depth
    # field, and the exact prefix boundary. Compare to a cold standalone run.
    command=args(['--particle','4','5','--symmetry','2','6'],512,8,23,'off')
    command+=['--paired-depth','6','--orientation-prefix-output','256']
    if a.backend=='cuda':command+=['--po-trace-backend','cpu']
    cold=run('warm_reference',command)
    left,right=socket.socketpair();left.settimeout(240)
    with (root/'worker_service.log').open('w')as stream:
        worker=subprocess.Popen([binary,'--fullauto-worker-server',str(right.fileno())],pass_fds=[right.fileno()],stdout=stream,stderr=stream)
    right.close()
    try:
        first=run('worker_first',command,left);second=run('worker_warm',command,left)
        assert 'Prepared memory reuse:'in second['log']
        errs={}
        for f in Path(cold['output']).glob('result*.dat'):
            errs[f.name]=difference(load(f),load(Path(second['output'])/f.name))
        assert max(errs.values())<=a.tolerance,errs
        warm=dict(case='persistent_memory',errors=errs,cold_chunks=first['phase']['pipeline_chunks'],
            warm_chunks=second['phase']['pipeline_chunks'],cold_processing=first['phase']['processing'],
            warm_processing=second['phase']['processing'])
        assert warm['warm_chunks']<warm['cold_chunks'],warm
        records.append(warm);print(json.dumps(warm),flush=True)
        left.sendall(struct.pack('=I',0));worker.wait(timeout=20)
    finally:
        left.close()
        if worker.poll()is None:worker.terminate();worker.wait(timeout=20)
    (root/'result.json').write_text(json.dumps(dict(passed=True,tolerance=a.tolerance,records=records),indent=2))

if __name__=='__main__':main()
