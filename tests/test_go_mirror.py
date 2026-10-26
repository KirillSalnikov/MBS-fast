#!/usr/bin/env python3
"""Bounded GO mirror tests: all 16 entries, poles, threads, prefixes and controller.

The independent paired audit traces both partners instead of applying parity.
No production campaign is started. Output folder must be new.
"""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess

def matrix(path):
    rows=[list(map(float,line.split())) for line in path.read_text().splitlines()[1:] if line.strip()]
    assert rows and all(math.isfinite(v) for r in rows for v in r),path
    return rows

def difference(a,b):
    assert len(a)==len(b)
    scale=max(abs(r[2]) for r in b)
    assert scale>0
    return max(abs(x-y) for r,s in zip(a,b) for x,y in zip(r[2:],s[2:]))/scale

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--old-binary',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();binary=args.binary.resolve();root=args.output.resolve()
    root.mkdir(parents=True,exist_ok=False)
    records=[]
    def run(cmd,name,expected=(0,)):
        folder=root/name;folder.mkdir()
        (folder/'command.json').write_text(json.dumps(list(map(str,cmd)),indent=2)+'\n')
        with (folder/'stdout.log').open('w') as log:
            p=subprocess.Popen(list(map(str,cmd)),stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:p.wait(timeout=45)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid,signal.SIGKILL);p.wait();raise
        assert p.returncode in expected,(name,p.returncode,(folder/'stdout.log').read_text()[-2000:])
        return p.returncode
    def native(particle,name,flags=(),count=128,offset=0,threads=1,seed=11,index=('1.3116','0'),sym=('2','6'),exe=None,depth=4):
        out=root/name/'result'
        cmd=[exe or binary,'--method','go','--backend','cpu','--particle',*particle,
             '--ri',*index,'--wavelength-um','.532','--max-reflections',str(depth),
             '--symmetry',*sym,'--sobol-seed',str(count),str(seed),'--sobol-offset',str(offset),
             '--threads',str(threads),'--scattering-grid','0','180','1','720',
             '--go-fullauto-worker','--cutoff-profile','off','--trace-max-beams','0',
             '--close','--output',out,*flags]
        run(cmd,name)
        return matrix(out/'result.dat')
    particles=[['1','100','70'],['1','20','70'],
               ['10','100','70',str(math.degrees(math.atan(.2/.7)))],
               ['10','100','70',str(math.degrees(math.atan(.8/.7)))],
               ['10','100','70','55.0079']]
    # Paired partners exercise the reflection law, not independent sampling noise.
    for k,particle in enumerate(particles):
        for seed in (11,97):
            tag='p%d_s%d'%(k,seed)
            reduced=native(particle,tag+'_mirror',['--mirror-gamma'],seed=seed)
            paired=native(particle,tag+'_paired',['--haar-mirror-audit'],seed=seed)
            error=difference(reduced,paired)
            assert error<2e-6,(tag,error)
            assert all(r[2+j]==0 for r in reduced for j in (2,3,6,7,8,9,12,13))
            assert len(reduced)==723 # 721 finite bins plus two integrated delta rows.
            records.append(dict(case=tag,all16_and_poles_scaled_global_M11=error))
    # Deeper paths and a full rotational domain remain supported.
    for k,(index,sym,depth) in enumerate(((('1.3116','0'),('1','1'),8),)):
        a=native(particles[3],'extra%d_mirror'%k,['--mirror-gamma'],count=64,index=index,sym=sym,depth=depth)
        b=native(particles[3],'extra%d_pair'%k,['--haar-mirror-audit'],count=64,index=index,sym=sym,depth=depth)
        error=difference(a,b);assert error<2e-6,error
        records.append(dict(case='extra%d'%k,all16_and_poles_scaled_global_M11=error))
    whole=native(particles[3],'threads1',['--mirror-gamma'])
    parallel=native(particles[3],'threads4',['--mirror-gamma'],threads=4)
    assert difference(parallel,whole)<1e-12
    parts=[native(particles[3],'part%d'%i,['--mirror-gamma'],count=32,offset=32*i) for i in range(4)]
    average=[[sum(p[r][j] for p in parts)/4 for j in range(18)] for r in range(len(whole))]
    assert difference(average,whole)<1e-12
    if args.old_binary:
        for k,particle in enumerate((particles[0],particles[3])):
            a=native(particle,'off_new%d'%k)
            b=native(particle,'off_old%d'%k,exe=args.old_binary.resolve())
            assert difference(a,b)<1e-12,(k,difference(a,b))
    # Fullauto must propagate mirror through workers and preserve exact cache identity.
    out=root/'controller_result'
    cmd=[binary,'--method','go','--backend','cpu','--particle',*particles[3],
         '--ri','1.3116','0','--wavelength-um','.532','--symmetry','2','6',
         '--max-reflections','4','--scattering-grid','0','180','1','36',
         '--fullauto','.01','--fullauto-initial','64','--max-orientations','64',
         '--fullauto-max-rounds','1','--fullauto-cpu-workers','2','--threads','2',
         '--fullauto-mirror','on','--fullauto-work-queue','on',
         '--cutoff-profile','off','--fullauto-cutoffs','off','--close','--output',out]
    run(cmd,'controller',(3,))
    state=json.loads((out/'fullauto_status.json').read_text());assert state['mirror_gamma']
    jobs=list((out/'jobs').glob('*/command.args'));assert len(jobs)>=8
    for path in jobs:assert '--mirror-gamma' in path.read_text().splitlines()
    def starts():
        return sum(json.loads(s)['event']=='worker_start' for s in (out/'fullauto_history.jsonl').read_text().splitlines())
    before=starts();run(cmd,'resume',(3,));assert starts()==before
    changed=list(cmd);changed[changed.index('--fullauto-mirror')+1]='off'
    run(changed,'reject_changed_symmetry',(2,))
    automatic=list(cmd);at=automatic.index('--max-reflections');del automatic[at:at+2]
    automatic+=['--fullauto-depth-start','2','--fullauto-depth-max','2','--fullauto-depth-pilot','16']
    automatic[automatic.index('--output')+1]=root/'automatic_depth_result'
    run(automatic,'automatic_depth',(0,3,4))
    auto_jobs=list((root/'automatic_depth_result/jobs').glob('*/command.args'));assert auto_jobs
    for path in auto_jobs:assert '--mirror-gamma' in path.read_text().splitlines()
    auto_state=json.loads((root/'automatic_depth_result/fullauto_status.json').read_text())
    assert auto_state['mirror_gamma']
    # A shape without proven symmetry must not silently enable reconstruction.
    bad=list(cmd);at=bad.index('--particle');bad[at+1:at+5]=['2','100','70']
    bad[bad.index('--output')+1]=root/'unsupported_result'
    run(bad,'reject_unsupported',(2,))
    bad[bad.index('--fullauto-mirror')+1]='auto'
    run(bad,'auto_fallback',(3,))
    fallback=json.loads((root/'unsupported_result/fullauto_status.json').read_text())
    assert not fallback['mirror_gamma']
    bad=list(cmd);bad[bad.index('--symmetry')+2]='12';bad[bad.index('--output')+1]=root/'invalid_domain_result'
    run(bad,'reject_invalid_domain',(2,))
    bad=list(cmd);bad[bad.index('--ri')+2]='.001';bad[bad.index('--output')+1]=root/'absorbing_result'
    run(bad,'reject_unverified_absorption',(2,))
    report=dict(passed=True,paired_cases=records,threads_and_prefixes_passed=True,
                fullauto_propagation_and_resume_passed=True,automatic_depth_mirror_passed=True,
                unsupported_geometry_rejected=True,
                mirror_off_regression_passed=bool(args.old_binary))
    (root/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
