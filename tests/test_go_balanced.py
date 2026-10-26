#!/usr/bin/env python3
"""Bounded GO cutoff comparison, parallel/mirror regression and fullauto audit.
Fresh output directory required. No production run; uses Python stdlib only.
"""
import argparse
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import time


def read_matrix(path):
    rows = [list(map(float, s.split())) for s in path.read_text().splitlines()[1:] if s.strip()]
    assert len(rows) == 723 and all(math.isfinite(x) for r in rows for x in r)
    return rows


def difference(a, b):
    peak = max(r[2] for r in a)
    error = max(abs(x-y) for r,s in zip(a,b) for x,y in zip(r[2:],s[2:]))/peak
    positive = [(r,s) for r,s in zip(a,b) if r[2]>0]
    pointwise = max(abs(x-y)/r[2] for r,s in positive for x,y in zip(r[2:],s[2:]))
    significant = [(r,s) for r,s in positive if r[2]>peak*1e-4]
    resolved = max(abs(x-y)/r[2] for r,s in significant for x,y in zip(r[2:],s[2:]))
    def integral(rows):
        return sum(r[1]*r[2] for r in rows[:-2])+sum(r[2] for r in rows[-2:])
    return dict(all16_scaled_global_M11=error,
                all16_max_pointwise_scaled_M11=pointwise,
                all16_pointwise_above_1e_4_peak=resolved,
                relative_integral_change=(integral(b)-integral(a))/integral(a))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--count',type=int,default=256)
    parser.add_argument('--repeats',type=int,default=2)
    args=parser.parse_args()
    root=args.output.resolve();root.mkdir(parents=True,exist_ok=False)
    binary=args.binary.resolve()
    report={'binary':str(binary),'count':args.count,'repeats':args.repeats,
            'scope':'paired fixed Sobol sample, not a production accuracy certificate','cases':[]}
    def save(): (root/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    def run(command, name, expected=(0,)):
        folder=root/name;folder.mkdir()
        (folder/'command.json').write_text(json.dumps(list(map(str,command)),indent=2)+'\n')
        start=time.monotonic()
        with (folder/'run.log').open('w') as log:
            proc=subprocess.Popen(list(map(str,command)),stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            try:proc.wait(timeout=90)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid,signal.SIGKILL);proc.wait();raise
        elapsed=time.monotonic()-start
        assert proc.returncode in expected,(name,proc.returncode,(folder/'run.log').read_text()[-1500:])
        return elapsed
    def native(mu,depth,profile,name,threads=4,flags=('--mirror-gamma',),size=100):
        particle=['1',str(size),str(.7*size)] if mu==0 else ['10',str(size),str(.7*size),str(math.degrees(math.atan(mu/.7)))]
        out=root/name/'result'
        command=[binary,'--method','go','--backend','cpu','--particle',*particle,
                 '--ri','1.3116','0','--wavelength-um','.532','--max-reflections',str(depth),
                 '--symmetry','2','6','--sobol-seed',str(args.count),'11','--threads',str(threads),
                 '--scattering-grid','0','180','1','720','--go-fullauto-worker',
                 '--cutoff-profile',profile,'--trace-max-beams','0',*flags,'--close','--output',out]
        elapsed=run(command,name)
        return read_matrix(out/'result.dat'),elapsed
    for mu in (0,.2,.5,.8,.98,.999997035336085):
        for depth in (6,10):
            times={'off':[],'balanced':[]}; pair={}
            for repeat in range(args.repeats):
                for profile in (('off','balanced') if repeat%2==0 else ('balanced','off')):
                    name='mu_%g_n%d_%s_r%d'%(mu,depth,profile,repeat)
                    pair[profile],seconds=native(mu,depth,profile,name)
                    times[profile].append(seconds)
            delta=difference(pair['off'],pair['balanced'])
            assert delta['all16_scaled_global_M11']<.00125,delta
            assert abs(delta['relative_integral_change'])<.00125,delta
            log=(root/('mu_%g_n%d_balanced_r0'%(mu,depth))/'run.log').read_text()
            rejected=int(re.search(r'Rejected output beams: (\d+)',log)[1])
            record=dict(mu=mu,depth=depth,seconds=times,
                        speedup=statistics.median(times['off'])/statistics.median(times['balanced']),
                        rejected_output_beams=rejected,**delta)
            report['cases'].append(record);save()
            print(json.dumps(record),flush=True)
    assert sum(r['rejected_output_beams'] for r in report['cases'])>0
    single,_=native(.5,8,'balanced','one_thread',threads=1)
    parallel,_=native(.5,8,'balanced','four_threads')
    assert difference(single,parallel)['all16_scaled_global_M11']<1e-12
    pair,_=native(.5,8,'balanced','mirror_pair',flags=('--haar-mirror-audit',))
    assert difference(parallel,pair)['all16_scaled_global_M11']<2e-6
    scaled,_=native(.5,8,'balanced','size_double',size=200)
    for r in scaled:
        for j in range(2,len(r)): r[j]/=4
    assert difference(parallel,scaled)['all16_scaled_global_M11']<2e-6
    report['parallel_mirror_scaling_passed']=True;save()
    out=root/'controller_result'
    command=[binary,'--method','go','--backend','cpu','--particle','1','100','70',
             '--ri','1.3116','0','--wavelength-um','.532','--max-reflections','6',
             '--symmetry','2','6','--fullauto-mirror','on','--fullauto','.5',
             '--fullauto-initial','32','--max-orientations','512','--fullauto-max-rounds','5',
             '--fullauto-cpu-workers','4','--threads','1','--scattering-grid','0','180','1','2',
             '--cutoff-profile','balanced','--fullauto-cutoffs','off','--fullauto-resources','off',
             '--close','--output',out]
    run(command,'controller')
    state=json.loads((out/'fullauto_status.json').read_text())
    audit=json.loads((out/'physics_validation.json').read_text())
    assert state['status']=='converged' and state['selected_cutoff']=='balanced' and state['mirror_gamma']
    assert audit['passed'] and audit['cutoffs_checked'] and audit['cutoff_budget']==.5/8
    assert state['max_cutoff_residual95']==audit['cutoffs_off_residual95_scaled_M11']
    jobs=[p.read_text().splitlines() for p in (out/'jobs').glob('*/command.args')]
    assert jobs and all('--mirror-gamma' in j for j in jobs)
    assert set(j[j.index('--cutoff-profile')+1] for j in jobs)=={'off','balanced'}
    history=out/'fullauto_history.jsonl'
    def starts(): return sum(json.loads(s)['event']=='worker_start' for s in history.read_text().splitlines())
    before=starts();run(command,'resume');assert starts()==before
    for flag in ('--fullauto-final-checks','--fullauto-physics-audit'):
        run(command+[flag,'off'],'reject_'+flag.lstrip('-'),expected=(2,))
    changed=list(command);changed[changed.index('--cutoff-profile')+1]='off'
    run(changed,'reject_changed_profile',expected=(2,))
    automatic=list(command);pos=automatic.index('--max-reflections');del automatic[pos:pos+2]
    automatic+=['--fullauto-depth-start','2','--fullauto-depth-max','2','--fullauto-depth-pilot','16']
    automatic[automatic.index('--output')+1]=root/'automatic_result'
    run(automatic,'automatic',expected=(0,3,4))
    state=json.loads((root/'automatic_result/fullauto_status.json').read_text())
    assert state['selected_cutoff']=='balanced' and state['mirror_gamma']
    jobs=list((root/'automatic_result/jobs').glob('*/command.args'));assert jobs
    for job in jobs:
        argv=job.read_text().splitlines();assert '--mirror-gamma' in argv
    report['fullauto_audit']=audit;report['controller_resume_validation_passed']=True
    report['passed']=True;save();print('PASS',root/'summary.json',flush=True)

if __name__=='__main__':main()
