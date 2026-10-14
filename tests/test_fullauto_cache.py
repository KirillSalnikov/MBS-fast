#!/usr/bin/env python3
"""Compare frozen estimators before/after cache indexing, including recovery."""
import argparse,json,os,subprocess,time
from pathlib import Path
import numpy as np
from test_fullauto_arbitrary import asymmetric_tetrahedron

def main():
    p=argparse.ArgumentParser();p.add_argument('--before',type=Path,required=True);p.add_argument('--after',type=Path,required=True)
    p.add_argument('--worker',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=True);tet=out/'asymmetric.particle';asymmetric_tetrahedron(tet)
    reports=[]
    for mode in ['coherent','incoherent','paired','controls','asymmetric']:
        pair=[];times={}
        for label,probe in [('before',a.before),('after',a.after)]:
            dest=out/(mode+'_'+label);shape=['--particle','1','.2','.1'] if mode!='asymmetric' else ['--particle-file',str(tet),'--resize-dmax-um','.2']
            cmd=[str(probe.resolve()),'--method','po','--backend','cpu']+shape+['--ri','1.3116','0','--wavelength-um','1','--max-reflections','4',
                '--threads','1','--fullauto','.02','--fullauto-analytic-controls','on' if mode=='controls' else 'off','--fullauto-mirror','off',
                '--fullauto-resources','off','--cutoff-profile','off','--phi-points','4','--scattering-grid','0','180','4','8','--fullauto-initial','256',
                '--fullauto-orientation-block','128','--fullauto-depth-sharing','on' if mode=='paired' else 'off','--fullauto-work-queue','on','--output',str(dest)]
            if mode=='incoherent':cmd+=['--incoherent']
            env=dict(os.environ,TASK_MBS_WORKER=str(a.worker.resolve()))
            (out/(dest.name+'.args')).write_text('\n'.join(cmd)+'\n');start=time.monotonic()
            with (out/(dest.name+'.log')).open('w')as log:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
            times[label]=time.monotonic()-start;pair.append(np.loadtxt(dest/'frozen_matrix.dat',skiprows=1))
            if label=='after':
                markers={str(f):f.stat().st_mtime_ns for f in (dest/'jobs').glob('*/complete.txt')};saved=(dest/'frozen_matrix.dat').read_bytes()
                with (out/(dest.name+'_resume.log')).open('w')as log:subprocess.run(cmd,env=dict(env,TASK_MBS_WARM_CACHE='1'),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
                assert saved==(dest/'frozen_matrix.dat').read_bytes()
                assert markers=={str(f):f.stat().st_mtime_ns for f in (dest/'jobs').glob('*/complete.txt')}
                if mode=='coherent':
                    damaged=next((dest/'jobs').glob('j*/result/result.dat'));original=damaged.read_bytes()
                    try:
                        damaged.write_bytes(original+b'changed\n')
                        rejected=subprocess.run(cmd,env=dict(env,TASK_MBS_WARM_CACHE='1'),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=600)
                        assert rejected.returncode and 'cached job changed' in rejected.stdout
                    finally:damaged.write_bytes(original)
        error=float(np.max(abs(pair[0][:,2:]-pair[1][:,2:])/pair[0][:,2,None]));assert error<1e-12,(mode,error)
        reports.append(dict(mode=mode,all16_difference_scaled_M11=error,seconds=times,resume_same_matrix_and_no_workers=True))
        (out/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
