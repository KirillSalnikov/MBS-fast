#!/usr/bin/env python3
"""Compare saved-prefix replay before/after, with real optical data and restart."""
import argparse,json,os,subprocess
from pathlib import Path
import numpy as np
from test_fullauto_arbitrary import asymmetric_tetrahedron

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--before',type=Path,required=True);p.add_argument('--after',type=Path,required=True)
p.add_argument('--worker',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
tet=out/'asymmetric.particle';asymmetric_tetrahedron(tet)
env=dict(os.environ,TASK_MBS_WORKER=str(a.worker.resolve()));results=[]
for mode in ['coherent','incoherent','controls','asymmetric']:
    matrices=[];profiles={}
    for label,probe in [('before',a.before),('after',a.after)]:
        dest=out/(mode+'_'+label)
        shape=['--particle-file',str(tet),'--resize-dmax-um','.2'] if mode=='asymmetric' else ['--particle','1','.2','.1']
        command=[str(probe.resolve()),'--method','po','--backend','cpu']+shape+[
            '--ri','1.3116','0','--wavelength-um','1','--max-reflections','4','--threads','1',
            '--fullauto','.02','--fullauto-initial','1024','--fullauto-orientation-block','64',
            '--fullauto-analytic-controls','on' if mode=='controls' else 'off','--fullauto-mirror','off',
            '--fullauto-resources','off','--cutoff-profile','off','--phi-points','4',
            '--scattering-grid','0','180','4','64','--output',str(dest),'--close']
        if mode=='incoherent':command+=['--incoherent']
        (out/(dest.name+'.args')).write_text('\n'.join(command)+'\n')
        with (out/(dest.name+'.log')).open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
        profiles[label]=json.loads((dest/'replay_profile.json').read_text());matrices.append(np.loadtxt(dest/'frozen_matrix.dat',skiprows=1))
        if label=='after':
            assert profiles[label]['new_derived_jobs']==8,profiles[label]
            data=(dest/'frozen_matrix.dat').read_bytes();markers={str(f):f.stat().st_mtime_ns for f in (dest/'jobs').glob('*/complete.txt')}
            with (out/(dest.name+'_resume.log')).open('w')as log:subprocess.run(command,env=dict(env,TASK_MBS_WARM_CACHE='1'),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
            assert data==(dest/'frozen_matrix.dat').read_bytes()
            assert markers=={str(f):f.stat().st_mtime_ns for f in (dest/'jobs').glob('*/complete.txt')}
            if mode=='coherent':
                damaged=next((dest/'jobs').glob('j*/result/result.dat'));saved=damaged.read_bytes()
                try:
                    damaged.write_bytes(saved+b'changed\n')
                    failed=subprocess.run(command,env=dict(env,TASK_MBS_WARM_CACHE='1'),stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=600)
                    assert failed.returncode and 'cached job changed' in failed.stdout,failed.stdout[-2000:]
                finally:damaged.write_bytes(saved)
    difference=float(np.max(np.abs(matrices[0][:,2:]-matrices[1][:,2:])/matrices[0][:,2,None]))
    assert difference<1e-12,(mode,difference)
    result=dict(mode=mode,all16_difference_scaled_M11=difference,profiles=profiles,
                cache_replay_speedup=profiles['before']['seconds']/profiles['after']['seconds'],
                scope='Measured cached angular-subset assembly only; excludes native optics, startup and full campaign',
                resume_preserves_blocks=True,integral_sources_complete=True)
    results.append(result);(out/'comparison.json').write_text(json.dumps(results,indent=2)+'\n');print(mode,'PASS',json.dumps(result),flush=True)
