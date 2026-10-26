#!/usr/bin/env python3
"""Exercise whole-orientation CUDA fallback with a generated 40-sided prism."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--binary',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    n=40;bottom=[(.2*math.cos(2*math.pi*i/n),.2*math.sin(2*math.pi*i/n),-.15) for i in range(n)]
    top=[(x,y,.15) for x,y,z in bottom]
    faces=[list(reversed(bottom)),top]
    faces += [[bottom[i],bottom[(i+1)%n],top[(i+1)%n],top[i]] for i in range(n)]
    particle=out/'prism40.particle'
    particle.write_text('0\n0\n180 360\n\n'+'\n\n'.join('\n'.join('%.17g %.17g %.17g'%v for v in f) for f in faces)+'\n')
    results=[]
    env=dict(os.environ,MBS_GPU_MULTI='0',MBS_GPU_MULTI_MAX='1',OMP_NUM_THREADS='2')
    for incoherent in [False,True]:
        matrices={};logs={}
        for backend in ['cpu','cuda']:
            name=('incoherent_' if incoherent else 'coherent_')+backend
            dest=out/name
            cmd=[str(a.binary.resolve()),'--method','po','--backend',backend,
                 '--particle-file',str(particle),'--ri','1.3116','0','--wavelength-um','.532',
                 '--max-reflections','4','--sobol-seed','8','42','--haar-alpha',
                 '--symmetry','1','1','--scattering-grid','0','180','4','24',
                 '--cutoff-profile','off','--threads','2','--output',str(dest),'--close']
            if incoherent:cmd.append('--incoherent')
            with (out/(name+'.log')).open('w') as log:
                subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
            logs[backend]=(out/(name+'.log')).read_text()
            data=(dest/(name+'.dat')).read_text().splitlines()[1:]
            matrices[backend]=[list(map(float,r.split())) for r in data if r.strip()]
        assert 'GPU aperture fallback:' in logs['cuda'],'Fallback not exercised'
        cpu,gpu=matrices['cpu'],matrices['cuda'];assert len(cpu)==len(gpu)==25
        error=0
        for c,g in zip(cpu,gpu):
            assert len(c)==len(g)==18 and c[:2]==g[:2]
            assert all(math.isfinite(v) for v in c+g) and c[2]>0
            error=max(error,max(abs(x-y)/c[2] for x,y in zip(c[2:],g[2:])))
        assert error<1e-6,error
        results.append({'incoherent':incoherent,'max_all16_scaled_M11':error,'fallback_exercised':True})
        print('PASS',results[-1],flush=True)
    (out/'summary.json').write_text(json.dumps(results,indent=2))

if __name__=='__main__':main()
