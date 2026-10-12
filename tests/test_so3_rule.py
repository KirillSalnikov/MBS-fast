#!/usr/bin/env python3
"""Verify weighted Euler tracing against identical native Sobol samples."""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=ROOT/'bin/mbs_po')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    exporter=out/'export_sobol_rule'
    subprocess.run(['g++','-O2','-std=c++11','-I'+str(ROOT/'src'),
                    str(ROOT/'tests/export_sobol_rule.cpp'),str(ROOT/'src/math/Sobol.cpp'),
                    '-o',str(exporter)],check=True)
    rule=out/'sobol.rule'
    with rule.open('w') as stream:
        subprocess.run([str(exporter),'256','11'],stdout=stream,check=True)
    base=[str(a.binary.resolve()),'--method','po','--backend','cpu',
          '--particle','1','.2','.1','--symmetry','1','1','--ri','1.3116','0',
          '--wavelength-um','1','--max-reflections','2','--threads','2',
          '--scattering-grid','0','180','4','16']
    reports=[]
    for control in [False,True]:
        arrays=[]
        for name,options in [('sobol',['--haar-alpha','--sobol-seed','256','11']),
                             ('rule',['--so3-rule',str(rule)])]:
            dest=out/(name+'_'+str(control))
            args=base+options+['--output',str(dest),'--close']
            if control:args+=['--analytic-facet-average','--analytic-shadow-control','circular']
            (dest.with_suffix('.args')).write_text('\n'.join(args)+'\n')
            with dest.with_suffix('.log').open('w') as log:
                subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
            matrix=dest/(dest.name+'.dat')
            assert matrix.exists(),'Orientation mode returned success without a matrix'
            arrays.append(np.loadtxt(matrix,skiprows=1))
        assert np.array_equal(arrays[0][:,:2],arrays[1][:,:2])
        difference=float(np.max(abs(arrays[0][:,2:]-arrays[1][:,2:])/arrays[0][:,2,None]))
        assert difference<1e-11,(control,difference)
        reports.append(dict(analytic_controls=control,all16_difference_scaled_M11=difference))
    invalid=['10 20 30 -1\n','10 20 30 .5\n','190 20 30 1\n','10 20 30 1 extra\n']
    for i,text in enumerate(invalid):
        path=out/('invalid_'+str(i)+'.rule');path.write_text(text)
        r=subprocess.run(base+['--so3-rule',str(path),'--output',str(out/'invalid')],
                         stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        assert r.returncode!=0,r.stdout
    (out/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))

if __name__=='__main__':main()
