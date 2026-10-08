#!/usr/bin/env python3
"""Native full-field checks for frozen analytic coefficients and CLI errors."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np

ROOT=Path(__file__).resolve().parents[1]


with tempfile.TemporaryDirectory(prefix='mbs_control_weights_') as work:
    work=Path(work);theta=np.array([0.,17.,45.,90.,170.,180.])
    np.savetxt(work/'theta.txt',theta,fmt='%.17g')
    command=[str(ROOT/'cpu/bin/mbs_po_mpi'),'--method','po','--backend','cpu',
             '--particle','1','4','4','--refractive-index','1.3116','0',
             '--wavelength-um','.532','--max-reflections','12','--cutoff-profile','off',
             '--sobol-seed','128','19','--haar-alpha','--analytic-facet-average',
             '--analytic-shadow-control','circular','--analytic-mean-cache',str(work/'means.cache'),
             '--theta-grid-file',str(work/'theta.txt'),'--phi-points','4','--threads','4','--close']
    def run(name,weights=None,extra=()):
        output=work/name
        cmd=command+['--output',str(output)]+list(extra)
        if weights:cmd+=['--analytic-control-weights',str(weights)]
        p=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        return p,output
    def file(name,coefficients):
        path=work/(name+'.tsv')
        np.savetxt(path,np.column_stack([theta,coefficients]),fmt='%.17g',comments='',
                   header='theta_deg reflection_weight shadow_weight')
        return path
    baseline,path=run('baseline',extra=['--profile-phases'])
    assert baseline.returncode==0,baseline.stderr
    base=np.loadtxt(path/'baseline.dat',skiprows=1)
    assert (path/'baseline_phase_timings.tsv').exists()
    for name,coef in [('ones',np.ones((6,2))),('zero',np.zeros((6,2))),
                      ('mixed',np.tile([2.,-.5],(6,1)))]:
        p,out=run(name,file(name,coef));assert p.returncode==0,p.stderr
        data=np.loadtxt(out/(name+'.dat'),skiprows=1)
        summary=np.atleast_1d(np.genfromtxt(out/(name+'_analytic_facets.tsv'),names=True,delimiter='\t'))
        expected=summary['raw_M11']+coef[:,0]*(summary['reflection_mean']-summary['sampled_reflection'])+coef[:,1]*(summary['shadow_mean']-summary['sampled_shadow'])
        np.testing.assert_allclose(data[:,2],expected,rtol=2e-13,atol=1e-12)
        error=np.max(abs(data[:,3:]-base[:,3:])/np.maximum(abs(base[:,2,None]),1e-100))
        assert error<1e-12,error
        if name=='ones':assert np.max(abs(data[:,2]-base[:,2])/abs(base[:,2]))<1e-12
    invalid=file('invalid',np.full((6,2),5.))
    p,_=run('bad',invalid);assert p.returncode!=0 and 'weight' in p.stderr
    invalid.write_text('theta_deg reflection_weight shadow_weight\n0 1 1\n')
    p,_=run('short',invalid);assert p.returncode!=0 and 'weight' in p.stderr
    p,_=run('cpu_pipeline',extra=['--orientation-pipeline'])
    assert p.returncode!=0 and 'CUDA' in p.stderr
print('PASS: independent frozen coefficient reconstruction; all15 other Mueller entries retained; malformed weights rejected; profiling emitted')
