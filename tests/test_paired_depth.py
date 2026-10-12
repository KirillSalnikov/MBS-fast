#!/usr/bin/env python3
"""Check coherent shared-depth matrices and measured wall time."""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
from test_fullauto_arbitrary import asymmetric_tetrahedron,rotated_box
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--backend',default='cpu');p.add_argument('--count',type=int,default=2048);p.add_argument('--threads',type=int,default=4);p.add_argument('--case',default='all');a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    tet=out/'tetrahedron.particle';box=out/'box.particle';asymmetric_tetrahedron(tet);rotated_box(box)
    cases=[('prism',['--particle','1','10','6'],'1.3116','0'),('asymmetric',['--particle-file',str(tet)],'1.3116','0'),('absorbing',['--particle-file',str(box)],'1.45','.02')]
    reports=[]
    for case,shape,real,imag in cases:
        if a.case!='all' and case!=a.case:continue
        cache=out/(case+'_means.cache');matrices={};seconds={}
        base=[str(a.binary.resolve()),'--method','po','--backend',a.backend]+shape+['--ri',real,imag,'--wavelength-um','.532','--symmetry','1','1','--threads',str(a.threads),'--sobol-seed',str(a.count),'11','--haar-alpha','--scattering-grid','0','25','8','128','--analytic-facet-average','--analytic-shadow-control','circular','--analytic-mean-cache',str(cache),'--beam-cutoff','.001','--profile-phases','--close']
        for mode,extra in [('warm',['--max-reflections','12']),('lower',['--max-reflections','8']),('upper',['--max-reflections','12']),('joint',['--max-reflections','12','--paired-depth','8'])]:
            dest=out/(case+'_'+mode);command=base+extra+['--output',str(dest)];(out/(dest.name+'.args')).write_text('\n'.join(command)+'\n');start=time.monotonic()
            with (out/(dest.name+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)
            seconds[mode]=time.monotonic()-start;matrices[mode]=np.loadtxt(dest/(dest.name+'.dat'),skiprows=1)
        low=np.loadtxt(out/(case+'_joint')/(case+'_joint_depth8.dat'),skiprows=1)
        errors=[float(np.max(abs(low[:,2:]-matrices['lower'][:,2:])/matrices['lower'][:,2,None])),float(np.max(abs(matrices['joint'][:,2:]-matrices['upper'][:,2:])/matrices['upper'][:,2,None]))]
        assert max(errors)<(1e-6 if a.backend=='gpu' else 1e-10),(case,errors)
        reports.append(dict(case=case,backend=a.backend,all16_lower_difference_scaled_M11=errors[0],all16_upper_difference_scaled_M11=errors[1],wall_seconds=seconds,paired_speedup=(seconds['upper']+seconds['lower'])/seconds['joint']))
        (out/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
