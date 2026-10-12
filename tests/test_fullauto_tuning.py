#!/usr/bin/env python3
"""Native tuning vs unchanged full-SO(3) reference, with independent matrix audits."""
import argparse
import json
import subprocess
import time
from pathlib import Path
import numpy as np
from audit_fullauto import run as audit
from test_fullauto_arbitrary import asymmetric_tetrahedron, rotated_box

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--backend',default='cpu')
    p.add_argument('--case',default='all')
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    asymmetric_tetrahedron(out/'asymmetric.particle');rotated_box(out/'absorbing.particle')
    cases=[('prism',['--particle','1','.2','.1'],'0','4'),
           ('asymmetric',['--particle-file',str(out/'asymmetric.particle')],'0','3'),
           ('absorbing',['--particle-file',str(out/'absorbing.particle')],'.02','3'),
           ('concave',['--particle-file',str(ROOT/'examples/particles/concave_hexagonal.particle'),'--resize-dmax-um','.3'],'0','4')]
    reports=[]
    for name,shape,imag,depth in cases:
        if a.case not in ('all',name):continue
        for mode in ('auto','off'):
            dest=out/(name+'_'+mode)
            command=[str(a.binary.resolve()),'--method','po','--backend',a.backend]+shape+[
                '--ri','1.3116' if imag=='0' else '1.45',imag,'--wavelength-um','1','--max-reflections',depth,
                '--fullauto','.01','--fullauto-resources',mode,'--fullauto-rotation-symmetry',mode,
                '--fullauto-cutoffs','auto' if mode=='auto' else 'off','--fullauto-mirror','off',
                '--phi-points','4','--scattering-grid','0','25','4','32','--fullauto-pilot','64',
                '--fullauto-kernel-count','256','--fullauto-initial','256','--fullauto-min-correction','64',
                '--fullauto-theta-start','5','--max-orientations','131072','--output',str(dest)]
            if mode=='off':command+=['--cutoff-profile','off','--threads','2']
            if a.backend=='cuda':command+=['--fullauto-gpus','auto']
            (out/(name+'_'+mode+'.args')).write_text('\n'.join(command)+'\n')
            start=time.monotonic()
            with (out/(name+'_'+mode+'.log')).open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
            elapsed=time.monotonic()-start
            report=audit(dest,ROOT/'scripts/gpu_campaign.py',dense=True)
            before={str(f):f.stat().st_mtime_ns for f in (dest/'jobs').glob('*/complete.txt')}
            with (out/(name+'_'+mode+'_resume.log')).open('w') as log:
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
            assert before=={str(f):f.stat().st_mtime_ns for f in (dest/'jobs').glob('*/complete.txt')}
            state=json.loads((dest/'fullauto_status.json').read_text())
            if name=='asymmetric':assert state['rotation_beta_factor']==state['rotation_gamma_factor']==1
            if mode=='off':assert state['worker_threads']==2 and state['worker_chunk']==0
            report.update(case=name,mode=mode,seconds=elapsed,resources=json.loads((dest/'resource_calibration.json').read_text()),
                          rotations=json.loads((dest/'rotation_verification.json').read_text()),cutoffs=json.loads((dest/'cutoff_calibration.json').read_text()),resume_preserved_jobs=True)
            reports.append(report);(out/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n')
        pair=[r for r in reports if r['case']==name]
        ref=np.loadtxt(out/(name+'_off')/'mueller_fullauto.dat',skiprows=1)
        candidate=np.loadtxt(out/(name+'_auto')/'mueller_fullauto.dat',skiprows=1)
        error=float(np.max(abs(candidate[:,2:]-ref[:,2:])/ref[:,2,None]))
        ratios=float(np.max(abs(candidate[:,3:]/candidate[:,2,None]-ref[:,3:]/ref[:,2,None])))
        bound=sum(r['statistical95_scaled_M11']+r['guard_residual95_scaled_M11'] for r in pair)
        assert error<=bound,(name,error,bound)
        for report in pair:report.update(all16_difference_scaled_M11=error,normalized_difference=ratios,
                                         auto_vs_reference_speedup=pair[1]['seconds']/pair[0]['seconds'])
        (out/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n')
    print(json.dumps(reports,indent=2))

if __name__=='__main__':main()
