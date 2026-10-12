#!/usr/bin/env python3
"""All-Mueller and dense-grid gates for local budgets on varied particles."""
import argparse,json,subprocess,time
from pathlib import Path
import numpy as np
from audit_fullauto import run as audit
from test_fullauto_arbitrary import asymmetric_tetrahedron,rotated_box

ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--backend',default='cpu');p.add_argument('--case',default='all');a=p.parse_args()
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=True);tet=out/'tetrahedron.particle';box=out/'box.particle';asymmetric_tetrahedron(tet);rotated_box(box)
    cases=[('prism',['--particle','1','.2','.1'],'0','4'),('asymmetric',['--particle-file',str(tet)],'0','3'),
           ('absorbing',['--particle-file',str(box)],'.02','3'),('concave',['--particle-file',str(ROOT/'examples/particles/concave_hexagonal.particle'),'--resize-dmax-um','.3'],'0','4')]
    reports=[]
    for case,shape,imag,depth in cases:
        if a.case!='all' and a.case!=case:continue
        for policy in ['auto','off']:
            dest=out/(case+'_'+policy);command=[str(a.binary.resolve()),'--method','po','--backend',a.backend]+shape+[
                '--ri','1.3116' if imag=='0' else '1.45',imag,'--wavelength-um','1','--max-reflections',depth,'--threads','2',
                '--fullauto','.01','--fullauto-local-orientations',policy,'--fullauto-mirror','off','--phi-points','4','--scattering-grid','0','25','4','32',
                '--fullauto-pilot','64','--fullauto-initial','256','--fullauto-min-correction','64','--fullauto-theta-start','5','--fullauto-max-rounds','32',
                '--fullauto-orientation-block','256','--max-orientations','131072','--cutoff-profile','off','--output',str(dest)]
            (out/(case+'_'+policy+'.args')).write_text('\n'.join(command)+'\n');start=time.monotonic()
            with (out/(case+'_'+policy+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
            seconds=time.monotonic()-start;report=audit(dest,ROOT/'scripts/gpu_campaign.py',dense=True)
            records={str(f):f.stat().st_mtime_ns for f in (dest/'jobs').glob('*/complete.txt')}
            with (out/(case+'_'+policy+'_resume.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
            assert records=={str(f):f.stat().st_mtime_ns for f in (dest/'jobs').glob('*/complete.txt')},'Resume recomputed finished work'
            budget=np.loadtxt(dest/'orientation_budget.csv',delimiter=',',skiprows=1)[:,1:]
            report.update(case=case,local=policy,seconds=seconds,orientation_row_work_per_seed=float(budget.sum()),
                          global_peak_row_work_per_seed=float(budget.max(0).sum()*len(budget)),completed_jobs_untouched_on_resume=True)
            reports.append(report);(out/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n')
    for case in {r['case'] for r in reports}:
        pair=[r for r in reports if r['case']==case];assert len(pair)==2
        matrices=[np.loadtxt(out/(case+'_'+r['local'])/'mueller_fullauto.dat',skiprows=1) for r in pair]
        error=float(np.max(abs(matrices[0][:,2:]-matrices[1][:,2:])/matrices[0][:,2,None]))
        bound=sum(r['statistical95_scaled_M11']+r['guard_residual95_scaled_M11'] for r in pair)
        assert error<=bound,(case,error,bound)
        for r in pair:r['difference_to_other_policy_scaled_M11']=error
    (out/'comparison.json').write_text(json.dumps(reports,indent=2)+'\n');print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
