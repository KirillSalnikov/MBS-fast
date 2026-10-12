#!/usr/bin/env python3
"""Check exact Sobol windows, checkpoint block reuse and full-matrix audits."""
import argparse
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from audit_fullauto import run as audit

ROOT=Path(__file__).resolve().parents[1]


def execute(command, log):
    log.write_text('\n'.join(command)+'\n')
    with log.with_suffix('.stdout').open('w') as stream:
        subprocess.run(command,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=1800)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=ROOT/'bin/mbs_po')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    base=[str(a.binary.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1',
          '--ri','1.3116','0','--wavelength-um','1','--max-reflections','2','--threads','2']
    whole=[]
    for controls in [False,True]:
        control=['--analytic-facet-average','--analytic-shadow-control','circular'] if controls else []
        arrays=[]
        for index,(count,offset) in enumerate([(256,0),(64,0),(64,64),(64,128),(64,192)]):
            dest=out/('window_'+str(controls)+'_'+str(index))
            command=base+['--haar-alpha','--sobol-seed',str(count),'11','--sobol-offset',str(offset),
                          '--scattering-grid','0','180','4','16','--output',str(dest),'--close']+control
            execute(command,out/(dest.name+'.args'))
            arrays.append(np.loadtxt(dest/(dest.name+'.dat'),skiprows=1))
        combined=np.stack(arrays[1:]).mean(axis=0)
        difference=float(np.max(abs(combined[:,2:]-arrays[0][:,2:])/arrays[0][:,2,None]))
        assert difference<1e-11,(controls,difference)
        whole.append(dict(controls=controls,all16_difference_scaled_M11=difference))
    records=[]
    for policy in ['off','on']:
        dest=out/('fullauto_'+policy)
        command=base+['--fullauto','.03','--fullauto-mirror','off','--phi-points','4',
                      '--scattering-grid','0','25','4','16','--fullauto-pilot','32',
                      '--fullauto-initial','256','--fullauto-min-correction','32',
                      '--fullauto-theta-start','5','--fullauto-max-rounds','16',
                      '--fullauto-orientation-block','64','--fullauto-incremental',policy,
                      '--max-orientations','16384','--output',str(dest)]
        execute(command,out/(dest.name+'.args'))
        report=audit(dest,ROOT/'scripts/gpu_campaign.py',dense=True)
        state=json.loads((dest/'fullauto_status.json').read_text())
        records.append(dict(policy=policy,state=state,audit=report))
        if policy=='on':
            marker={str(p):p.stat().st_mtime_ns for p in (dest/'jobs').glob('*/complete.txt')}
            execute(command,out/'resume.args')
            assert marker=={str(p):p.stat().st_mtime_ns for p in (dest/'jobs').glob('*/complete.txt')}
    matrices=[np.loadtxt(out/('fullauto_'+p)/'mueller_fullauto.dat',skiprows=1) for p in ['off','on']]
    difference=float(np.max(abs(matrices[0][:,2:]-matrices[1][:,2:])/matrices[0][:,2,None]))
    same_counts=records[0]['state']['counts_per_seed']==records[1]['state']['counts_per_seed']
    if same_counts:
        assert difference<1e-10,difference
    else:
        # Block startup changes the measured level costs and therefore the
        # allocation. Both estimators must independently satisfy the budget.
        bound=sum(r['audit']['statistical95_scaled_M11']+
                  r['audit']['independent_dense']['residual_plus95_scaled_M11'] for r in records)
        assert difference<=bound,(difference,bound)
    result=dict(windows=whole,controllers=records,all16_controller_difference_scaled_M11=difference,
                completed_blocks_untouched_on_resume=True)
    result['same_final_counts']=same_counts
    (out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
