#!/usr/bin/env python3
"""Compare actual incoherent upper/lower/prefix products on identical poses."""
import argparse,json,os,subprocess,time
from pathlib import Path
from test_droxtal_resident import load,difference,integral_difference

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--backend',choices=['cpu','cuda'],default='cpu')
p.add_argument('--tolerance',type=float,default=1e-7)
a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=True)
binary=str(a.binary.resolve());repo=Path(__file__).resolve().parents[1]
grid=root/'theta.csv';grid.write_text('160\n163.75909223159732\n167.13352685050799\n170\n179.5\n')
common=['--method','po','--incoherent','--backend',a.backend,'--ri','1.3116','0',
        '--wavelength-um','.532','--threads','2','--haar-alpha','--theta-grid-file',str(grid),
        '--phi-points','8','--sobol-seed','128','23','--sobol-offset','256',
        '--orientation-chunk','64','--profile-phases','--close']
if a.backend=='cuda':common+=['--orientation-pipeline']
env=os.environ.copy();env['OMP_WAIT_POLICY']='PASSIVE';env['MBS_HOST_MEM_BUDGET_MB']='2048'
records=[]
def run(name,shape,depth,cutoff,extra=()):
    dest=root/name;dest.mkdir();args=common+shape+['--max-reflections',str(depth),'--cutoff-profile',cutoff]+list(extra)+['--output',str(dest/'result')]
    (dest/'command.args').write_text('\n'.join([binary]+args)+'\n');start=time.monotonic()
    with (dest/'stdout.log').open('w')as log:
        child=subprocess.run([binary]+args,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=240)
    assert child.returncode==0,(name,child.returncode,dest/'stdout.log')
    return dest/'result',time.monotonic()-start
cases=[('column501',['--particle','1','501.2','155.8','--symmetry','2','6']),
       ('column631',['--particle','1','631','174.8','--symmetry','2','6']),
       ('droxtal100',['--particle','4','50','--symmetry','2','6']),
       ('asymmetric',['--particle-file',str(repo/'tests/data/resident_asymmetric.particle'),
                      '--resize-dmax-um','25','--symmetry','1','1'])]
for name,shape in cases:
    for cutoff in ['off','balanced']:
        label=name+'_'+cutoff
        upper,tu=run(label+'_upper',shape,12,cutoff,['--orientation-prefix-output','64'])
        lower,tl=run(label+'_lower',shape,10,cutoff,['--orientation-prefix-output','64'])
        pair,tp=run(label+'_paired',shape,12,cutoff,['--paired-depth','10','--orientation-prefix-output','64'])
        errors={}
        for prefix in ['', '_prefix64']:
            for tag,reference,suffix in [('upper',upper,''),('lower',lower,'_depth10')]:
                stem='result'+prefix
                errors[tag+prefix]=difference(load(reference/(stem+'.dat')),load(pair/(stem+suffix+'.dat')))
                errors[tag+prefix+'_integrals']=integral_difference(reference/(stem+'_integrals.tsv'),pair/(stem+suffix+'_integrals.tsv'))
        record=dict(case=label,errors=errors,tolerance=a.tolerance,passed=max(errors.values())<=a.tolerance,
                    separate_seconds=tu+tl,paired_seconds=tp,speedup=(tu+tl)/tp)
        records.append(record);(root/'comparison.json').write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps(record),flush=True)
        assert record['passed'],record
