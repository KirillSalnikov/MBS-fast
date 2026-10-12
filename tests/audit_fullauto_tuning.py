#!/usr/bin/env python3
"""Independently reconstruct the fresh-seed automatic-cutoff bound."""
import argparse,json,math
from pathlib import Path
import numpy as np

T=2.3646242510102993
def option(args,name,count=1):
    if name not in args:return None
    i=args.index(name)+1;return args[i:i+count]

def run(root):
    state=json.loads((root/'fullauto_status.json').read_text())
    budget=sum(state[k] for k in ['statistical_budget','interpolation_budget','mirror_budget','analytic_mean_budget','rotation_budget','cutoff_budget','execution_budget'])
    assert abs(budget-json.loads((root/'precision_contract.json').read_text())['interpolation_budget']*4)<1e-12
    saved=json.loads((root/'automatic_cutoff_validation.json').read_text())
    if 'profile' not in saved:return dict(passed=True,cutoff_audit_enabled=False)
    first=1600000000+state['round']*1024;seeds=[first+17*i for i in range(8)]
    samples={'off':{},saved['profile']:{}}
    for folder in (root/'jobs').iterdir():
        if not (folder/'complete.txt').exists():continue
        args=(folder/'command.args').read_text().splitlines()
        sobol=option(args,'--sobol-seed',2);depth=option(args,'--max-reflections');profile=option(args,'--cutoff-profile')
        if not sobol or not depth or not profile:continue
        count,seed=map(int,sobol)
        if seed not in seeds or count!=saved['count_per_seed'] or int(depth[0])!=state['target_depth'] or profile[0] not in samples:continue
        assert '--analytic-facet-average' not in args
        matrix=np.loadtxt(folder/'result/result.dat',skiprows=1,ndmin=2)
        assert matrix.shape[0]==saved['theta_count'];samples[profile[0]][seed]=matrix[:,2:]
    a,b=[np.array([samples[key][seed] for seed in seeds]) for key in [saved['profile'],'off']]
    am,bm=a.mean(0),b.mean(0);delta=a-b
    raw=(abs(delta.mean(0))+T*delta.std(0,ddof=1)/math.sqrt(8))/bm[:,0,None]
    ra,rb=am/am[:,0,None],bm/bm[:,0,None]
    influence=(a-ra[None]*a[:,:,0,None])/am[None,:,0,None]-(b-rb[None]*b[:,:,0,None])/bm[None,:,0,None]
    ratios=abs(ra-rb)+T*influence.std(0,ddof=1)/math.sqrt(8)
    bound=float(max(raw.max(),ratios[:,1:].max()))
    assert abs(bound-saved['bias95'])<1e-12,(bound,saved)
    assert (bound<=saved['budget'])==saved['passed']
    result=dict(passed=True,cutoff_audit_enabled=True,reconstructed_bias95=bound,budget=saved['budget'],count_per_seed=saved['count_per_seed'],theta_count=saved['theta_count'],all16_and_ratios_checked=True)
    (root/'independent_cutoff_audit.json').write_text(json.dumps(result,indent=2)+'\n');return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();print(json.dumps(run(a.root),indent=2))
