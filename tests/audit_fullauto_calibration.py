#!/usr/bin/env python3
"""Recompute frozen all-element control fits and phi scores with NumPy."""
import argparse,json
from pathlib import Path
import numpy as np
from audit_fullauto import option
from audit_fullauto_v2 import intervals,matrix
TRAIN=[1009,1013,1019,1021,1031,1033,1039,1049]
VALID=[2003,2011,2017,2027,2029,2039,2053,2063]
def main():
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();root=a.root.resolve();cal=json.loads((root/'calibration.json').read_text())
    weights=matrix(root/'polarization_weights.tsv');angles=weights[:,0];jobs=[]
    for folder in (root/'jobs').iterdir():
        if not (folder/'complete.txt').exists():continue
        args=(folder/'command.args').read_text().splitlines()
        if option(args,'--analytic-control-weights') is not None or '--haar-mirror-audit' in args:continue
        data=matrix(folder/'result/result.dat')
        if len(data)!=len(angles) or not np.allclose(data[:,0],angles,rtol=0,atol=1e-10):continue
        if int(option(args,'--max-reflections'))!=cal['calibration_depth']:continue
        n,seed=map(int,option(args,'--sobol-seed',2))
        if n!=cal['pilot_count'] or seed not in TRAIN+VALID:continue
        jobs.append((int(option(args,'--phi-points')),seed,data[:,2:],matrix(folder/'result/result_analytic_facets.tsv'),float((folder/'complete.txt').read_text().split()[1])))
    reports=[]
    for candidate in cal['phi_candidates']:
        phi=candidate['phi'];train=[next(j for j in jobs if j[0]==phi and j[1]==seed) for seed in TRAIN];valid=[next(j for j in jobs if j[0]==phi and j[1]==seed) for seed in VALID]
        coefficients=np.zeros((len(angles),16,2));coefficients[:,0]=1;adjusted=np.stack([j[2] for j in valid]);accepted=np.zeros(16,int)
        for r in range(len(angles)):
            controls=np.array([j[3][r,4:6] for j in train]);centred=controls-controls.mean(0);scale=centred.std(0,ddof=1)
            for element in range(16 if cal['polarization_controls'] else 1):
                y=np.array([j[3][r,1] if element==0 else j[2][r,element] for j in train]);active=scale>max(abs(y.mean()),1e-100)*1e-12
                z=np.zeros_like(controls);z[:,active]=centred[:,active]/scale[active]
                cy=z.T@(y-y.mean())/7;b=coefficients[r,element].copy()
                if active.all():b=np.linalg.solve(z.T@z/7+.1*np.eye(2),cy)/scale
                elif active.any():b[active]=cy[active]/1.1/scale[active]
                b=np.clip(b,-4,4);v=np.array([j[3][r] for j in valid]);raw=v[:,1] if element==0 else np.array([j[2][r,element] for j in valid]);delta=v[:,2:4]-v[:,4:6]
                default=raw+(delta.sum(1) if element==0 else 0);proposed=raw+delta@b
                if proposed.var(ddof=1)<.7*default.var(ddof=1):coefficients[r,element]=b;adjusted[:,r,element]=proposed;accepted[element]+=1
                else:adjusted[:,r,element]=default
        _,ci,ratio=intervals(adjusted);score=sum(j[4] for j in valid)*max(float(ci.max()),float(ratio.max()) if cal['score_includes_ratios'] else 0)**2
        assert np.isclose(score,candidate['score'],rtol=1e-10,atol=1e-16),(phi,score,candidate['score'])
        assert accepted.tolist()==candidate['accepted_control_rows']
        if phi==cal['selected_phi']:
            error=float(np.max(abs(coefficients.reshape(len(angles),32)-weights[:,1:])))
            assert error<1e-10,error
        reports.append(dict(phi=phi,python_score=score,native_score=candidate['score'],accepted_control_rows=accepted.tolist()))
    assert min(reports,key=lambda x:x['python_score'])['phi']==cal['selected_phi']
    result=dict(passed=True,working_grid_recheck=cal['working_grid_recheck'],depth=cal['calibration_depth'],selected_phi=cal['selected_phi'],all16_control_coefficients_checked=True,candidates=reports)
    (root/'independent_calibration_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
