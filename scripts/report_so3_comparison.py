#!/usr/bin/env python3
"""Compare complete randomized SO(3) ensembles; high-N Sobol is a reference,
not an exact solution. All errors are scaled by reference M11 per angle.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.stats import t


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--backend',choices=['cpu','gpu'],default='gpu')
    a=p.parse_args();groups={}
    for marker in (a.root/'results'/a.backend).glob('*/complete.json'):
        r=json.loads(marker.read_text());j=r['job']
        data=np.loadtxt(marker.parent/'matrix/matrix.dat',skiprows=1)
        groups.setdefault((j['case'],j['method'],j['level']),[]).append((j['replica'],r,data))
    summaries=[]
    for key,items in sorted(groups.items()):
        if len(items)!=8:continue
        items=sorted(items,key=lambda x:x[0]);stack=np.stack([x[2][:,2:] for x in items])
        mean=stack.mean(0);variance=stack.var(0,ddof=1)/len(items)
        scale=np.maximum(abs(mean[:,0]),1e-300)
        reference=groups.get((key[0],'sobol','reference'),[])
        row=dict(case=key[0],method=key[1],level=key[2],orientations_per_replica=items[0][1]['job']['orientations'],
                 replicas=len(items),total_seconds=sum(x[1]['seconds'] for x in items),
                 mean_seconds=float(np.mean([x[1]['seconds'] for x in items])),
                 statistical95_scaled_M11=float(np.max(t.ppf(.975,7)*np.sqrt(variance)/scale[:,None])))
        if len(reference)==8:
            refs=np.stack([x[2][:,2:] for x in reference]);ref=refs.mean(0);refvar=refs.var(0,ddof=1)/8
            scale=np.maximum(abs(ref[:,0]),1e-300)
            row['statistical95_scaled_M11']=float(np.max(t.ppf(.975,7)*np.sqrt(variance)/scale[:,None]))
            row['difference_to_reference_scaled_M11']=float(np.max(abs(mean-ref)/scale[:,None]))
            row['reference_statistical95_scaled_M11']=float(np.max(t.ppf(.975,7)*np.sqrt(refvar)/scale[:,None]))
            row['difference_by_element']={label:float(np.max(abs(mean[:,col]-ref[:,col])/scale))
                                          for label,col in [('M11',0),('M12',1),('M22',5)]}
            row['normalized_ratio_difference']={label:float(np.max(abs(mean[:,col]/mean[:,0]-ref[:,col]/ref[:,0])))
                                                for label,col in [('M12/M11',1),('M22/M11',5)]}
        summaries.append(row)
    report=dict(completed_ensembles=summaries,
                precision_contract='Pointwise Student-t 95% dispersion of eight independent scrambles/global Haar rotations; no simultaneous guarantee. High-N Sobol has finite uncertainty; no ADDA accuracy claim.',
                error_contract='Max across selected theta and all 16 Mueller elements, divided by reference M11 at each theta when available; otherwise by current mean M11. Selected nodes only, not a continuous-angle guarantee.',
                timings='Sum of measured process wall times, including setup and output; compare on the same host/backend.')
    (a.root/('comparison_'+a.backend+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
