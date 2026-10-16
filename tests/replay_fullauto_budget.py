"""Replay real saved fullauto uncertainties through the native acceptance code.

Input: bounds.npz and production_bounds.npz produced by the saved-data diagnostic.
No simulated optical data, new orientations, or server jobs are used here.
"""
import argparse
import json
from pathlib import Path
import subprocess
import numpy as np


def fixture(path, d, indices=None):
    angles=d['angles'];requested=d['requested']
    ids=np.searchsorted(requested,angles)
    assert np.allclose(requested[ids],angles,rtol=0,atol=1e-10)
    ci=d['direct_stat'];ri=d['direct_ratios'];oci=d['output_stat'];ori=d['output_ratios']
    if 'interpolation_matrix' in d:
        raw=d['interpolation_matrix'];ratios=d['interpolation_ratios']
    else:
        # Round19 archived only the maximum interpolation bound per angle.
        # Reusing that maximum for every component is conservative.
        raw=np.repeat(d['interpolation'][:,None],16,axis=1);ratios=raw.copy();ratios[:,0]=0
    if indices is not None:
        requested=angles[indices];oci=oci[ids[indices]];ori=ori[ids[indices]]
        angles=angles[indices];ci=ci[indices];ri=ri[indices];raw=raw[indices];ratios=ratios[indices];ids=np.arange(len(indices))
    out_raw=np.full_like(oci,.005);out_ratio=np.full_like(ori,.005);out_ratio[:,0]=0
    out_raw[ids]=raw;out_ratio[ids]=ratios
    reserve=.00548828125
    direct=np.maximum(ci+raw,ri+ratios)+reserve
    output=np.maximum(oci+out_raw,ori+out_ratio)+reserve
    failures=int((np.maximum(direct.max(1),output[ids].max(1))>.02).sum())
    with path.open('w') as f:
        f.write(f'{len(requested)} {len(angles)}\n')
        np.savetxt(f,requested[None],fmt='%.17g');np.savetxt(f,ids[None],fmt='%d')
        for array in (ci,ri,oci,ori,raw,ratios):np.savetxt(f,array,fmt='%.17g')
        f.write(f'{direct.max():.17g} {output.max():.17g} {failures}\n')
    return dict(checked_angles=len(angles),output_angles=len(requested),direct_max=float(direct.max()),output_max=float(output.max()),failed_checked_angles=failures)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True);report={}
    for name,source in [('independent','bounds.npz'),('production','production_bounds.npz')]:
        d=np.load(args.data/source);path=args.output/(name+'.txt');report[name]=fixture(path,d)
        result=subprocess.run([str(args.probe.resolve()),str(path.resolve()),str((args.output/name).resolve())],check=True,stdout=subprocess.PIPE,text=True)
        print(result.stdout,end='')
    (args.output/'replay.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
