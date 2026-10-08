#!/usr/bin/env python3
"""Calibrate azimuth/controls on disjoint pilots and validate kernel tuning."""
import argparse
import fcntl
import json
from pathlib import Path

import numpy as np

from gpu_campaign import (JobRunner, TRAIN_SEEDS, VALIDATION_SEEDS, atomic_json,
                          confidence, estimate_cost_score, fit_control_weights,
                          reweight, write_weights)


def calibrate(runner, pilot_count=32768, phi_candidates=(4,8,16), kernel_count=8192, fine_depth=12):
    result=dict(status='active',physical=runner.physical,binary_sha256=runner.binary_hash,
                training_seeds=TRAIN_SEEDS,validation_seeds=VALIDATION_SEEDS,
                pilot_count=pilot_count,fine_depth=fine_depth,kernel_candidates=[],phi_candidates=[])
    output=runner.work/'calibration.json'
    # Tune the costly fine kernel against identical orientations and all rows.
    reference=None;kernel_results=[]
    for pipeline,block,warp in [(False,64,'auto'),(True,64,'auto'),
                                (True,128,'auto'),(True,256,'auto'),(True,64,'thread')]:
        job=runner.job(runner.gpus[-1],kernel_count,4001,fine_depth,'off',4,
                       pipeline=pipeline,block=block,warp=warp,profile=True)
        data=np.loadtxt(job['data'],skiprows=1)
        if reference is None:reference=data
        error=float(np.max(abs(data[:,2:]-reference[:,2:])/np.maximum(abs(reference[:,2,None]),1e-100)))
        row=dict(pipeline=pipeline,block=block,warp=warp,seconds=job['seconds'],
                 max_mueller_difference_scaled_M11=error,passed=error<1e-7,
                 job=job,phase_timings=job['phase_timings'])
        kernel_results.append(row);result['kernel_candidates']=kernel_results
        atomic_json(output,result);print('kernel',pipeline,block,warp,job['seconds'],error,flush=True)
    valid=[r for r in kernel_results if r['passed']]
    if not valid:raise RuntimeError('no GPU kernel candidate passed field validation')
    winner=min(valid,key=lambda r:r['seconds'])
    options=dict(pipeline=winner['pipeline'],block=winner['block'],warp=winner['warp'])
    result['selected_kernel']=options
    for phi in phi_candidates:
        train,summaries_train,train_time,_=runner.stage(pilot_count,phi=phi,seeds=TRAIN_SEEDS,**options)
        valid,summaries_valid,valid_time,_=runner.stage(pilot_count,phi=phi,seeds=VALIDATION_SEEDS,**options)
        coefficients,fit=fit_control_weights(train,valid,summaries_train,summaries_valid)
        adjusted=reweight(valid,summaries_valid,coefficients)
        baseline_metric,_,_=confidence(valid)
        fitted_metric,_,_=confidence(adjusted)
        weights=runner.work/f'control_weights_phi{phi}.tsv'
        write_weights(weights,valid[0,:,0],coefficients)
        record=dict(phi=phi,stage_seconds=valid_time,training_seconds=train_time,
                    baseline_precision=baseline_metric,weighted_precision=fitted_metric,
                    time_variance_score=estimate_cost_score(adjusted,valid_time),
                    weights=str(weights),fit=fit)
        result['phi_candidates'].append(record);atomic_json(output,result)
        print('phi',phi,'seconds',valid_time,'weighted95',fitted_metric['max_mueller95_scaled_M11'],
              'fitted rows',fit['rows_fitted'],flush=True)
    selected=min(result['phi_candidates'],key=lambda r:r['time_variance_score'])
    result.update(status='calibrated',selected_phi=selected['phi'],
                  weights=selected['weights'] if runner.controls=='facets' else None,
                  scope='Pilot time*worst normalized95interval² ranking; production error not established by pilot.',
                  physical=runner.physical)
    atomic_json(output,result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--case',choices=['test1','test2','test3'],default='test1')
    p.add_argument('--height-um',type=float);p.add_argument('--diameter-um',type=float)
    p.add_argument('--refractive-index',type=float,nargs=2,default=[1.3116,0.])
    p.add_argument('--wavelength-um',type=float,default=.532)
    p.add_argument('--fine-depth',type=int,default=12)
    p.add_argument('--analytic-controls',choices=['facets','off'],default='facets')
    p.add_argument('--theta-grid-file',type=Path);p.add_argument('--mean-cache',type=Path)
    p.add_argument('--gpus',type=int,nargs='+',default=[0,1,2,3])
    p.add_argument('--threads',type=int,default=16);p.add_argument('--orientation-chunk',type=int,default=256)
    p.add_argument('--pilot-count',type=int,default=32768);p.add_argument('--kernel-count',type=int,default=8192)
    p.add_argument('--phi-candidates',type=int,nargs='+',default=[4,8,16]);p.add_argument('--cuda-lib-dir',type=Path)
    args=p.parse_args()
    if min(args.pilot_count,args.kernel_count,*args.phi_candidates)<=0:p.error('counts must be positive')
    runner=JobRunner(args.binary,args.output,case=args.case,height=args.height_um,diameter=args.diameter_um,
                     index=args.refractive_index,wave=args.wavelength_um,theta_file=args.theta_grid_file,
                     gpus=args.gpus,threads=args.threads,chunk=args.orientation_chunk,
                     cuda_lib_dir=args.cuda_lib_dir,mean_cache=args.mean_cache,controls=args.analytic_controls)
    lock=(runner.work/'controller.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    result=calibrate(runner,args.pilot_count,args.phi_candidates,args.kernel_count,args.fine_depth)
    print(json.dumps({k:v for k,v in result.items() if k not in ['kernel_candidates','phi_candidates']},indent=2))


if __name__=='__main__':main()
