#!/usr/bin/env python3
"""Adaptive2/3+level full-Mueller estimator with independent calibration."""
import argparse
import fcntl
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np

from gpu_campaign import (JobRunner, PRODUCTION_SEEDS, TRAIN_SEEDS, VALIDATION_SEEDS,
                          allocation, atomic_json, confidence, fingerprint)


def difference(upper, lower):
    if not np.array_equal(upper[:, :, :2],lower[:, :, :2]):
        raise ValueError('paired levels have different angular coordinates')
    result=upper.copy();result[:,:,2:]-=lower[:,:,2:]
    return result


def combine(levels):
    result=levels[0].copy()
    for delta in levels[1:]:
        if not np.array_equal(delta[:,:,0],result[:,:,0]):
            raise ValueError('multilevel angular grids differ')
        result[:,:,2:]+=delta[:,:,2:]
    return result


def validate_level_order(depths, cutoffs):
    """Increasing depth or stricter cutoff at an unchanged final depth."""
    if len(depths)!=len(cutoffs) or min(depths)<=0 or depths!=sorted(depths):
        raise ValueError('depths must be positive and nondecreasing')
    for previous,current,left,right in zip(depths,depths[1:],cutoffs,cutoffs[1:]):
        if previous==current:
            before=0. if left=='off' else float(left)
            after=0. if right=='off' else float(right)
            if before<=after:
                raise ValueError('equal depths require a strictly smaller cutoff')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--case',choices=['test1','test2','test3'],default='test1')
    p.add_argument('--height-um',type=float);p.add_argument('--diameter-um',type=float)
    p.add_argument('--refractive-index',type=float,nargs=2,default=[1.3116,0.]);p.add_argument('--wavelength-um',type=float,default=.532)
    p.add_argument('--theta-grid-file',type=Path);p.add_argument('--mean-cache',type=Path)
    p.add_argument('--analytic-controls',choices=['facets','off'],default='facets')
    p.add_argument('--gpus',type=int,nargs='+',default=[0,1,2,3]);p.add_argument('--threads',type=int,default=16)
    p.add_argument('--orientation-chunk',type=int,default=256);p.add_argument('--cuda-lib-dir',type=Path)
    p.add_argument('--calibration',type=Path);p.add_argument('--calibrate',action='store_true')
    p.add_argument('--pilot-count',type=int,default=32768);p.add_argument('--kernel-count',type=int,default=8192)
    p.add_argument('--phi-points',type=int,default=4)
    p.add_argument('--gpu-block-size',type=int,choices=[64,128,256],default=64)
    p.add_argument('--gpu-warp-beams',choices=['auto','warp','thread'],default='auto')
    p.add_argument('--orientation-pipeline',action='store_true')
    p.add_argument('--depths',type=int,nargs='+',default=[8,10,12])
    p.add_argument('--level-cutoffs',nargs='+',default=['.001','off','off'])
    p.add_argument('--initial-counts',type=int,nargs='+',default=[131072,8192,8192])
    p.add_argument('--min-correction-count',type=int,default=8192)
    p.add_argument('--max-count',type=int,default=67108864)
    p.add_argument('--max-rounds',type=int,default=12)
    p.add_argument('--max-growth',type=int,choices=[2,4,8],default=4)
    p.add_argument('--relative-error',type=float,default=.03)
    p.add_argument('--allocation',choices=['cost','double'],default='cost')
    p.add_argument('--seeds',type=int,nargs=8,default=PRODUCTION_SEEDS)
    p.add_argument('--dry-run',action='store_true')
    a=p.parse_args()
    if not (len(a.depths)==len(a.level_cutoffs)==len(a.initial_counts)>=2):p.error('depths, cutoffs and counts must have equal length>=2')
    if a.depths!=sorted(a.depths) or min(a.depths)<=0:p.error('depths must be positive and nondecreasing')
    if a.level_cutoffs[-1]!='off':p.error('final level must have cutoff off')
    if min(*a.initial_counts,a.max_count,a.min_correction_count,a.max_rounds,a.phi_points)<=0:p.error('counts must be positive')
    if a.max_count>=2147483647 or max(a.initial_counts)>a.max_count or a.min_correction_count>a.max_count:
        p.error('initial/minimum counts must fit max-count and signed32bit count limit')
    if not 0<a.relative_error<1:p.error('relative error must lie in(0,1)')
    if len(set(a.seeds))!=8 or any(s<0 for s in a.seeds):p.error('eight distinct nonnegative seeds required')
    if set(a.seeds)&set(TRAIN_SEEDS+VALIDATION_SEEDS):p.error('production seeds must be disjoint from calibration pilots')
    for cutoff in a.level_cutoffs:
        if cutoff!='off' and not 0<float(cutoff)<1:p.error('cutoffs must be off or in(0,1)')
    try:validate_level_order(a.depths,a.level_cutoffs)
    except ValueError as error:p.error(str(error))
    runner=JobRunner(a.binary,a.output,case=a.case,height=a.height_um,diameter=a.diameter_um,
                     index=a.refractive_index,wave=a.wavelength_um,theta_file=a.theta_grid_file,
                     gpus=a.gpus,threads=a.threads,chunk=a.orientation_chunk,
                     cuda_lib_dir=a.cuda_lib_dir,mean_cache=a.mean_cache,controls=a.analytic_controls)
    options=dict(pipeline=a.orientation_pipeline,block=a.gpu_block_size,warp=a.gpu_warp_beams)
    calibration=None;weights=None;phi=a.phi_points
    if a.calibration:
        calibration=json.loads(a.calibration.read_text())
        if calibration['status']!='calibrated' or calibration['physical']!=runner.physical or calibration['binary_sha256']!=runner.binary_hash:
            p.error('calibration must match physical inputs and binary')
        if set(a.seeds)&set(calibration['training_seeds']+calibration['validation_seeds']):p.error('production/calibration seed overlap')
        phi=calibration['selected_phi'];options=calibration['selected_kernel']
        weights=Path(calibration['weights']) if calibration['weights'] is not None else None
    if a.dry_run:
        print(json.dumps(dict(physical=runner.physical,depths=a.depths,cutoffs=a.level_cutoffs,
                             initial_counts=a.initial_counts,phi=phi,options=options),indent=2));return
    lock=(runner.work/'controller.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if a.calibrate:
        if a.calibration:p.error('use either --calibrate or --calibration')
        from benchmark_gpu_optimized import calibrate
        calibration=calibrate(runner,a.pilot_count,(4,8,16),a.kernel_count,a.depths[-1])
        phi=calibration['selected_phi'];options=calibration['selected_kernel']
        weights=Path(calibration['weights']) if calibration['weights'] is not None else None
    config=dict(physical=runner.physical,binary_sha256=runner.binary_hash,depths=a.depths,
                cutoffs=a.level_cutoffs,seeds=a.seeds,phi=phi,options=options,
                weights_sha256=None if weights is None else hashlib.sha256(weights.read_bytes()).hexdigest(),
                tolerance=a.relative_error,allocation=a.allocation)
    state=dict(status='active',configuration=config,configuration_sha256=fingerprint(config),levels=[],
               controller_pid=os.getpid(),started=time.time(),sampling_policy=dict(max_growth=a.max_growth))
    status=runner.work/'status.json'
    if status.exists():
        previous=json.loads(status.read_text())
        if previous['configuration_sha256']!=state['configuration_sha256']:raise ValueError('resume inputs changed')
        if previous['status']=='converged':print('already converged',previous['output']);return
    atomic_json(status,state)
    counts=a.initial_counts[:];previous_mean=None;streak=0
    previous_coarse=None;previous_coarse_count=0;coarse_verified=False
    try:
        for round_number in range(a.max_rounds):
            levels=[];costs=[];wall=0.
            for index,n in enumerate(counts):
                upper,_,upper_time,_=runner.stage(n,a.depths[index],a.level_cutoffs[index],phi,
                                                  seeds=a.seeds,weights=weights,profile=True,**options)
                if index:
                    lower,_,lower_time,_=runner.stage(n,a.depths[index-1],a.level_cutoffs[index-1],phi,
                                                      seeds=a.seeds,weights=weights,profile=True,**options)
                    value=difference(upper,lower);level_time=upper_time+lower_time
                else:value=upper;level_time=upper_time
                levels.append(value);costs.append(level_time/(8*n));wall+=level_time
            data=combine(levels);metric,mean,ci=confidence(data,previous_mean,a.relative_error)
            coarse_change=None
            if previous_coarse is not None and counts[0]>previous_coarse_count:
                # Hold all current corrections fixed when checking the main
                # level: two changing corrections cannot hide a coarse error.
                coarse_change=float(np.max(abs(levels[0].mean(axis=0)[:,2:]-previous_coarse.mean(axis=0)[:,2:])/mean[:,2,None]))
                coarse_verified=coarse_verified or coarse_change<=a.relative_error
            predecessor_counts=None
            if metric['confidence_pass'] and not coarse_verified and counts[0]>=2*a.initial_counts[0]:
                # A very sparse warmup can be far from an already precise
                # estimate. Verify its nearest nested predecessor, rather than
                # forcing another large expensive forward sample merely to
                # compare it against an underresolved warmup.
                predecessor_count=counts[0]//2
                predecessor,_,_,_=runner.stage(predecessor_count,a.depths[0],a.level_cutoffs[0],phi,
                                               seeds=a.seeds,weights=weights,profile=True,**options)
                predecessor_counts=[predecessor_count,*counts[1:]]
                predecessor_mean=combine([predecessor,*levels[1:]]).mean(axis=0)
                metric,mean,ci=confidence(data,predecessor_mean,a.relative_error)
                coarse_change=metric['max_nested_change']
                coarse_verified=bool(metric['nested_pass'])
            streak=streak+1 if metric['confidence_pass'] and metric['nested_pass'] and coarse_verified else 0
            record=dict(round=round_number,counts=counts[:],stage_time_upper_bound=wall,
                        cost_per_orientation=costs,pass_streak=streak,
                        coarse_refinement_verified=coarse_verified,coarse_refinement_change=coarse_change,
                        refinement_baseline_counts=predecessor_counts,**metric)
            state['levels'].append(record);atomic_json(status,state);print(json.dumps(record),flush=True)
            if streak>=2:
                output=runner.work/'mueller_multilevel.dat'
                np.savetxt(output,mean,fmt='%.17g',comments='',
                    header='ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44')
                np.savetxt(runner.work/'per_angle_accuracy.csv',np.column_stack([mean[:,0],mean[:,2],ci[:,0],ci.max(axis=1)]),
                    delimiter=',',header='theta_deg,M11,estimated95_M11_relative,estimated95_max_mueller_scaled_M11',comments='')
                state.update(status='converged',output=str(output),final_precision=metric,completed=time.time(),
                             target_depth=a.depths[-1],target_cutoff=a.level_cutoffs[-1]);atomic_json(status,state);return
            if a.allocation=='cost' and not metric['confidence_pass']:
                proposed=allocation(levels,costs,counts,mean[:,2],a.relative_error*.85)
                counts=[max(n,min(a.max_count,a.max_growth*n,forecast),a.min_correction_count if i else 1)
                        for i,(n,forecast) in enumerate(zip(counts,proposed))]
                if counts==record['counts']:counts=[min(a.max_count,2*n) for n in counts]
            elif a.allocation=='cost':
                # Once the main estimate has passed its own N refinement,
                # verify the dominant physical correction without repeatedly
                # oversampling an already precise expensive coarse estimate.
                scores=[float(np.max(x[:,:,2:].var(axis=0,ddof=1)/mean[:,2,None]**2)) for x in levels]
                correction_indices=[i for i in range(1,len(scores)) if scores[i]>0]
                index=max(correction_indices,key=lambda i:scores[i]) if coarse_verified and correction_indices else 0
                counts[index]=min(a.max_count,2*counts[index])
            else:counts=[min(a.max_count,2*n) for n in counts]
            if counts==record['counts']:break
            previous_mean=mean
            previous_coarse=levels[0];previous_coarse_count=record['counts'][0]
        state['status']='requires_further_sampling';atomic_json(status,state)
    except Exception as error:
        state.update(status='failed',error=str(error));atomic_json(status,state);raise


if __name__=='__main__':main()
