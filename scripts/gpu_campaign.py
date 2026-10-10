"""Reusable four-GPU jobs, frozen control calibration and multilevel statistics."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import math
import os
import queue
import subprocess
import time

import numpy as np

T95_7 = 2.3646242510102993
PRODUCTION_SEEDS = [11, 23, 37, 53, 71, 89, 107, 131]
TRAIN_SEEDS = [1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049]
VALIDATION_SEEDS = [2003, 2011, 2017, 2027, 2029, 2039, 2053, 2063]


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def confidence(samples, previous=None, tolerance=.03):
    """Per-seed combined estimators; retain covariance between every level."""
    if samples.shape[0] != 8 or samples.ndim != 3 or samples.shape[2] != 18:
        raise ValueError('require eight scrambles, angular rows and18 columns')
    if not np.isfinite(samples).all() or not np.all(samples[:, :, :2] == samples[0, :, :2]):
        raise ValueError('nonfinite or mismatched angular grids')
    mean = samples.mean(axis=0)
    if np.any(mean[:, 2] <= 0):
        raise ValueError('estimated M11 must be positive')
    ci = T95_7 * samples[:, :, 2:].std(axis=0, ddof=1) / math.sqrt(8) / mean[:, 2, None]
    change = None if previous is None else abs(mean[:, 2:] - previous[:, 2:]) / mean[:, 2, None]
    metrics = dict(max_M11_pointwise95=float(ci[:, 0].max()),
                   max_mueller95_scaled_M11=float(ci.max()),
                   confidence_pass=bool(np.all(ci <= tolerance)),
                   max_nested_change=None if change is None else float(change.max()),
                   nested_pass=change is not None and bool(np.all(change <= tolerance)))
    return metrics, mean, ci


def estimate_cost_score(data, seconds):
    metric, _, _ = confidence(data)
    return seconds * metric['max_mueller95_scaled_M11'] ** 2


def reweight(data, summaries, coefficients):
    """Reconstruct M11 from raw Y and its two controls; leave all other entries."""
    result = data.copy()
    wr, ws = coefficients.T
    for i, summary in enumerate(summaries):
        result[i, :, 2] = (summary['raw_M11']
                           + wr * (summary['reflection_mean'] - summary['sampled_reflection'])
                           + ws * (summary['shadow_mean'] - summary['sampled_shadow']))
    return result


def fit_control_weights(training, validation, summaries_train, summaries_validation):
    """Independent pilot with regularized two-control covariance per angle.

    The validation pilot selects between fixed b=1 and the training fit.
    Production must use a third disjoint set of scrambles. No production
    samples enter either fitting or model selection.
    """
    rows = training.shape[1]
    coefficients = np.ones((rows, 2))
    for row in range(rows):
        y = np.array([s['raw_M11'][row] for s in summaries_train])
        c = np.array([[s['sampled_reflection'][row], s['sampled_shadow'][row]] for s in summaries_train])
        c -= c.mean(axis=0)
        scale = np.sqrt((c*c).sum(axis=0)/7)
        active = scale > max(abs(y.mean()), 1e-100)*1e-12
        if not active.any():
            continue
        x = c[:, active] / scale[active]
        cov = x.T @ x / 7
        cy = x.T @ (y-y.mean()) / 7
        beta = np.linalg.solve(cov + .1*np.eye(len(cy)), cy) / scale[active]
        candidate = np.ones(2)
        candidate[active] = np.clip(beta, -4, 4)
        coefficients[row] = candidate
    proposed = reweight(validation, summaries_validation, coefficients)
    default = reweight(validation, summaries_validation, np.ones_like(coefficients))
    candidate_var = proposed[:, :, 2].var(axis=0, ddof=1)
    default_var = default[:, :, 2].var(axis=0, ddof=1)
    # A substantial held-out gain is required before departing from b=1.
    use = candidate_var < .7 * default_var
    coefficients[~use] = 1.
    return coefficients, dict(rows_fitted=int(use.sum()), rows=rows,
                             validation_selection='candidate variance<0.7*default; production seeds disjoint',
                             regularization=.1, coefficient_bound=4.)


def write_weights(path, theta, coefficients):
    np.savetxt(path, np.column_stack([theta, coefficients]), fmt='%.17g', comments='',
               header='theta_deg reflection_weight shadow_weight')


def allocation(level_samples, costs, current_counts, reference_m11, tolerance=.03):
    """Classical V/c allocation as a forecast; actual stopping uses confidence.

    Owen QMC variances need not scale exactly as1/N. Allocation therefore
    proposes budgets only. Correlation is retained by the final combined
    per-seed estimator and both consecutive refinement checks.
    """
    variance = np.stack([x[:, :, 2:].var(axis=0, ddof=1)*n
                         for x, n in zip(level_samples, current_counts)])
    costs = np.asarray(costs, dtype=float)
    if np.any(costs <= 0) or not np.isfinite(variance).all():
        raise ValueError('invalid level cost or variance')
    target = (tolerance*np.asarray(reference_m11)/T95_7)**2*8
    root = np.sqrt(variance*costs[:, None, None])
    ni = np.sqrt(variance/costs[:, None, None])*root.sum(axis=0)[None]/target[None, :, None]
    return [max(1, int(2**math.ceil(math.log2(max(float(x.max()), 1))))) for x in ni]


class JobRunner:
    def __init__(self, binary, output, case='test1', height=None, diameter=None,
                 index=(1.3116, 0.), wave=.532, theta_file=None, gpus=(0,1,2,3),
                 threads=16, chunk=256, cuda_lib_dir=None, mean_cache=None, controls='facets', scheduling='waves'):
        self.binary = Path(binary).resolve()
        self.work = Path(output).resolve()
        self.work.mkdir(parents=True, exist_ok=True)
        self.case = case
        presets={'test1':(446.7,138.9,267,0.,25.),
                 'test2':(891.3,207.8,531,0.,25.),
                 'test3':(316.2,123.8,227,170.,180.)}
        if case not in presets:raise ValueError('unknown test case')
        preset_height,preset_diameter,count,theta_min,theta_max=presets[case]
        self.height = height if height is not None else preset_height
        self.diameter = diameter if diameter is not None else preset_diameter
        self.index, self.wave = list(index), wave
        self.gpus, self.threads, self.chunk = list(gpus), threads, chunk
        self.cuda_lib_dir = cuda_lib_dir
        if controls not in ('facets','off'):
            raise ValueError('controls must be facets or off')
        self.controls=controls
        if scheduling not in ('waves','queue'):raise ValueError('unknown GPU scheduling mode')
        self.scheduling=scheduling
        if not self.gpus or len(set(self.gpus)) != len(self.gpus) or min(self.gpus) < 0:
            raise ValueError('require distinct nonnegative GPUs')
        if min(self.height, self.diameter, self.wave, threads, chunk) <= 0:
            raise ValueError('sizes and counts must be positive')
        if not all(math.isfinite(x) for x in [self.height,self.diameter,self.wave,*self.index]):
            raise ValueError('physical inputs must be finite')
        self.theta_file = Path(theta_file).resolve() if theta_file else self.work/'theta.txt'
        if theta_file is None:
            expected = ''.join(f'{theta_min+(theta_max-theta_min)*i/(count-1):.14g}\n' for i in range(count))
            if self.theta_file.exists() and self.theta_file.read_text() != expected:
                raise ValueError('existing theta grid differs from preset')
            self.theta_file.write_text(expected)
        self.theta = np.atleast_1d(np.loadtxt(self.theta_file))
        self.mean_cache = Path(mean_cache).resolve() if mean_cache else self.work/'physical_means.cache'
        self.mean_cache.parent.mkdir(parents=True, exist_ok=True)
        self.physical = dict(height_um=self.height,diameter_um=self.diameter,
                             index=self.index,wavelength_um=self.wave,
                             theta_sha256=hashlib.sha256(self.theta_file.read_bytes()).hexdigest())
        if controls!='facets':self.physical['analytic_controls']=controls
        self.binary_hash = hashlib.sha256(self.binary.read_bytes()).hexdigest()

    def job(self, gpu, n, seed, depth, cutoff, phi=4, pipeline=False, block=64,
            warp='auto', weights=None, profile=False):
        if min(n,depth,phi)<=0 or n>=2147483647 or seed<0 or block not in [64,128,256] or warp not in ['auto','warp','thread']:
            raise ValueError('invalid job counts or kernel configuration')
        if hashlib.sha256(self.binary.read_bytes()).hexdigest()!=self.binary_hash:
            raise RuntimeError('binary changed during campaign; use an immutable build directory')
        config = dict(**self.physical, binary_sha256=self.binary_hash, N=n,seed=seed,
                      depth=depth,cutoff=cutoff,phi=phi,pipeline=pipeline,block=block,warp=warp,
                      weights_sha256=None if weights is None else hashlib.sha256(Path(weights).read_bytes()).hexdigest(),
                      threads=self.threads,chunk=self.chunk,profile=profile)
        digest = fingerprint(config)
        name=f'{self.case}_N{n}_s{seed}_d{depth}_phi{phi}_{digest[:12]}'
        folder=self.work/'jobs'/name
        folder.mkdir(parents=True,exist_ok=True)
        metadata=folder/'job.json'
        output=folder/name
        if metadata.exists():
            old=json.loads(metadata.read_text())
            if old.get('configuration_sha256') != digest:
                raise ValueError('job cache input mismatch')
            if old.get('status')=='finished':
                data_path=Path(old['data'])
                if hashlib.sha256(data_path.read_bytes()).hexdigest()!=old['data_sha256']:
                    raise ValueError('cached matrix changed')
                return old
            if old.get('status')=='running':
                try:os.kill(old['pid'],0)
                except ProcessLookupError:pass
                else:raise RuntimeError('cached job is still running; avoid duplicate launch')
        command=[str(self.binary),'--method','po','--backend','cuda','--particle','1',str(self.height),str(self.diameter),
                 '--refractive-index',*[str(x) for x in self.index],'--wavelength-um',str(self.wave),
                 '--max-reflections',str(depth),'--sobol-seed',str(n),str(seed),'--haar-alpha',
                 '--theta-grid-file',str(self.theta_file),
                 '--phi-points',str(phi),'--threads',str(self.threads),'--orientation-chunk',str(self.chunk),
                 '--output',str(output),'--close']
        if self.controls=='facets':
            command += ['--analytic-facet-average','--analytic-shadow-control','circular',
                        '--analytic-mean-cache',str(self.mean_cache)]
        elif weights:
            raise ValueError('control weights require analytic controls')
        command += ['--cutoff-profile','off'] if cutoff=='off' else ['--beam-cutoff',str(cutoff)]
        if pipeline:command.append('--orientation-pipeline')
        if weights:command += ['--analytic-control-weights',str(Path(weights).resolve())]
        if profile:command.append('--profile-phases')
        env={k:v for k,v in os.environ.items() if not k.startswith('MBS_')}
        env.update(CUDA_VISIBLE_DEVICES=str(gpu),MBS_GPU_BLOCK=str(block))
        if warp!='auto':env['MBS_GPU_WARP_BEAMS']='1' if warp=='warp' else '0'
        if profile:env['MBS_GPU_TIMING']='1'
        if self.cuda_lib_dir:env['LD_LIBRARY_PATH']=str(self.cuda_lib_dir)+':'+env.get('LD_LIBRARY_PATH','')
        job=dict(configuration=config,configuration_sha256=digest,command=command,
                 gpu=gpu,status='running',started=time.time())
        start=time.monotonic()
        with (folder/'stdout.log').open('w') as log:
            process=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            job['pid']=process.pid;atomic_json(metadata,job)
            rc=process.wait()
        path=output/(name+'.dat')
        job.update(seconds=time.monotonic()-start,returncode=rc,status='finished' if rc==0 else 'failed',
                   data=str(path),summary=str(output/(name+'_analytic_facets.tsv')),
                   phase_timings=str(output/(name+'_phase_timings.tsv')))
        if rc==0:
            data=np.loadtxt(path,skiprows=1)
            if data.shape!=(len(self.theta),18) or not np.isfinite(data).all():
                raise ValueError('invalid output matrix: '+name)
            log_text=(folder/'stdout.log').read_text()
            for label in ['Hard tree-limit hits:','Orientations still incomplete:']:
                line=next((line for line in log_text.splitlines() if line.startswith(label)),None)
                if line is None or int(line.split(':',1)[1].strip())!=0:
                    raise ValueError('missing or nonzero limit diagnostic: '+name)
            job['data_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
        atomic_json(metadata,job)
        if rc:raise RuntimeError('GPU job failed; see '+str(folder/'stdout.log'))
        return job

    def stage(self, n, depth=8, cutoff='.001', phi=4, seeds=PRODUCTION_SEEDS, **options):
        if len(seeds)!=8 or len(set(seeds))!=8:
            raise ValueError('require eight distinct independent scrambles')
        jobs=[];wall=0.
        if self.scheduling=='queue':
            pending=queue.Queue();ordered=[None]*len(seeds)
            for index,seed in enumerate(seeds):pending.put((index,seed))
            def worker(gpu):
                busy=0.
                while True:
                    try:index,seed=pending.get_nowait()
                    except queue.Empty:return busy
                    result=self.job(gpu,n,seed,depth,cutoff,phi,**options)
                    ordered[index]=result;busy+=result['seconds']
            with concurrent.futures.ThreadPoolExecutor(max_workers=len(self.gpus)) as pool:
                busy=list(pool.map(worker,self.gpus))
            jobs=ordered
            # Cached jobs can be returned instantly by a single worker.
            # Preserve their recorded assignment when forecasting level cost.
            recorded_busy={gpu:0. for gpu in self.gpus}
            for result in jobs:
                slot=result.get('assigned_slot',result.get('gpu'))
                if slot in recorded_busy:recorded_busy[slot]+=result['seconds']
            wall=max(recorded_busy.values()) if any(recorded_busy.values()) else max(busy)
        else:
            for start in range(0,len(seeds),len(self.gpus)):
                slots=range(min(len(self.gpus),len(seeds)-start))
                with concurrent.futures.ThreadPoolExecutor(max_workers=len(self.gpus)) as pool:
                    batch=list(pool.map(lambda slot:self.job(self.gpus[slot],n,seeds[start+slot],depth,cutoff,phi,**options),slots))
                jobs += batch;wall += max(job['seconds'] for job in batch)
        data=np.stack([np.loadtxt(j['data'],skiprows=1) for j in jobs])
        if self.controls=='facets':
            summaries=[np.atleast_1d(np.genfromtxt(j['summary'],names=True,delimiter='\t')) for j in jobs]
        else:
            # Zero controls represent the unchanged full-integral estimator.
            dtype=[(key,float) for key in ('raw_M11','reflection_mean','shadow_mean',
                                         'sampled_reflection','sampled_shadow')]
            summaries=[]
            for matrix in data:
                summary=np.zeros(len(self.theta),dtype=dtype)
                summary['raw_M11']=matrix[:,2];summaries.append(summary)
        return data,summaries,wall,jobs
