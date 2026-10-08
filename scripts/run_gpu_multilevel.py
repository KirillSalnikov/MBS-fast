from pathlib import Path
import concurrent.futures,hashlib,json,math,os,subprocess,time,fcntl
import numpy as np

import argparse
parser=argparse.ArgumentParser(description='GPU full-coherent Mueller multilevel averaging with physical analytic facet controls.')
parser.add_argument('--binary',type=Path,required=True)
parser.add_argument('--case',choices=['test1','test2'],default='test2')
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--theta-grid-file',type=Path)
parser.add_argument('--height-um',type=float)
parser.add_argument('--diameter-um',type=float)
parser.add_argument('--refractive-index',type=float,nargs=2,default=[1.3116,0.])
parser.add_argument('--wavelength-um',type=float,default=.532)
parser.add_argument('--phi-points',type=int,default=4)
parser.add_argument('--gpus',type=int,nargs=4,default=[0,1,2,3])
parser.add_argument('--threads',type=int,default=16)
parser.add_argument('--orientation-chunk',type=int,default=256)
parser.add_argument('--seeds',type=int,nargs=8,default=[11,23,37,53,71,89,107,131])
parser.add_argument('--coarse-counts',type=int,nargs='+',default=[131072,524288,2097152,8388608,33554432,67108864])
parser.add_argument('--coarse-reflections',type=int,default=8)
parser.add_argument('--fine-reflections',type=int,default=12)
parser.add_argument('--min-correction-count',type=int,default=32768)
parser.add_argument('--coarse-cutoff',default='.001')
parser.add_argument('--cuda-lib-dir',type=Path)
parser.add_argument('--dry-run',action='store_true')
args=parser.parse_args()
height=args.height_um if args.height_um is not None else (446.7 if args.case=='test1' else 891.3)
diameter=args.diameter_um if args.diameter_um is not None else (138.9 if args.case=='test1' else 207.8)
if any(not math.isfinite(v) for v in [height,diameter,args.wavelength_um,*args.refractive_index]):parser.error('physical values must be finite')
if args.min_correction_count<1:parser.error('minimum correction count must be positive')
if min(height,diameter,args.wavelength_um,args.phi_points,args.threads,args.orientation_chunk)<=0:parser.error('sizes, counts and wavelength must be positive')
if args.refractive_index[0]<1 or args.refractive_index[1]<0:parser.error('require real index>=1 and imaginary index>=0')
if len(set(args.gpus))!=4 or min(args.gpus)<0 or len(set(args.seeds))!=8:parser.error('use four distinct GPU indices and eight distinct scrambles')
if len(args.coarse_counts)<2 or any(n<=0 or n>=2147483647 for n in args.coarse_counts) or args.coarse_counts!=sorted(set(args.coarse_counts)):parser.error('coarse counts must be positive, increasing and distinct')
if not args.fine_reflections>=args.coarse_reflections>0:parser.error('fine reflection depth must be at least the positive coarse depth')
WORK=args.output.resolve();WORK.mkdir(parents=True,exist_ok=True)
BIN=args.binary.resolve();SEEDS=args.seeds
if args.theta_grid_file:theta_file=args.theta_grid_file.resolve()
else:
    theta_file=WORK/'theta.txt';rows=267 if args.case=='test1' else 531
    if not theta_file.exists():theta_file.write_text(''.join('{:.14g}\n'.format(25*i/(rows-1)) for i in range(rows)))
rows=len(np.atleast_1d(np.loadtxt(theta_file)));mean_cache=WORK/'physical_means.cache'
configuration=dict(height_um=height,diameter_um=diameter,index=args.refractive_index,wavelength_um=args.wavelength_um,phi=args.phi_points,seeds=SEEDS,coarse_reflections=args.coarse_reflections,fine_reflections=args.fine_reflections,coarse_cutoff=args.coarse_cutoff,theta_sha256=hashlib.sha256(theta_file.read_bytes()).hexdigest())
configuration_sha256=hashlib.sha256(json.dumps(configuration,sort_keys=True).encode()).hexdigest()
if args.dry_run:print(json.dumps(configuration,indent=2));raise SystemExit(0)
status=WORK/'status.json'
state=dict(status='active',controller_pid=os.getpid(),started=time.time(),jobs=[],levels=[],test=args.case,configuration_sha256=configuration_sha256,configuration=configuration,
    phi=args.phi_points,binary_sha256=hashlib.sha256(BIN.read_bytes()).hexdigest())
if status.exists():
    previous=json.load(open(status));assert previous['binary_sha256']==state['binary_sha256'] and previous['configuration_sha256']==configuration_sha256
    state['jobs']=previous['jobs'];state['levels']=previous['levels']

def save():
    temporary=status.with_suffix('.tmp');temporary.write_text(json.dumps(state,indent=2));temporary.replace(status)

def execute(gpu_slot,N,seed,depth,cutoff):
    gpu=args.gpus[gpu_slot]
    name=f'{args.case}_N{N}_seed{seed}_n{depth}_{cutoff}_phi{args.phi_points}';target=WORK/name
    matches=[j for j in state['jobs'] if j['name']==name and j.get('returncode')==0]
    if matches:
        job=matches[-1];assert Path(job['output']).exists();return job
    command=[str(BIN),'--method','po','--backend','cuda','--particle','1',str(height),str(diameter),
        '--refractive-index',str(args.refractive_index[0]),str(args.refractive_index[1]),'--wavelength-um',str(args.wavelength_um),'--max-reflections',str(depth),
        '--sobol-seed',str(N),str(seed),'--haar-alpha','--analytic-facet-average','--analytic-shadow-control','circular',
        '--analytic-mean-cache',str(mean_cache),
        '--theta-grid-file',str(theta_file),'--phi-points',str(args.phi_points),'--threads',str(args.threads),
        '--orientation-chunk',str(args.orientation_chunk),'--progress-interval','10','--output',str(target),'--close']
    command+=['--cutoff-profile','off'] if cutoff=='off' else ['--beam-cutoff',cutoff]
    env={k:v for k,v in os.environ.items() if not k.startswith('MBS_')}
    env['CUDA_VISIBLE_DEVICES']=str(gpu)
    if args.cuda_lib_dir:env['LD_LIBRARY_PATH']=str(args.cuda_lib_dir)+':'+env.get('LD_LIBRARY_PATH','')
    start=time.monotonic();job=dict(name=name,N=N,seed=seed,gpu=gpu,depth=depth,cutoff=cutoff,command=command,status='running')
    path=WORK/(name+'.job.json')
    with (WORK/(name+'.log')).open('w') as log:
        process=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT);job['pid']=process.pid
        path.write_text(json.dumps(job,indent=2));rc=process.wait()
    job.update(returncode=rc,status='finished' if rc==0 else 'failed',seconds=time.monotonic()-start,
        output=str(target/(name+'.dat')));path.write_text(json.dumps(job,indent=2));assert rc==0,name
    return job

def stage(N,depth=None,cutoff=None):
    if depth is None:depth=args.coarse_reflections
    if cutoff is None:cutoff=args.coarse_cutoff
    jobs=[]
    for start in [0,4]:
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            futures=[pool.submit(execute,gpu,N,SEEDS[start+gpu],depth,cutoff) for gpu in range(4)]
            for future in concurrent.futures.as_completed(futures):
                job=future.result();jobs.append(job)
                if not any(j['name']==job['name'] for j in state['jobs']):state['jobs'].append(job)
                save();print('finished',job['name'],job['seconds'],flush=True)
    data=np.array([np.loadtxt(next(j['output'] for j in jobs if j['seed']==seed),skiprows=1) for seed in SEEDS])
    assert data.shape==(8,rows,18) and np.isfinite(data).all()
    assert np.all(data[:,:,0]==data[0,:,0])
    return data

def metric(data,previous=None):
    mean=data.mean(axis=0);se=data[:,:,2:].std(axis=0,ddof=1)/math.sqrt(8)
    scale=np.maximum(abs(mean[:,2,None]),1e-100);ci=2.365*se/scale
    change=None if previous is None else abs(mean[:,2:]-previous[:,2:])/scale
    return dict(max_M11_pointwise95=float(ci[:,0].max()),max_mueller95_scaled_M11=float(ci.max()),
        confidence_pass=bool(np.all(ci<=.03)),max_nested_change=None if change is None else float(change.max()),
        nested_pass=previous is not None and bool(np.all(change<=.03))),mean

def combine(coarse,fine,small):
    assert np.array_equal(coarse[:,:,0],fine[:,:,0]) and np.array_equal(fine[:,:,0],small[:,:,0])
    result=coarse.copy();result[:,:,2:]+=fine[:,:,2:]-small[:,:,2:];return result

lock=(WORK/'campaign.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);save()
try:
    previous_data=None;previous_mean=None
    for N in args.coarse_counts:
        coarse=stage(N);coarse_metric,coarse_mean=metric(coarse,previous_mean)
        coarse_metric.update(kind='coarse',N=N);state['levels'].append(coarse_metric);save();print('coarse',coarse_metric,flush=True)
        if previous_data is not None and coarse_metric['confidence_pass']:
            C=max(args.min_correction_count,N//64);small=stage(C);fine=stage(C,args.fine_reflections,'off')
            before_metric,before_mean=metric(combine(previous_data,fine,small))
            current_metric,current_mean=metric(combine(coarse,fine,small),before_mean)
            streak=1 if current_metric['confidence_pass'] and current_metric['nested_pass'] else 0
            state['levels'].append(dict(kind='multilevel',coarse_N=N,correction_N=C,pass_streak=streak,**current_metric));save()
            for refinedC in [4*C,16*C,64*C]:
                small=stage(refinedC);fine=stage(refinedC,args.fine_reflections,'off');data=combine(coarse,fine,small)
                assessment,mean=metric(data,current_mean)
                streak=streak+1 if assessment['confidence_pass'] and assessment['nested_pass'] else 0
                record=dict(kind='multilevel',coarse_N=N,correction_N=refinedC,pass_streak=streak,**assessment)
                state['levels'].append(record);save();print('multilevel',record,flush=True)
                current_mean=mean
                if streak>=2:
                    np.savetxt(WORK/'mueller_multilevel.dat',mean,fmt='%.17g',comments='',
                        header='ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44')
                    state.update(status='converged',completed=time.time(),coarse_N=N,correction_N=refinedC,
                        output=str(WORK/'mueller_multilevel.dat'),final_precision=assessment,
                        target=dict(height_um=height,diameter_um=diameter,index=args.refractive_index,wavelength_um=args.wavelength_um,
                            max_reflections=args.fine_reflections,cutoff='off',phi=args.phi_points,full_Haar=True))
                    save();break
            if state['status']=='converged':break
        previous_data=coarse;previous_mean=coarse_mean
    else:state['status']='requires_further_sampling';save()
except Exception as error:
    state.update(status='failed',error=str(error));save();raise
