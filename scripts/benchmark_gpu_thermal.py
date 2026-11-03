#!/usr/bin/env python3
"""Minute-scale MBS workload with excluded thermal warmup and GPU-bound gates."""
import argparse
import json
import os
import re
import statistics
import threading
import time
from pathlib import Path
import benchmark_hardware as b

SUITE=Path(__file__).resolve().parents[1]/'benchmarks/gpu_thermal_v1'

def main():
    suite=json.loads((SUITE/'suite.json').read_text())
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['reference','run'])
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--prefix',type=Path)
    p.add_argument('--profile',choices=list(suite['profiles']),default=suite['default_profile'])
    p.add_argument('--reference',type=Path,default=SUITE/'reference/fp64.json')
    p.add_argument('--warmup-seconds',type=float,default=120)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--threads',type=int,default=8)
    p.add_argument('--device',type=int,default=0)
    p.add_argument('--timeout',type=int,default=1200)
    a=p.parse_args()
    if a.repeats<3 or a.warmup_seconds<120 or a.threads!=8 or a.device<0:
        p.error('Use >=3 repeats, >=120 s excluded warmup, exactly 8 threads and a nonnegative device')
    a.output=a.output.expanduser().absolute()
    if a.output.exists():p.error('Use a fresh output directory')
    if a.action=='reference':a.profile='gpu-fp64'
    a.binary=a.python=None;a.label='RTX GPU thermal workload'
    build=b.build_info(a,suite);hardware=b.hardware(a,build)
    if hardware['before'] and (hardware['before'].get('utilization.gpu') or 0)>10:
        raise RuntimeError('GPU busy before test; retry when idle')
    a.output.mkdir(parents=True);env=b.environment(a,build);case=suite['cases'][0]
    case_hash=b.digest(dict(suite=suite['id'],case=case,threads=a.threads))
    report=dict(schema=1,suite=suite['id'],profile=a.profile,threads=a.threads,
                time=time.time(),runner_sha256=b.sha(Path(__file__)),
                build={k:v for k,v in build.items() if k!='exports'},hardware=hardware,
                case=case,warmup_seconds_required=a.warmup_seconds,repeats=a.repeats,
                timing_protocol='fresh complete GPU-dominant MBS processes, >=120 s excluded thermal warmup, fixed 2048 orientations/2048 azimuths/181 theta points',
                cases=[],passed=False)
    def invoke(label):
        directory=a.output/label;command=b.command_for(case,build,a.profile,a.threads,directory)
        guard=[];stop=threading.Event()
        def monitor_processes():
            while not stop.is_set():
                value=b.capture([b.smi_path(),'-i',str(a.device),'--query-compute-apps=pid','--format=csv,noheader,nounits'])
                entry=dict(time=time.time(),processes=[],forbidden=[])
                if value is None:entry['error']='GPU compute process query unavailable'
                for token in (value or '').splitlines():
                    pid=token.strip()
                    if not pid.isdigit():
                        entry['error']='Unrecognized GPU process query output'
                        continue
                    try:
                        status=Path('/proc/'+pid+'/stat').read_text().split()
                        executable=str(Path('/proc/'+pid+'/exe').resolve())
                        ours=executable==build['native_binary'] and int(status[3])==os.getpid()
                        entry['processes'].append(dict(pid=int(pid),state=status[2],start_ticks=status[21],benchmark=ours))
                        if not ours:entry['forbidden'].append(int(pid))
                    except FileNotFoundError:
                        continue
                guard.append(entry);stop.wait(1)
        thread=threading.Thread(target=monitor_processes,daemon=True);thread.start()
        try:run=b.timed_process(command,env,directory,a.device,a.timeout)
        finally:stop.set();thread.join(timeout=35)
        b.write(directory/'compute_process_guard.json',guard)
        if not guard or any(v.get('error') or v['forbidden'] for v in guard):
            raise RuntimeError('Concurrent/unidentified GPU compute process; hardware score rejected')
        run['compute_process_guard']=guard
        matrix,metrics,work=b.read_result(case,directory,a.profile)
        grids=set(re.findall(r'nAz=(\d+) nZen=(\d+)',(directory/'stdout.log').read_text()))
        if grids!={('2048','180')}:raise ValueError('Native angular GPU work does not match the frozen case')
        if b.sha(build['native_binary'])!=build['binary_sha256']:raise RuntimeError('Binary changed during test')
        run.update(command=command,metrics=metrics,work=work,matrix_sha256=b.digest(matrix))
        return run,matrix
    if a.action=='reference':
        run,matrix=invoke('reference_generation')
        b.write(a.output/'fp64.json',dict(schema=1,suite=suite['id'],case_sha256=case_hash,
                matrix=matrix,matrix_sha256=b.digest(matrix),build=report['build'],generation=run,
                note='Frozen FP64 numerical curve; not a physical-convergence certificate'))
        print('REFERENCE',a.output/'fp64.json',flush=True);return
    reference=json.loads(a.reference.read_text())
    if reference['case_sha256']!=case_hash or b.digest(reference['matrix'])!=reference['matrix_sha256']:
        raise ValueError('Frozen reference/workload mismatch')
    limits=suite['profiles'][a.profile]['validation']
    record=dict(id=case['id'],comparison_key=b.digest([case_hash,a.profile,a.threads,reference['matrix_sha256'],report['timing_protocol'],a.warmup_seconds]),reference_sha256=b.sha(a.reference),runs=[],passed=False)
    report['cases'].append(record)
    warm=0.;index=0;measured_seconds=0.
    while warm<a.warmup_seconds or index<a.repeats or measured_seconds<180:
        warming=warm<a.warmup_seconds
        label='warmup_'+str(len(record['runs'])) if warming else 'repeat_'+str(index+1)
        run,matrix=invoke(label);error=b.matrix_error(matrix,reference['matrix'])
        if any(error[k]>v for k,v in limits.items()):raise ValueError('Full Mueller matrix gate failed: '+str(error))
        run.update(name=label,warmup=warming,validation=error)
        run['native_gpu_fraction']=run['metrics']['gpu_kernels_s']/run['wall_s']
        utilization=[s['utilization.gpu'] for s in run['telemetry'][2:-2] if s.get('utilization.gpu') is not None]
        if not utilization:raise RuntimeError('GPU utilization telemetry is required')
        run['median_gpu_utilization']=statistics.median(utilization)
        record['runs'].append(run)
        if warming:warm+=run['wall_s']
        else:
            if run['native_gpu_fraction']<.90 or run['median_gpu_utilization']<90:
                b.write(a.output/'results.json',report)
                raise RuntimeError('Measured workload is not GPU dominated; not eligible as a GPU score')
            index+=1;measured_seconds+=run['wall_s']
        b.write(a.output/'results.json',report)
        print(label,f"{run['wall_s']:.3f}s CUDA fraction={run['native_gpu_fraction']:.4f} GPU median={run['median_gpu_utilization']}%",flush=True)
    measured=[r for r in record['runs'] if not r['warmup']]
    if sum(r['wall_s'] for r in measured)<180:raise RuntimeError('Measured interval must span at least three minutes')
    record['statistics']={k:b.statistics_of([r[k] for r in measured]) for k in ['wall_s','cpu_user_s','cpu_system_s','peak_child_rss_MiB','native_gpu_fraction','median_gpu_utilization']}
    for key in measured[0]['metrics']:record['statistics'][key]=b.statistics_of([r['metrics'][key] for r in measured])
    record.update(passed=True,timing_stable=record['statistics']['wall_s']['relative_MAD']<=.1)
    report.update(passed=True,warmup_elapsed_s=warm,repeats=index,warmups=sum(r['warmup'] for r in record['runs']))
    b.write(a.output/'results.json',report);print('VALIDATED',a.output/'results.json',flush=True)

if __name__=='__main__':main()
