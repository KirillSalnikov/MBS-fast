#!/usr/bin/env python3
"""Compare bounded GPU beam/facet scheduling on identical optical samples.

No campaign files or processes are modified. Every run saves its command,
environment, native outputs and timings; fallback cannot count as GPU speedup.
"""
import argparse, json, os, re, statistics, subprocess, time
from pathlib import Path
from test_droxtal_resident import load, difference, integral_difference

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary',type=Path,required=True)
    parser.add_argument('--before',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--tolerance',type=float,default=1e-7)
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--chunk',type=int,default=0)
    parser.add_argument('--modes',nargs='+',default=['before','legacy','compact','facets4','facets8','facets8_pooled'])
    a=parser.parse_args();root=a.output.resolve();root.mkdir(parents=True)
    repo=Path(__file__).resolve().parents[1]
    grid=root/'theta.txt';grid.write_text('0\n0.1\n170\n173\n176\n179\n179.9\n180\n')
    narrow=root/'narrow.txt';narrow.write_text('160\n160.01\n')
    common=['--method','po','--backend','cuda','--wavelength-um','.532','--ri','1.3116','0',
        '--threads','2','--haar-alpha','--theta-grid-file',str(grid),'--phi-points','8',
        '--orientation-pipeline','--profile-phases','--trace-cache-mib','0','--close','--allow-experimental-environment']
    if a.chunk:common+=['--orientation-chunk',str(a.chunk)]
    shapes={'column':['--particle','1','398.1','138.9','--symmetry','2','6'],
        'droxtal':['--particle','4','44.565','--symmetry','2','6'],
        'asymmetric':['--particle-file',str(repo/'tests/data/resident_asymmetric.particle'),'--resize-dmax-um','25','--symmetry','1','1']}
    records=[]
    def run(label,mode,command):
        dest=root/(label+'_'+mode);dest.mkdir()
        env=os.environ.copy();env.update(OMP_WAIT_POLICY='PASSIVE',MBS_HOST_MEM_BUDGET_MB='2048')
        env['MBS_GPU_TRACE_CONVEX_SCHEDULE']='legacy'if mode=='before'else mode.replace('_pooled','')
        env['MBS_GPU_TRACE_CONVEX_WORKSPACE_MIB']='64'if mode.endswith('_pooled')else'0'
        binary=a.before if mode=='before'else a.binary
        cmd=[str(binary)]+command+['--po-trace-backend','resident','--output',str(dest/'result')]
        (dest/'command.args').write_text('\n'.join(cmd)+'\n')
        (dest/'environment.json').write_text(json.dumps({k:v for k,v in env.items()if k.startswith(('MBS_','OMP_','CUDA_'))},indent=2))
        start=time.monotonic()
        with(dest/'stdout.log').open('w')as stream:
            result=subprocess.run(cmd,stdout=stream,stderr=subprocess.STDOUT,env=env,timeout=300)
        assert result.returncode==0,(label,mode,dest/'stdout.log')
        phases={}
        for line in (dest/'result/result_phase_timings.tsv').read_text().splitlines()[1:]:
            name,value,_=line.split();phases[name]=float(value)
        text=(dest/'stdout.log').read_text()
        accepted=sum(map(int,re.findall(r'GPU convex resident direct: accepted=(\d+)',text)))
        fallback=sum(map(int,re.findall(r'CPU_fallback=(\d+)',text)))
        return dict(directory=str(dest),seconds=time.monotonic()-start,phase=phases,accepted=accepted,fallback=fallback)
    def compare(label,command,order,count):
        runs={mode:run(label,mode,command)for mode in order};errors={}
        reference=Path(runs['before']['directory'])/'result'
        for mode in order:
            assert runs[mode]['accepted']==count and runs[mode]['fallback']==0,(label,mode,runs[mode])
            if mode=='before':continue
            current=Path(runs[mode]['directory'])/'result';err={}
            for f in sorted(reference.glob('result*.dat')):err[f.name]=difference(load(f),load(current/f.name))
            for f in sorted(reference.glob('result*_integrals.tsv')):err[f.name]=integral_difference(f,current/f.name)
            assert err and max(err.values())<=a.tolerance,(label,mode,err)
            errors[mode]=err
        record=dict(case=label,runs=runs,errors=errors);records.append(record)
        (root/'comparison.json').write_text(json.dumps(records,indent=2))
        print(json.dumps(dict(case=label,max_error=max(v for e in errors.values()for v in e.values()),
            processing={k:v['phase']['processing']for k,v in runs.items()},
            kernel={k:v['phase']['gpu_resident_trace_kernel']for k,v in runs.items()})),flush=True)
    for name,shape in shapes.items():
        depth=14 if name=='droxtal'else 8
        for incoherent in [False,True]:
            for cutoff in ['off','balanced']:
                cmd=common+shape+['--max-reflections',str(depth),'--paired-depth',str(depth-2),
                    '--orientation-prefix-output','64','--sobol-seed','128','23','--cutoff-profile',cutoff]
                if incoherent:cmd+=['--incoherent']
                compare(name+('_incoherent_'if incoherent else '_coherent_')+cutoff,cmd,
                    a.modes if cutoff=='off'else list(reversed(a.modes)),128)
    for name,shape in [('column',shapes['column']),('droxtal',shapes['droxtal']),
        ('column631',['--particle','1','631','174.8','--symmetry','2','6'])]:
        count=512 if name=='droxtal'else 4096;depth=14 if name=='droxtal'else(6 if name=='column631'else 8)
        for repeat in range(a.repeats):
            cmd=common+shape+['--max-reflections',str(depth),'--sobol-seed',str(count),'23','--cutoff-profile','off','--incoherent']
            cmd[cmd.index('--theta-grid-file')+1]=str(narrow)
            cmd[cmd.index('--phi-points')+1]='64'if name=='column631'else'16'
            compare('timing_'+name+'_'+str(repeat),cmd,
                a.modes if repeat%2==0 else list(reversed(a.modes)),count)
    summary={}
    for name in ['column','droxtal','column631']:
        selected=[r for r in records if r['case'].startswith('timing_'+name+'_')]
        if not selected:continue
        summary[name]={phase:{mode:statistics.median(r['runs'][mode]['phase'][phase]for r in selected)
            for mode in a.modes}for phase in ['processing','trace_prepare','gpu_resident_trace_kernel']}
    (root/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
