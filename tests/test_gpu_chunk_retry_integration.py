#!/usr/bin/env python3
"""Actual FullAutoWorkers controller, injected statuses, real CPU optical output."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--probe',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
source=Path(__file__).resolve().parents[1]/'src/FullAuto.cpp'
code=source.read_text();carry=code.index('auto workingGpuChunkCaps=std::move(e.gpuChunkRetryCaps);')
reset=code.index('controller.reset(new Engine(args,config));controller->Init(',carry)
restore=code.index('controller->gpuChunkRetryCaps=std::move(workingGpuChunkCaps);',reset)
refinement=code.index('if(controller->refinementWindow)PrepareRefinementWindow(*controller);',restore)
assert carry<reset<restore<refinement
subprocess.run([str(a.probe.resolve()),'--promotion-lifecycle'],check=True,timeout=15)
worker=root/'fault_worker.py';shutil.copy2(Path(__file__).with_name('gpu_chunk_fault_worker.py'),worker);worker.chmod(0o755)
report=[];reference=None
for memory in (0,1):
 for name,sequence,chunk,success,retries in [('reference','',1024,True,0),('once','76',1024,True,1),('small_explicit','76',64,True,1),('implicit','76',0,True,1),('twice','76,76',1024,False,1),('memory_then_chunk','75,76',1024,False,0),('chunk_then_memory','76,75',1024,False,1),('unknown','77',1024,False,0),('physics','1',1024,False,0),('eof','eof',1024,False,0)]+([('partial_frame','partial76',1024,False,0)] if memory else []):
  case=root/f'{memory}_{name}';marker=root/(case.name+'.marker')
  env=dict(os.environ,TASK_REAL_WORKER=str(a.binary.resolve()),TASK_GPU_RETRY_MARKER=str(marker),TASK_GPU_RETRY_SEQUENCE=sequence)
  if success and retries:env['TASK_EXPECT_GPU_BACKOFF']='1'
  with (root/(case.name+'.log')).open('w') as out:
   result=subprocess.run([str(a.probe.resolve()),str(case),str(worker),str(memory),str(chunk)],env=env,stdout=out,stderr=subprocess.STDOUT,timeout=45)
  assert result.returncode==(0 if success else 1),(name,memory,result.returncode,(root/(case.name+'.log')).read_text()[-1800:])
  history=[json.loads(line) for line in (case/'fullauto_history.jsonl').read_text().splitlines()]
  events=[e for e in history if e['event']=='worker_gpu_chunk_retry'];assert len(events)==retries,(name,memory,events)
  attempts=[json.loads(line) for line in marker.with_suffix('.requests.jsonl').read_text().splitlines()] if marker.with_suffix('.requests.jsonl').exists() else []
  assert len(attempts)<=2,(name,'unbounded retry')
  if retries:
   event=events[0]['data'];assert 0<event['chunk_cap']<=min(128,chunk or 128)
   assert event['same_sampling'] and event['failed_outputs_quarantined'] and event['backoff_scope']=='current_controller_worker_identity'
   assert event['failed_attempt_seconds']>0 and event['total_attempt_seconds']>=event['failed_attempt_seconds']
  if success:
   values=[[float(v) for v in line.split()] for line in (case/'matrix.txt').read_text().splitlines()]
   if reference is None:reference=values
   assert values==reference,(name,memory,'matrix changed')
   if retries:
    assert int((case/'successful_worker.pid').read_text())!=attempts[0]['pid'],(name,'failed worker reused')
    complete=[e['data'] for e in history if e['event']=='worker_gpu_chunk_retry_complete'];assert len(complete)==1
    assert complete[0]['total_attempt_seconds']>=complete[0]['successful_attempt_seconds']
    final=(case/'first.command').read_text().splitlines()[1:]
    original=attempts[0]['words'][3:]
    def physics(words):
     words=list(words)
     while '--orientation-chunk' in words:
      i=words.index('--orientation-chunk');del words[i:i+2]
     return words
    assert physics(final)==physics(original),(name,'scientific argv changed')
    for path in (case/'first.command',case/'next.command'):
     words=path.read_text().splitlines();assert int(words[words.index('--orientation-chunk')+1])<=min(128,chunk or 128)
    if memory:
     files=(case/'virtual_paths.txt').read_text().splitlines()
     assert any('.result.failed/result.dat' in f for f in files) and any('.stdout.log' in f and 'gpu-chunk-retry_' in f for f in files)
    else:
     failed=list((case/'jobs').glob('*/gpu-chunk-retry_*.result.failed/result.dat'))
     assert len(failed)==1 and 'POISON INCOMPLETE' in failed[0].read_text()
  report.append({'memory':bool(memory),'case':name,'status':result.returncode,'gpu_retries':len(events),'failed_attempts':len(attempts),'passed':True})
(root/'result.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'passed':True,'cases':len(report),'output':str(root)}))
