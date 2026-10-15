#!/usr/bin/env python3
"""Real optical holdouts: local vs uniform, NumPy reconstruction and SIGKILL resume."""
import argparse,json,os,signal,subprocess,time
from pathlib import Path
import numpy as np
from audit_fullauto_v2 import matrix,grid,reconstruct,intervals,local_spline,common_bounds,interpolation_bounds
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe',type=Path,required=True);p.add_argument('--worker',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=True)
 base=[str(a.probe.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1','--ri','1.3116','0','--wavelength-um','1','--max-reflections','4','--threads','1','--fullauto','.1','--fullauto-analytic-controls','off','--fullauto-resources','off','--fullauto-cutoffs','off','--fullauto-rotation-symmetry','off','--fullauto-mirror','off','--phi-points','4','--scattering-grid','0','25','4','64','--fullauto-initial','512','--fullauto-min-correction','128','--fullauto-orientation-block','128','--max-orientations','65536','--close']
 records=[]
 def run(name,policy,interrupt=False):
  out=root/name;cmd=base+['--output',str(out)];env=dict(os.environ,TASK_MBS_WORKER=str(a.worker.resolve()),TASK_MBS_AUDIT_POLICY=policy)
  (root/(name+'.args')).write_text('\n'.join(cmd)+'\n');started=time.monotonic();saved={}
  with (root/(name+'.log')).open('w') as log:
   if interrupt:
    process=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    deadline=time.monotonic()+300
    while time.monotonic()<deadline:
     assert process.poll() is None,'probe ended before interruption'
     if (out/'independent_audit_plan.json').exists():
      saved={str(p):p.stat().st_mtime_ns for p in (out/'jobs').glob('j*/complete.txt')
             if int((p.parent/'sampling.index').read_text().splitlines()[0].split()[4])>=1000000000}
      if saved:break
     time.sleep(.03)
    else:raise AssertionError('no committed holdout before timeout')
    plan=json.loads((out/'independent_audit_plan.json').read_text());os.killpg(process.pid,signal.SIGKILL);process.wait();time.sleep(.15)
    subprocess.run(cmd,env=dict(env,TASK_MBS_WARM_CACHE='1'),stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    assert all(Path(p).stat().st_mtime_ns==stamp for p,stamp in saved.items()),'finished holdout recomputed'
    resumed=json.loads((out/'independent_audit_plan.json').read_text());assert resumed['resumed_pending_audit'] and resumed['audit_serial']==plan['audit_serial']
   else:subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
  state=json.loads((out/'fullauto_status.json').read_text());report=json.loads((out/'independent_final_validation.json').read_text());assert report['passed']
  assert set(report['seeds']).isdisjoint(state['seeds'])
  angles=grid(out/'independent_theta.csv');support=grid(out/'theta_support.csv');budget=np.loadtxt(out/'independent_orientation_budget.csv',delimiter=',',skiprows=1,ndmin=2)
  assert len(angles)==65 and np.allclose(angles,budget[:,0]);ns=budget[:,1:].astype(int)
  levels,seen,_=reconstruct(out,state,angles,ns,out/'independent_final_sources.tsv',report['seeds']);samples=levels.sum(0)
  mean,raw,ratios=intervals(samples);spline=local_spline(support,samples[:,np.searchsorted(angles,support)]);_,out_raw,out_ratios=intervals(spline(grid(out/'requested_theta.csv')))
  raw95=max(raw.max(),out_raw.max());ratio95=max(ratios.max(),out_ratios.max())
  assert abs(raw95-report['max_pointwise95_scaled_M11'])<1e-11 and abs(ratio95-report['max_ratio95'])<1e-11
  mode=json.loads((out/'interpolation_validation_mode.json').read_text())
  if mode['common_prefix']:bound,checked=common_bounds(out,state,'validation',angles,support,report['seeds'],ns);seen.update(checked)
  else:bound=interpolation_bounds(spline(angles),samples,True)
  assert abs(bound.max()-report['max_interpolation_residual95'])<1e-10
  assert max(raw95,ratio95)<=state['statistical_budget'] and bound.max()<=state['interpolation_budget']
  if policy=='local':assert len(np.unique(ns[:,0]))>1 and mode['common_prefix']
  records.append(dict(case=name,policy=policy,seconds=time.monotonic()-started,native_audit_seconds=json.loads((out/'audit_seconds.json').read_text())['seconds'],all16_reconstructed=True,confidence95_reconstructed=True,common_prefix_reconstructed=mode['common_prefix'],maximum95=float(max(raw95,ratio95)),interpolation95=float(bound.max()),checked_angles=len(angles),preserved_holdout_blocks=len(saved),jobs_validated=len(seen)))
  (root/'results.json').write_text(json.dumps(records,indent=2))
  return out,samples
 local,ls=run('local','local');uniform,us=run('uniform','uniform')
 assert (local/'frozen_production.dat').read_bytes()==(uniform/'frozen_production.dat').read_bytes(),'production estimator changed'
 lm,lc,lr=intervals(ls);um,uc,ur=intervals(us);difference=np.max(abs(lm-um)/lm[:,0,None]);assert difference<lc.max()+uc.max()
 auto,_=run('auto','auto')
 assert not json.loads((auto/'independent_audit_cost_gate.json').read_text())['selected_local']
 assert not json.loads((auto/'independent_final_validation.json').read_text())['local_orientation_budgets']
 run('interrupted','local',True)
 (root/'results.json').write_text(json.dumps(dict(cases=records,local_uniform_difference_scaled_M11=float(difference),same_frozen_production=True),indent=2));print((root/'results.json').read_text())
if __name__=='__main__':main()
