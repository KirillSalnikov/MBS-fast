#!/usr/bin/env python3
"""Exercise explicit overrides and incoherent physical matrices at runtime."""
import argparse,json,subprocess,time,os
from pathlib import Path
import numpy as np
from audit_fullauto import run as audit
from test_fullauto_arbitrary import rotated_box
ROOT=Path(__file__).resolve().parents[1]

def main():
 p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--backend',default='cpu');a=p.parse_args()
 out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
 particle=out/'rotated_box.particle';rotated_box(particle)
 base=[str(a.binary.resolve()),'--method','po','--backend',a.backend,'--particle-file',str(particle),'--ri','1.31','0','--wavelength-um','1','--max-reflections','2','--fullauto','.01','--symmetry','1','1','--fullauto-mirror','off','--fullauto-resources','auto','--fullauto-pilot','64','--fullauto-initial','256','--fullauto-min-correction','64','--threads','2','--orientation-chunk','32','--fullauto-orientation-block','256','--max-orientations','65536','--scattering-grid','0','180','7','32','--fullauto-theta-start','5','--cutoff-profile','off','--close']
 if a.backend=='cuda':base+=['--fullauto-gpus','0','--orientation-pipeline','--allow-experimental-environment']
 records=[]
 def calculate(name,extra,phi):
  dest=out/name;cmd=base+extra+['--output',str(dest)]
  (out/(name+'.args')).write_text('\n'.join(cmd)+'\n');started=time.monotonic()
  with (out/(name+'.log')).open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=900)
  state=json.loads((dest/'fullauto_status.json').read_text());assert state['status']=='converged'
  assert state['phi']==phi and state['fixed_phi']==phi
  assert state['worker_threads']==2 and state['worker_chunk']==32 and state['orientation_block']==256
  if name=='embedded_phi':assert state['pass_streak']==1
  assert state['target_depth']==2 and state['selected_cutoff']=='requested'
  assert state['rotation_beta_factor']==state['rotation_gamma_factor']==1 and not state['mirror_gamma']
  if a.backend=='cuda':
   assert state['pipeline']
   if os.environ.get('MBS_GPU_BLOCK'):assert state['block']==int(os.environ['MBS_GPU_BLOCK'])
   if os.environ.get('MBS_GPU_WARP_BEAMS'):assert state['warp']==('warp' if os.environ['MBS_GPU_WARP_BEAMS']=='1' else 'thread')
  for path in (dest/'jobs').glob('*/command.args'):
   args=path.read_text().splitlines()
   assert int(args[args.index('--phi-points')+1])==phi
   assert int(args[args.index('--threads')+1])==2
   assert int(args[args.index('--orientation-chunk')+1])==32
   if '--incoherent' in extra:assert '--incoherent' in args
   if '--analytic-mean-cache' in extra and '--analytic-facet-average' in args:
    assert args[args.index('--analytic-mean-cache')+1]==extra[extra.index('--analytic-mean-cache')+1]
   if a.backend=='cuda':assert '--orientation-pipeline' in args
  checked=audit(dest,ROOT/'scripts/gpu_campaign.py',dense=True)
  records.append(dict(case=name,seconds=time.monotonic()-started,state=state,audit=checked));(out/'results.json').write_text(json.dumps(records,indent=2))
  return np.loadtxt(dest/'mueller_fullauto.dat',skiprows=1)
 cache=str(out/'specified_mean.cache')
 calculate('embedded_phi',['--analytic-mean-cache',cache,'--analytic-shadow-control','circular','--checkpoint','--convergence-passes','1'],7)
 assert Path(cache).exists()
 coherent=calculate('coherent_reference',['--phi-points','3','--fullauto-analytic-controls','off'],3)
 incoherent=calculate('incoherent',['--phi-points','3','--incoherent','--fullauto-final-audit','on'],3)
 state=records[-1]['state'];assert not state['coherent'] and not state['analytic_controls']
 sharing=out/'incoherent/depth_sharing_verification.json'
 if sharing.exists():assert not json.loads(sharing.read_text()).get('enabled',False)
 # Incoherent per-beam Mueller summation must change a real scattering
 # matrix; equivalence here would expose a silently ignored mode flag.
 delta=float(np.max(abs(incoherent[:,2:]-coherent[:,2:])/coherent[:,2,None]));assert delta>.01,delta
 records[-1]['coherent_vs_incoherent_difference_scaled_M11']=delta
 (out/'results.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
if __name__=='__main__':main()
