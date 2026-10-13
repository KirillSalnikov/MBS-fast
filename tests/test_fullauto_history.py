#!/usr/bin/env python3
"""Check that real adaptive decisions and all per-angle budgets survive logging."""
import argparse,json,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--epsilon',type=float,default=.01);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
 grid=out/'theta.grid';grid.write_text('\n'.join(str(v) for v in np.linspace(0,25,17))+'\n')
 dest=out/'calculation';command=[str(a.binary.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1','--ri','1.3116','0','--wavelength-um','1','--max-reflections','2','--threads','2','--fullauto',str(a.epsilon),'--fullauto-resources','off','--fullauto-cutoffs','off','--cutoff-profile','off','--symmetry','1','1','--fullauto-mirror','off','--fullauto-work-queue','on','--fullauto-orientation-block','128','--fullauto-pilot','32','--fullauto-initial','512','--fullauto-min-correction','64','--fullauto-theta-start','5','--max-orientations','32768','--theta-grid-file',str(grid),'--output',str(dest)]
 (out/'command.args').write_text('\n'.join(command)+'\n')
 with (out/'calculation.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=600)
 entries=[json.loads(line) for line in (dest/'fullauto_history.jsonl').read_text().splitlines()]
 types={e['event'] for e in entries};assert {'worker_start','worker_complete','allocation_checkpoint','round_decision','phi_candidate','round_diagnostics'}<=types,types
 final=json.loads((dest/'final_parameters.json').read_text());assert final['physical_target_depth']==2 and final['worker_threads']==2 and final['orientation_block']==128
 assert final['selected_state']['status']=='converged' and abs(final['interpolation_budget']-a.epsilon/4)<1e-14
 candidates={int(e['data']['phi']) for e in entries if e['event']=='phi_candidate'};assert candidates=={1,2,4,8,16,32,64},candidates
 rounds=[e for e in entries if e['event']=='round_diagnostics'];assert rounds
 for e in rounds:
  folder=Path(e['data']['folder']);table=np.genfromtxt(folder/'parameters_accuracy.csv',delimiter=',',names=True)
  matrix=np.loadtxt(folder/'mueller_evaluated.dat',skiprows=1,ndmin=2)
  assert np.allclose(table['theta_deg'],matrix[:,0]);assert all(name in table.dtype.names for name in ['N_level0','N_level1','estimated95_M12_over_M11','interpolation_residual95'])
 before=sum(e['event']=='worker_start' for e in entries)
 with (out/'resume.log').open('w') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=300)
 after=[json.loads(line) for line in (dest/'fullauto_history.jsonl').read_text().splitlines()];assert sum(e['event']=='worker_start' for e in after)==before
 result=dict(events_before_resume=len(entries),events_after_resume=len(after),phi_candidates=sorted(candidates),round_snapshots=len(rounds),worker_starts=before,no_worker_repeated_on_resume=True,final_parameters_valid=True,all_jsonl_lines_valid=True)
 (out/'history_test.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
