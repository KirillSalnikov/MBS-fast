#!/usr/bin/env python3
"""Compare all Mueller entries with a frozen multilevel/per-angle allocation."""
import argparse,json,os,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,default=ROOT/'bin/mbs_po')
    p.add_argument('--probe',type=Path,default=ROOT/'bin/fullauto_queue_probe')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    env=dict(os.environ,TASK_MBS_WORKER=str(a.binary.resolve()))
    base=[str(a.probe.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1','--ri','1.3116','0','--wavelength-um','1','--max-reflections','4','--threads','2','--fullauto','.01','--fullauto-analytic-controls','off','--fullauto-mirror','off','--fullauto-resources','off','--cutoff-profile','off','--phi-points','4','--scattering-grid','0','180','4','24','--fullauto-initial','2048','--fullauto-orientation-block','128','--fullauto-depth-sharing','off']
    reports=[]
    for mode in ['coherent','incoherent','paired','controls']:
        matrices=[];timings={}
        for queue in ['off','on']:
            dest=out/(mode+'_'+queue);command=base+['--fullauto-work-queue',queue,'--output',str(dest)]
            if mode=='incoherent':command+=['--incoherent']
            if mode=='paired':command[command.index('--fullauto-depth-sharing')+1]='on'
            if mode=='controls':command[command.index('--fullauto-analytic-controls')+1]='on'
            (out/(dest.name+'.args')).write_text('\n'.join(command)+'\n')
            with (out/(dest.name+'.log')).open('w') as log:subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
            matrices.append(np.loadtxt(dest/'frozen_matrix.dat',skiprows=1))
            timings[queue]=float((dest/'frozen_seconds.txt').read_text())
        error=float(np.max(abs(matrices[0][:,2:]-matrices[1][:,2:])/matrices[0][:,2,None]))
        assert error<1e-10,(mode,error)
        reports.append(dict(mode=mode,all16_difference_scaled_M11=error,seconds=timings))
    probe=subprocess.run([str(a.probe.resolve()),'--test-file-hash',str(out/'hash_fixture.dat')],stdout=subprocess.PIPE,universal_newlines=True,check=True)
    result=dict(matrices=reports,file_hash=json.loads(probe.stdout))
    (out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
