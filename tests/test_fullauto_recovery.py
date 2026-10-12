#!/usr/bin/env python3
"""Kill native controller mid-block, verify durable exact recovery."""
import argparse,json,os,signal,subprocess,time
from pathlib import Path
import numpy as np
from audit_fullauto import run as audit
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    base=[str(a.binary.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1','--ri','1.3116','0','--wavelength-um','1','--max-reflections','4','--threads','2','--fullauto','.01','--fullauto-mirror','off','--phi-points','4','--scattering-grid','0','25','4','32','--fullauto-pilot','64','--fullauto-initial','4096','--fullauto-min-correction','256','--fullauto-theta-start','5','--fullauto-orientation-block','256','--max-orientations','131072','--cutoff-profile','off']
    target=out/'interrupted';command=base+['--output',str(target)]
    with (out/'interrupted.log').open('w') as log:
        process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        deadline=time.monotonic()+180
        while time.monotonic()<deadline:
            if process.poll() is not None:raise RuntimeError('Controller ended before interruption')
            if (target/'orientation_progress.json').exists():break
            time.sleep(.05)
        else:raise RuntimeError('No durable block checkpoint')
        markers={str(p):p.stat().st_mtime_ns for p in (target/'jobs').glob('*/complete.txt')}
        progress=json.loads((target/'orientation_progress.json').read_text())
        os.killpg(process.pid,signal.SIGKILL);process.wait()
    time.sleep(.2)
    with (out/'resume.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)
    assert all(Path(p).stat().st_mtime_ns==mtime for p,mtime in markers.items()),'Completed work was recomputed'
    report=audit(target,ROOT/'scripts/gpu_campaign.py',dense=True)
    history=np.loadtxt(target/'refinements.tsv',ndmin=2);assert len(np.unique(history[:,0]))==len(history),'Duplicated rounds'
    marker=target/'mueller_fullauto.dat';original=marker.read_bytes();marker.write_bytes(original+b'corrupted\n')
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    marker.write_bytes(original);assert result.returncode!=0 and 'completed output changed' in result.stdout
    result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,check=True)
    assert 'no calculation repeated' in result.stdout
    report.update(interrupted_prefix=progress,completed_jobs_preserved=len(markers),forced_sigkill_recovery=True,corrupt_final_output_rejected=True)
    (out/'recovery.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
