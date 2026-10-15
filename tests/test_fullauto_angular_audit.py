#!/usr/bin/env python3
"""Recover a real independent angular audit without changing its frozen plan."""
import argparse,json,os,signal,subprocess,time
from pathlib import Path
from audit_fullauto_v2 import run as audit


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=True);out=root/'calculation'
    command=[str(a.binary.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1',
        '--ri','1.3116','0','--wavelength-um','1','--max-reflections','4','--threads','1',
        '--fullauto','.02','--fullauto-resources','off','--fullauto-cutoffs','off','--cutoff-profile','off',
        '--fullauto-mirror','off','--symmetry','1','1','--phi-points','4','--scattering-grid','0','25','4','16',
        '--fullauto-workers','persistent','--fullauto-cpu-workers','2','--fullauto-pilot','32',
        '--fullauto-initial','4096','--fullauto-min-correction','256','--fullauto-theta-start','5',
        '--fullauto-orientation-block','128','--max-orientations','65536','--fullauto-audit-angular','on',
        '--output',str(out),'--close']
    (root/'command.args').write_text('\n'.join(command)+'\n')
    with (root/'interrupted.log').open('w') as log:
        process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        deadline=time.monotonic()+180
        while time.monotonic()<deadline:
            if process.poll() is not None:raise RuntimeError('controller finished before audit interruption')
            marker=out/'independent_audit_progress.json'
            if marker.exists() and list((out/'jobs').glob('*/worker.pid')):break
            time.sleep(.01)
        else:raise RuntimeError('no angular audit checkpoint')
        reservation=(out/'independent_audit_reservation.tsv').read_bytes()
        plan=(out/'angular_audit_budget.csv').read_bytes()
        markers={str(f):f.read_bytes() for f in (out/'jobs').glob('*/complete.txt')}
        os.killpg(process.pid,signal.SIGKILL);process.wait()
    with (root/'resumed.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
    assert all(Path(f).read_bytes()==v for f,v in markers.items()),'completed optical blocks changed'
    assert (out/'independent_audit_reservation.tsv').read_bytes()==reservation,'pending independent seeds changed'
    assert (out/'angular_audit_budget.csv').read_bytes()==plan,'holdout changed frozen allocation'
    report=audit(out,None);assert report['independent_final_validation']['estimator']=='angular_residual_v1'
    assert report['independent_final_validation']['passed']
    before=(out/'mueller_fullauto.dat').read_bytes()
    again=subprocess.check_output(command,stderr=subprocess.STDOUT,timeout=30).decode()
    assert 'no calculation repeated' in again and (out/'mueller_fullauto.dat').read_bytes()==before
    report.update(forced_sigkill_recovery=True,completed_markers_preserved=len(markers),frozen_plan_preserved=True,
                  independent_seeds_preserved=True,completed_resume_without_recalculation=True)
    (root/'recovery_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
