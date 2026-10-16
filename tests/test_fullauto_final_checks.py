#!/usr/bin/env python3
"""Skip both final audits on a real unfinished automatic-depth PO checkpoint."""
import argparse, json, os, signal, subprocess, time
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=True);out=root/'calculation'
    command=[str(a.binary.resolve()),'--method','po','--backend','cpu','--particle','1','.2','.1',
        '--ri','1.3116','0','--wavelength-um','1','--incoherent','--threads','2',
        '--fullauto','.02','--fullauto-resources','off','--fullauto-cutoffs','off','--cutoff-profile','off',
        '--fullauto-mirror','off','--fullauto-rotation-symmetry','off','--symmetry','1','1',
        '--fullauto-analytic-controls','off','--phi-points','4','--fullauto-pilot','32',
        '--fullauto-depth-pilot','64','--fullauto-initial','512','--fullauto-min-correction','64',
        '--fullauto-theta-start','3','--max-orientations','32768','--scattering-grid','0','25','4','4',
        '--fullauto-work-queue','on','--fullauto-orientation-block','256','--output',str(out),'--close']
    (root/'command.args').write_text('\n'.join(command)+'\n')
    with (root/'interrupted.log').open('w') as log:
        process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            deadline=time.monotonic()+240
            while time.monotonic()<deadline:
                if process.poll() is not None:raise RuntimeError('controller finished before audit interruption')
                marker=out/'independent_audit_plan.json'
                if marker.exists():break
                time.sleep(.005)
            else:raise RuntimeError('no independent audit checkpoint')
        finally:
            if process.poll() is None:os.killpg(process.pid,signal.SIGKILL);process.wait()
    config=(out/'configuration.txt').read_bytes();depth=(out/'automatic_depth.tsv').read_bytes()
    markers={str(f):f.read_bytes() for f in (out/'jobs').glob('*/complete.txt')}
    prior=json.loads((out/'fullauto_status.json').read_text());assert prior['automatic_depth']
    resumed=command+['--fullauto-final-checks','off']
    with (root/'resumed.log').open('w') as log:
        subprocess.run(resumed,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=240)
    status=json.loads((out/'fullauto_status.json').read_text())
    assert status['status']=='converged' and status['automatic_depth'] and not status['final_checks']
    assert not status['independent_final_validation'] and not status['final_physics_audit']
    assert status['combined_error95']<=.02 and status['max_guard_residual95_scaled_M11']<=.005
    for name in ['independent_final_validation.json','physics_validation.json']:
        report=json.loads((out/name).read_text());assert not report['enabled'] and not report['passed']
    contract=json.loads((out/'precision_contract.json').read_text())
    assert not contract['physical_accuracy_audited'] and not contract['independent_audit_enabled']
    assert (out/'configuration.txt').read_bytes()==config and (out/'automatic_depth.tsv').read_bytes()==depth
    assert all(Path(f).read_bytes()==v for f,v in markers.items())
    events=[json.loads(x) for x in (out/'fullauto_history.jsonl').read_text().splitlines()]
    assert not any(e['event']=='worker_start' and e['data']['seed']>=1_000_000_000 for e in events if e['controller_pid']==status['controller_pid'])
    matrix=(out/'mueller_fullauto.dat').read_bytes();before=len([x for x in events if x['event']=='worker_start'])
    again=subprocess.check_output(resumed,stderr=subprocess.STDOUT,timeout=60).decode()
    assert 'no calculation repeated' in again and (out/'mueller_fullauto.dat').read_bytes()==matrix
    strict=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=30)
    assert strict.returncode!=0 and b'completed final-checks policy differs' in strict.stdout
    assert not json.loads((out/'final_checks_policy.json').read_text())['final_checks']
    report=dict(status='passed',production_error=status['combined_error95'],automatic_depth_retained=True,
                configuration_preserved=True,committed_receipts_preserved=len(markers),final_audits_disabled=True,
                completed_resume_without_recalculation=True,completed_policy_upgrade_rejected=True)
    (root/'verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report))


if __name__=='__main__':main()
