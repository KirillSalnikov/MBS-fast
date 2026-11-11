#!/usr/bin/env python3
"""Small real optical tests; deliberately loose tolerance tests stop policy, not physics."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

parser=argparse.ArgumentParser()
parser.add_argument('--binary',required=True,type=Path)
parser.add_argument('--output',required=True,type=Path)
args=parser.parse_args()
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('memory_tests',root/'tests/test_fullauto_memory_runtime.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
s=object.__new__(m.Suite);s.binary=args.binary.resolve();s.args=SimpleNamespace(rounds=3)
output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
receipts={}
for case in ['po_coherent','po_incoherent','go']:
    values={}
    for policy in ['default','off']+(['on'] if case=='go' else []):
        folder=output/(case+'_'+policy)
        cmd=s.command(case,folder,'memory',checkpoint=True)
        pos=cmd.index('--fullauto-final-checks');del cmd[pos:pos+2]
        if policy!='default':cmd+=['--fullauto-final-checks',policy]
        cmd[cmd.index('--fullauto')+1]='.8'
        cmd[cmd.index('--max-orientations')+1]='256'
        (output/(case+'_'+policy+'.args.json')).write_text(json.dumps(cmd,indent=2))
        with (output/(case+'_'+policy+'.log')).open('w') as log:
            process=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=180)
        assert process.returncode==0,(case,policy,process.returncode)
        rules=json.loads((folder/'final_checks_policy.json').read_text())
        status=json.loads((folder/'fullauto_status.json').read_text())
        assert rules['required_production_passes']==1
        assert rules['final_checks']==(policy=='on')
        assert rules['independent_audit_enabled']==(policy=='on')
        assert status['status']=='converged'
        audit=json.loads((folder/'independent_final_validation.json').read_text())
        assert audit['enabled']==(policy=='on')
        if policy=='on':assert audit['passed']
        values[policy]=m.matrix(folder/'mueller_fullauto.dat')
        receipts[case+'_'+policy]={'policy':rules,'status':status['status'],'round':status['round'],'audit':audit}
    assert len(values['default'])==len(values['off']),case
    assert all(abs(a-b)<=1e-11*max(1,abs(a),abs(b)) for x,y in zip(values['default'],values['off']) for a,b in zip(x,y)),case
(output/'receipt.json').write_text(json.dumps(receipts,indent=2))
print('PASS: real PO coherent/incoherent and GO default == explicit off; one successful pass; explicit GO audit on executes and passes')
