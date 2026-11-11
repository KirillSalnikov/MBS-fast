#!/usr/bin/env python3
import argparse,importlib.util,json,subprocess,time
from pathlib import Path
from types import SimpleNamespace
parser=argparse.ArgumentParser();parser.add_argument('--binary',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
r=Path(__file__).resolve().parents[1];spec=importlib.util.spec_from_file_location('memory_tests',r/'tests/test_fullauto_memory_runtime.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
s=object.__new__(m.Suite);s.binary=args.binary.resolve();s.args=SimpleNamespace(rounds=8)
out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);receipt={}
for count in [256,32]:
 p=out/('pilot'+str(count));cmd=s.command('po_coherent',p,'memory',checkpoint=True)
 cmd[cmd.index('--particle')+1]='2'
 cmd[cmd.index('--max-orientations')+1]='4096'
 pos=cmd.index('--scattering-grid');cmd[pos+1:pos+5]=['0','180','4','128']
 cmd[cmd.index('--fullauto-theta-audit-count')+1]=str(count)
 (out/('command'+str(count)+'.json')).write_text(json.dumps(cmd,indent=2))
 start=time.monotonic()
 with (out/('pilot'+str(count)+'.log')).open('w') as f:job=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=450)
 seconds=time.monotonic()-start
 state=json.loads((p/'fullauto_status.json').read_text());pilot=json.loads((p/'dense_theta_pilot.json').read_text())
 receipt[str(count)]={'seconds':seconds,'exit':job.returncode,'status':state,'pilot':pilot}
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2));print(count,seconds,job.returncode,state['status'],flush=True)
assert all(x['exit']==0 for x in receipt.values()),'Full convergence not established'
a=m.matrix(out/'pilot256/mueller_fullauto.dat');b=m.matrix(out/'pilot32/mueller_fullauto.dat');assert len(a)==len(b)
maximum=max(abs(x[j]-y[j])/max(x[2],1e-300) for x,y in zip(a,b) for j in range(2,18))
receipt['comparison']={'maximum_all16_difference_scaled_M11':maximum,'speedup':receipt['256']['seconds']/receipt['32']['seconds']}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2));print(receipt['comparison'])
