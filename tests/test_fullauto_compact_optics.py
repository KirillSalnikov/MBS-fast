import sys,importlib.util,json,subprocess,concurrent.futures,argparse
from pathlib import Path
from types import SimpleNamespace
R=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--output',required=True,type=Path)
parser.add_argument('--binary',required=True,type=Path)
args=parser.parse_args()
spec=importlib.util.spec_from_file_location('reg',R/'tests/test_fullauto_memory_runtime.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
O=args.output.resolve();O.mkdir(exist_ok=False,parents=True)
s=object.__new__(m.Suite);s.binary=args.binary.resolve();s.args=SimpleNamespace(rounds=2)
def run(case):
 results={}
 for layout in ['legacy','compact']:
  out=O/(case+'_'+layout);assert not out.exists()
  cmd=s.command(case,out,'memory',checkpoint=True)
  if layout=='compact':cmd+=['--fullauto-layout','compact']
  (O/(case+'_'+layout+'.args.json')).write_text(json.dumps(cmd,indent=2))
  with (O/(case+'_'+layout+'.log')).open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=180)
  assert p.returncode in (0,3),(case,layout,p.returncode)
  work=out/'work' if layout=='compact' else out
  state=json.loads((work/'fullauto_status.json').read_text())
  results[layout]=dict(code=p.returncode,status=state['status'],round=state['round'])
  assert (work/'fullauto_checkpoint.bin').exists()
  if layout=='compact':
   latest=(out/'latest').resolve(strict=True)
   assert (latest/'mueller_bins.dat').exists()
   accuracy=json.loads((latest/'accuracy.json').read_text());assert accuracy['method']==('go' if case=='go' else 'po')
   assert not (work/'public_snapshots').exists()
   before=m.matrix(latest/'mueller_bins.dat')
   with (O/(case+'_resume.log')).open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,timeout=180)
   assert p.returncode in (0,3),(case,'resume',p.returncode)
   after=m.matrix(out/'latest/mueller_bins.dat')
   assert len(after)==len(before)
   assert all(abs(v-w)<=1e-11*max(1,abs(w)) for row,old in zip(after,before) for v,w in zip(row,old))
  else:
   pointer=json.loads((out/'intermediate_current.json').read_text());source=Path(pointer['snapshot_folder'])
   ref=(source/'mueller_evaluated.dat').read_bytes()
   if case=='go':reference=ref
   else:
    # Full production table includes interpolation; compare all 16 columns against legacy annotated output.
    reference=list(__import__('csv').DictReader((source/'mueller_fullauto_annotated.csv').open()))
 if case=='go':assert (O/(case+'_compact/latest/mueller_bins.dat')).read_bytes()==reference
 else:
  rows=m.matrix(O/(case+'_compact/latest/mueller_bins.dat'));assert len(rows)==len(reference)
  for row,old in zip(rows,reference):
   for j in range(16):assert abs(row[j+2]-float(old['M'+str(j//4+1)+str(j%4+1)]))<=1e-11*max(1,abs(row[j+2]))
 return case,results
with concurrent.futures.ThreadPoolExecutor(2) as ex:
 for case,result in ex.map(run,['po_coherent','go','po_incoherent']):
  print(case,result,flush=True);(O/(case+'.receipt.json')).write_text(json.dumps(result,indent=2))
