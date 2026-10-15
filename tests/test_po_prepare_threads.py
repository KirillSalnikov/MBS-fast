#!/usr/bin/env python3
import argparse,json,subprocess,time,sys
from pathlib import Path
import numpy as np
repo=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description='Bitwise preparation-thread equivalence on convex, asymmetric, absorbing, concave and incoherent PO')
p.add_argument('--binary',type=Path,default=repo/'bin/mbs_po');p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(repo/'tests'));from test_fullauto_arbitrary import asymmetric_tetrahedron
asymmetric_tetrahedron(root/'asymmetric.particle')
base=[str(a.binary.resolve()),'--method','po','--backend','cpu','--threads','1','--sobol-seed','512','11','--haar-alpha','--wavelength-um','1','--ri','1.3116','0','--max-reflections','6','--scattering-grid','0','25','4','8','--orientation-chunk','128','--full-only','--cutoff-profile','off','--symmetry','1','1','--profile-phases','--close']
cases=[('prism',['--particle','1','1','.4','--paired-depth','4']),('asymmetric',['--particle-file',str(root/'asymmetric.particle')]),('absorbing',['--particle-file',str(root/'asymmetric.particle')]),('concave',['--particle-file',str(repo/'examples/particles/concave_hexagonal.particle'),'--resize-dmax-um','.3']),('incoherent',['--particle','1','.2','.1','--incoherent'])]
records=[]
for name,shape in cases:
 matrices=[]
 for threads in [1,4]:
  out=root/'native'/name/str(threads)/'result';out.parent.mkdir(parents=True,exist_ok=True)
  command=base+shape+['--po-prepare-threads',str(threads),'--output',str(out)]
  if name=='absorbing':command[command.index('--ri')+1:command.index('--ri')+3]=['1.3243','.0004']
  (out.parent/(str(threads)+'.args')).write_text('\n'.join(command)+'\n');start=time.monotonic()
  with (out.parent/(str(threads)+'.log')).open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
  matrices.append(((out/'result.dat').read_bytes(),(out/'result_depth4.dat').read_bytes() if name=='prism' else b''));logtext=(out.parent/(str(threads)+'.log')).read_text();assert 'PO preparation threads: '+str(threads) in logtext
  assert 'partial_angular_integral' in logtext and 'WARNING: the C_sca numerical error estimate' not in logtext
  if name=='prism':assert (out/'result_depth4.dat').exists()
  records.append(dict(case=name,threads=threads,seconds=time.monotonic()-start))
 assert matrices[0]==matrices[1],(name,'matrix changed')
 print(name,'all16 bitwise equal',flush=True)
(root/'native_cpu_results.json').write_text(json.dumps(dict(all16_bitwise_equal=True,cases=records),indent=2))
