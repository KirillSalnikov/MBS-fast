#!/usr/bin/env python3
"""Separate analytic CPU pool: real prism, exact matrices, cache and CLI checks."""
import argparse,json,subprocess,time,struct
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--setup-threads',type=int,default=4);p.add_argument('--theta-points',type=int,default=9);a=p.parse_args()
root=a.output.resolve();root.mkdir(parents=True,exist_ok=True);binary=str(a.binary.resolve())
# Actual current forward particle and sparse points from its requested range.
grid=root/'theta.csv';angles=[0,.004,1,3,6,10,15,20,25] if a.theta_points==9 else [25.*i/(a.theta_points-1) for i in range(a.theta_points)];grid.write_text('\n'.join(format(x,'.17g') for x in angles)+'\n')
def values(path):
 return b''.join(struct.pack('!d',float(x)) for line in path.read_text().splitlines()[1:] for x in line.split())
base=[binary,'--method','po','--backend','cpu','--particle','1','446.7','147.1','--ri','1.3116','0','--wavelength-um','.532','--max-reflections','2','--sobol-seed','4','71','--haar-alpha','--phi-points','4','--theta-grid-file',str(grid),'--symmetry','1','1','--threads','1','--analytic-facet-average','--analytic-shadow-control','circular','--cutoff-profile','off','--close']
records=[];matrices=[];means=[]
for name,threads,cache in [('serial',1,root/'serial.cache'),('parallel',a.setup_threads,root/'parallel.cache'),('warm',a.setup_threads,root/'serial.cache')]:
 dest=root/name;cmd=base+['--analytic-mean-threads',str(threads),'--analytic-mean-cache',str(cache),'--output',str(dest)]
 begin=time.monotonic();run=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,timeout=600);seconds=time.monotonic()-begin
 (root/(name+'.log')).write_text(run.stdout);assert run.returncode==0,run.stdout[-4000:]
 assert 'OpenMP threads: 1' in run.stdout
 if name!='warm':assert 'Analytic mean setup: rows='+str(len(angles))+' threads='+str(min(threads,len(angles))) in run.stdout
 else:assert 'missing=0' in run.stdout and 'Analytic mean setup:' not in run.stdout
 matrices.append(values(dest/(name+'.dat')));means.append(values(cache));records.append({'case':name,'seconds':seconds,'command':cmd})
assert matrices[0]==matrices[1]==matrices[2];assert means[0]==means[1]==means[2]
for invalid in ['0','-1']:
 run=subprocess.run(base+['--analytic-mean-threads',invalid,'--output',str(root/'invalid')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
 assert run.returncode!=0 and '--analytic-mean-threads must be positive' in run.stdout
(root/'comparison.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps({'bitwise_equal':True,'speedup':records[0]['seconds']/records[1]['seconds'],'runs':records},indent=2))
