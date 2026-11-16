#!/usr/bin/env python3
"""CUDA regression: paired-depth pipeline learns RAM before the next enqueue.
Run on an idle CUDA device: python3 tests/test_prepared_pipeline_budget.py BINARY
Creates temporary outputs; requires enough RAM for the 4096 MiB reference run.
"""
import json, math, os, subprocess, sys, tempfile
from pathlib import Path

def matrix(path):
    rows=[]
    for line in path.read_text().splitlines():
        try: row=list(map(float,line.split()))
        except ValueError: continue
        if len(row)==18: rows.append(row)
    assert len(rows)==257 and all(r[2]>0 and all(math.isfinite(v) for v in r) for r in rows)
    return rows

def compare(rows,reference):
    assert all(a[:2]==b[:2] for a,b in zip(rows,reference))
    absolute=max(abs(x-y)/b[2] for a,b in zip(rows,reference) for x,y in zip(a[2:],b[2:]))
    ratios=max(abs(x/a[2]-y/b[2]) for a,b in zip(rows,reference) for x,y in zip(a[2:],b[2:]))
    assert max(absolute,ratios)<3e-5,(absolute,ratios)
    return {'scaled_M11':absolute,'ratios':ratios}

def main():
    binary=str(Path(sys.argv[1]).resolve())
    with tempfile.TemporaryDirectory(prefix='mbs-pipeline-budget-') as tmp:
        root=Path(tmp);grid=root/'theta.csv';grid.write_text('\n'.join(str(i*180/256) for i in range(257))+'\n')
        results=[]
        for budget in (4096,866):
            out=root/str(budget);out.mkdir();target=out/'result'
            cmd=[binary,'--method','po','--backend','cuda','--particle','4','5','--ri','1.3243','0',
                 '--wavelength-um','0.355','--symmetry','2','6','--mirror-gamma','--haar-alpha',
                 '--sobol-seed','64','7001','--max-reflections','18','--paired-depth','16',
                 '--phi-points','8','--theta-grid-file',str(grid),'--threads','5','--po-prepare-threads','5',
                 '--po-trace-backend','cpu','--orientation-pipeline','--orientation-chunk','128',
                 '--cutoff-profile','off','--trace-cache-mib','0','--output',str(target),'--close']
            env=dict(os.environ,MBS_HOST_MEM_BUDGET_MB=str(budget),MBS_GPU_RECYCLE_PREPARED='1',
                     OMP_NUM_THREADS='5',OMP_THREAD_LIMIT='5',OPENBLAS_NUM_THREADS='1')
            log=out/'stdout.log'
            with log.open('w') as f: p=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=120)
            if p.returncode: raise RuntimeError('budget='+str(budget)+'\n'+log.read_text()[-6000:])
            upper=matrix(target/'result.dat');lower=list(out.rglob('*depth16*.dat'));assert len(lower)==1
            results.append((upper,matrix(lower[0])))
        print(json.dumps({'upper':compare(results[1][0],results[0][0]),'lower':compare(results[1][1],results[0][1])}))
if __name__=='__main__':main()
