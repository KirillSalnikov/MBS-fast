#!/usr/bin/env python3
"""Prepare a portable, physically matched eight-replica SO(3) benchmark."""
import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tests'))
from test_fullauto_arbitrary import asymmetric_tetrahedron
from generate_so3_rules import make_rule,write_rule


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--degrees',type=int,nargs='+',default=[15,31])
    p.add_argument('--reference',type=int,default=65536)
    p.add_argument('--reference-seed-offset',type=int,default=0)
    p.add_argument('--cases',nargs='+',choices=['prism_1um','asymmetric_1um','prism_H501','prism_H1000'])
    a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=True)
    (out/'inputs').mkdir(exist_ok=True)
    asymmetric_tetrahedron(out/'inputs/asymmetric_tetrahedron.particle')
    grids={'full':[0,.1,1,5,10,20,25,90,170,175,178,179,180],
           'forward':[0,.1,1,5,10,15,20,21.895,22,24.989,25]}
    for key,values in grids.items():
        (out/('inputs/theta_'+key+'.grid')).write_text('\n'.join(map(str,values))+'\n')
    rules={};moments=[]
    for method in ['lebedev','wigner']:
        for degree in a.degrees:
            angles,weights,error=make_rule(method,degree)
            folder=out/'rules'/(method+'_'+str(degree));folder.mkdir(parents=True,exist_ok=True)
            write_rule(folder/'unrotated.rule',angles,weights)
            for replica in range(8):
                path=folder/('replica_'+str(replica)+'.rule')
                write_rule(path,angles,weights,22001+replica*101)
            rules[(method,degree)]=len(weights)
            moments.append(dict(method=method,degree=degree,orientations=len(weights),checked_moment_error=error))
    cases=[('prism_1um',['--particle','1','1','1'],'full','6'),
           ('asymmetric_1um',['--particle-file','inputs/asymmetric_tetrahedron.particle'],'full','6'),
           ('prism_H501',['--particle','1','501.2','155.8'],'forward','12'),
           ('prism_H1000',['--particle','1','1000','220.1'],'forward','12')]
    jobs=[];seeds=[11,23,37,53,71,89,107,131]
    for case,shape,grid,depth in cases:
        if a.cases and case not in a.cases:continue
        base=shape+['--symmetry','1','1','--ri','1.3116','0','--wavelength-um','.532',
                    '--max-reflections',depth,'--beam-cutoff','.001','--phi-points','4',
                    '--theta-grid-file','inputs/theta_'+grid+'.grid',
                    '--analytic-facet-average','--analytic-shadow-control','circular',
                    '--analytic-mean-cache','cache_'+case,'--orientation-pipeline']
        # Interleave methods and replicas to avoid systematic timing drift.
        configs=[]
        for degree in a.degrees:
            sobol=1<<(rules[('wigner',degree)]-1).bit_length()
            configs.extend([('sobol',str(degree),sobol),
                            ('lebedev',str(degree),rules[('lebedev',degree)]),
                            ('wigner',str(degree),rules[('wigner',degree)])])
        configs.append(('sobol','reference',a.reference))
        for method,level,count in configs:
            for replica in range(8):
                job=dict(id='{}_{}_{}_r{}'.format(case,method,level,replica),case=case,method=method,
                         level=level,orientations=count,replica=replica,
                         seed=seeds[replica]+(a.reference_seed_offset if level=='reference' else 0),
                         args=base,theta_rows=len(grids[grid]))
                if method!='sobol':job['rule']='rules/{}_{}/replica_{}.rule'.format(method,level,replica)
                jobs.append(job)
    manifest=dict(jobs=jobs,rules=moments,created_by=str(Path(__file__).name),
                  rotations='ZYZ: Rz(alpha) Ry(beta) Rz(gamma); Lebedev sphere x alpha Fourier; Wigner Gauss(beta) x Fourier(alpha,gamma), integral zero mode.',
                  comparison='Common physical parameters, full SO(3), phi=4, identical direct theta for each case. Eight independent Owen scrambles or global Haar rotations. No symmetry folding.',
                  reference='Finite-N scrambled Sobol; report its dispersion, do not treat it as exact.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(jobs=len(jobs),rules=moments),indent=2))


if __name__=='__main__':main()
