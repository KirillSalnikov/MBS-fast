#!/usr/bin/env python3
"""Run with the test_particle_planarity_rounding C++ loader executable."""
import math,subprocess,sys,tempfile
from pathlib import Path

def geometry(scale,warp=0):
    e=.05
    v=[(10,0,0),(0,10,0),(-10,0,0),(0,-10,0),(0,0,-10),
       (e,0,10-e+warp),(0,e,10-e),(-e,0,10-e),(0,-e,10-e)]
    faces=[[0,1,6,5],[1,2,7,6],[2,3,8,7],[3,0,5,8],
           [1,0,4],[2,1,4],[3,2,4],[0,3,4],[5,6,7,8]]
    def rotate(p):
        x,y,z=p
        for angle,axis in [(.713,'x'),(.431,'y'),(.219,'z')]:
            c,s=math.cos(angle),math.sin(angle)
            if axis=='x':y,z=c*y-s*z,s*y+c*z
            elif axis=='y':x,z=c*x+s*z,-s*x+c*z
            else:x,y=c*x-s*y,s*x+c*y
        return x*scale,y*scale,z*scale
    v=list(map(rotate,v))
    return '0\n0\n180 360\n\n'+'\n\n'.join('\n'.join(' '.join(format(a,'.17g') for a in v[k]) for k in face) for face in faces)+'\n'
with tempfile.TemporaryDirectory() as directory:
    p=Path(directory)/'particle.dat'
    for scale in [.001,1,1000]:
        p.write_text(geometry(scale));r=subprocess.run([sys.argv[1],str(p)],capture_output=True,text=True)
        assert r.returncode==0,(scale,r.stderr)
        p.write_text(geometry(scale,1e-5));r=subprocess.run([sys.argv[1],str(p)],capture_output=True,text=True)
        assert r.returncode!=0 and 'not planar' in r.stderr,(scale,r.stdout,r.stderr)
print('PASS planar small facets at three scales; real nonplanarity rejected')
