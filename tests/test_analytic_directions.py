#!/usr/bin/env python3
"""Independent aperture references for arbitrary incidence and rotation on CPU/CUDA."""
import argparse
from io import StringIO
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.special import j1
from scipy.spatial.transform import Rotation

from test_analytic_facets import ROOT, R, polygon, transform


def run(probe, arguments, rows):
    stream=''.join(' '.join(map(str,row))+'\n' for row in rows)
    output=subprocess.check_output([str(probe),*map(str,arguments)],input=stream,text=True)
    return np.atleast_2d(np.loadtxt(StringIO(output)))


def check(probe, gpu=False):
    rng=np.random.default_rng(41387)
    q=Rotation.random(random_state=rng).as_quat()
    k=2*np.pi/.532
    worst=0.
    # Incidence on the facet spans normal, grazing, and unilluminated cases.
    mu=np.r_[1.,1-1e-12,1e-8,0.,-1e-8,-.3,rng.uniform(-1,1,90)]
    az=rng.uniform(-np.pi,np.pi,len(mu))
    u=np.column_stack((np.sqrt(1-mu*mu)*np.cos(az),np.sqrt(1-mu*mu)*np.sin(az),mu))
    reference=np.column_stack((-np.sin(az),np.cos(az),np.zeros(len(mu))))
    for degrees in [0.,1e-7,17.,45.,90.,135.,170.,179.9999999,180.]:
        theta=np.deg2rad(degrees)
        phi=rng.uniform(-np.pi,np.pi,len(mu))
        frame_rows=np.column_stack((u,reference,np.full(len(mu),theta),phi))
        frames=run(probe,['frame'],frame_rows)
        axes=frames[:,:9].reshape(-1,3,3)
        v=frames[:,9:]
        assert np.max(abs(axes@axes.transpose(0,2,1)-np.eye(3)))<2e-15
        assert np.max(abs(np.linalg.det(axes)-1))<2e-15
        assert np.max(abs(np.sum(u*v,axis=1)+np.cos(theta)))<1e-15
        # A simultaneous rotation preserves the azimuth reference convention.
        rotation=Rotation.from_quat(q).as_matrix()
        rotated=run(probe,['frame'],np.column_stack((u@rotation.T,reference@rotation.T,
            np.full(len(mu),theta),phi)))
        assert np.max(abs(rotated[:,9:]-v@rotation.T))<1e-15
        for shape,scale,index in [('rectangle',1.,1.31+0j),('triangle',3.,1.53+.0018j),('hex',1.,1.31+0j)]:
            p=polygon(shape,scale)
            area=abs(np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1)))/2
            for shadow in ['facets','circular','off']:
                arguments=['gpu-directions' if gpu else 'directions',shape,scale,
                    index.real,index.imag,shadow,degrees,*q]
                if gpu:
                    # phi=0 in the CUDA grid; encode the desired azimuth in
                    # each physical source frame's transverse reference.
                    transverse=(np.cos(phi)[:,None]*axes[:,0]+np.sin(phi)[:,None]*axes[:,1])
                    rows=np.column_stack((u,transverse))
                else:
                    # Unequal vector lengths must not alter directional physics.
                    rows=np.column_stack((u*rng.uniform(.2,4,(len(mu),1)),v*rng.uniform(.2,4,(len(mu),1))))
                actual=run(probe,arguments,rows)
                h,c=np.sin(theta/2),np.cos(theta/2)
                if degrees==180:c=0.
                if degrees==0:h=0.
                illuminated=mu>0
                nu=v[:,2]
                amplitude=transform(p,k*(u[:,0]+v[:,0]),k*(u[:,1]+v[:,1]))
                reflected=illuminated*h*h*np.maximum(0,c*c+mu*nu)*R(np.maximum(0,mu),index)*abs(amplitude)**2/.532**2
                if shadow=='off':sh=np.zeros(len(mu))
                elif shadow=='circular':
                    x=k*np.sqrt(area/(4*np.pi))*np.sin(theta)
                    envelope=1. if abs(x)<1e-12 else 2*j1(x)/x
                    sh=illuminated*c**4*envelope**2*mu**2*area**2/.532**2
                else:
                    tangent=v+np.cos(theta)*u
                    amplitude=transform(p,k*tangent[:,0],k*tangent[:,1])
                    sh=illuminated*c**4*mu**2*abs(amplitude)**2/.532**2
                expected=np.column_stack((reflected,sh))
                # Normalize by the finite aperture intensity scale: a relative
                # error at a diffraction zero would be meaningless.
                scale_intensity=area**2/.532**2
                error=np.max(abs(actual[:,:2]-expected))/scale_intensity
                assert error<3e-11,(degrees,shape,shadow,error)
                for first in range(2,actual.shape[1],2):
                    covariance=np.max(abs(actual[:,first:first+2]-actual[:,:2]))/scale_intensity
                    assert covariance<3e-11,(degrees,shape,shadow,first,covariance,
                        actual[np.argmax(np.max(abs(actual[:,first:first+2]-actual[:,:2]),axis=1))])
                    worst=max(worst,covariance)
                worst=max(worst,error)
    # Invalid directions/azimuth references fail, rather than silently changing
    # the incident polarization plane.
    reference=run(probe,['frame'],[[.6,0,.8,0,1,0,.7,1.3]])
    for magnitude in [1e-308,1e308]:
        extreme=run(probe,['frame'],[[.6*magnitude,0,.8*magnitude,0,1,0,.7,1.3]])
        assert np.max(abs(extreme-reference))<1e-15
    for row in [[0,0,0,1,0,0,0,0],[0,0,1,0,0,2,0,0],[0,0,1,1,0,0,-.1,0]]:
        result=subprocess.run([str(probe),'frame'],input=' '.join(map(str,row))+'\n',
            text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        assert result.returncode!=0
    print(f'PASS: {"CUDA/CPU" if gpu else "CPU"} arbitrary incidence, independent area integrals, '
          f'rotation covariance, poles and invalid inputs; max aperture-scaled error {worst:.3g}')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--probe',type=Path)
    parser.add_argument('--gpu',action='store_true')
    args=parser.parse_args()
    if args.probe:
        check(args.probe.resolve(),args.gpu)
        return
    if args.gpu:parser.error('--gpu requires a CUDA-built --probe')
    with tempfile.TemporaryDirectory(prefix='mbs_analytic_directions_') as directory:
        probe=Path(directory)/'probe'
        subprocess.run(['g++','-std=gnu++11','-O2','-march=native','-Wall','-Wextra',
            '-I'+str(ROOT/'src'),'-I'+str(ROOT/'src/math'),'-I'+str(ROOT/'src/handler'),
            str(ROOT/'tests/analytic_facet_probe.cpp'),str(ROOT/'src/AnalyticFacetAverage.cpp'),
            str(ROOT/'src/AnalyticBackscatter.cpp'),str(ROOT/'src/math/Sobol.cpp'),'-o',str(probe)],check=True)
        check(probe)


if __name__=='__main__':main()
