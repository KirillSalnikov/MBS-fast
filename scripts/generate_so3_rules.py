#!/usr/bin/env python3
"""Positive Haar rules; Lebedev x Fourier or Wigner-zero-mode projection.

The Wigner rule integrates D^j_mn exactly through j=degree: Fourier sums
remove m,n!=0 and Gauss-Legendre integrates the remaining P_j(cos beta).
Lebedev replaces this polar/azimuth product with spherical cubature; a full
Fourier alpha rule remains necessary. Global Haar rotations preserve the
algebraic order and provide independent randomized replicas for comparison.
"""
import argparse
import json
from pathlib import Path
import warnings
import numpy as np
from scipy.integrate import lebedev_rule
from scipy.special import roots_legendre, eval_legendre, sph_harm_y
from scipy.spatial.transform import Rotation


def make_rule(method,degree):
    if degree<1:raise ValueError('degree must be positive')
    count=degree+1
    periodic=2*np.pi*(np.arange(count)+.5)/count
    if method=='lebedev':
        xyz,w=lebedev_rule(degree)
        if np.min(w)<0:raise ValueError('This implementation uses positive rules; choose a positive Lebedev degree.')
        beta=np.arccos(np.clip(xyz[2],-1,1));gamma=np.mod(np.arctan2(xyz[1],xyz[0]),2*np.pi)
        sphere_weight=w/(4*np.pi)
    else:
        mu,w=roots_legendre((degree+2)//2)
        beta=np.repeat(np.arccos(mu),count);gamma=np.tile(periodic,len(mu))
        sphere_weight=np.repeat(w/2/count,count)
    angles=np.column_stack([np.tile(periodic,len(beta)),np.repeat(beta,count),np.repeat(gamma,count)])
    weights=np.repeat(sphere_weight,count)/count
    assert abs(weights.sum()-1)<1e-12 and np.min(weights)>=0
    # Check algebraic moments before applying any random rotation.
    maximum=0.
    for j in range(1,degree+1):
        maximum=max(maximum,abs(float(np.dot(sphere_weight,eval_legendre(j,np.cos(beta))))))
        if method=='lebedev':
            for m in sorted(set([0,1,j//2,j])):
                maximum=max(maximum,float(abs(np.dot(sphere_weight,sph_harm_y(j,m,beta,gamma)))))
    for m in range(1,degree+1):maximum=max(maximum,float(abs(np.exp(1j*m*periodic).mean())))
    assert maximum<2e-11,(method,degree,maximum)
    return angles,weights,maximum


def write_rule(path,angles,weights,seed=None):
    if seed is not None:
        generator=np.random.default_rng(seed)
        # Rotation.random uses a uniform distribution on SO(3).
        rotation=Rotation.random(random_state=generator)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',UserWarning)
            angles=(rotation*Rotation.from_euler('ZYZ',angles)).as_euler('ZYZ')
    columns=np.column_stack([np.rad2deg(angles[:,1]),np.mod(np.rad2deg(angles[:,2]),360),
                             np.mod(np.rad2deg(angles[:,0]),360),weights])
    np.savetxt(path,columns,fmt='%.17g',header='BETA_DEG GAMMA_DEG ALPHA_DEG WEIGHT')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--method',choices=['lebedev','wigner'],required=True)
    p.add_argument('--degree',type=int,required=True)
    p.add_argument('--replicas',type=int,default=8)
    p.add_argument('--seed',type=int,default=22001)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    angles,weights,error=make_rule(a.method,a.degree)
    files=[]
    write_rule(a.output/'unrotated.rule',angles,weights)
    for i in range(a.replicas):
        path=a.output/('replica_'+str(i)+'.rule');write_rule(path,angles,weights,a.seed+i*101);files.append(str(path.resolve()))
    report=dict(method=a.method,algebraic_degree=a.degree,orientations=len(weights),replicas=a.replicas,
                normalized_weight_sum=float(weights.sum()),max_checked_moment_error=error,files=files,
                wigner_operation='Integral projection onto D^0_00; not reconstruction of all spectral coefficients.',
                sources=['https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.lebedev_rule.html',
                         'https://doi.org/10.1109/LSP.2015.2490676'])
    (a.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
