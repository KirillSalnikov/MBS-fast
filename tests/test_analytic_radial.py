#!/usr/bin/env python3
"""Exact long-phase primitives versus independent oscillatory quadrature."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn

ROOT=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='mbs_radial_') as directory:
    probe=Path(directory)/'probe'
    subprocess.run(['g++','-std=gnu++11','-O3','-fopenmp','-I'+str(ROOT/'src'),
        str(ROOT/'tests/analytic_radial_probe.cpp'),str(ROOT/'src/AnalyticBackscatter.cpp'),
        '-o',str(probe)],check=True)
    points=[(degree,z)for z in (64.,128.,512.,2048.,8192.)for degree in (0,2,4,64,128,320)]
    stream=''.join('{} {}\n'.format(degree,z)for degree,z in points)
    result=np.loadtxt(subprocess.check_output([str(probe)],input=stream,text=True).splitlines())
    worst=0
    for (degree,z),actual in zip(points,result):
        value=quad(lambda s:s*spherical_jn(degree,s)*np.log(z/s),0,z,epsabs=1e-10,limit=5000)[0]
        derivative=quad(lambda s:s*spherical_jn(degree,s),0,z,epsabs=1e-10,limit=5000)[0]/z
        error=max(abs(actual[0]-value),abs(actual[1]-derivative))
        assert error<3e-8,(degree,z,actual,value,derivative,error)
        worst=max(worst,error)
    print('PASS: 30 long-phase primitives/derivatives, degrees0..320, z64..8192; max absolute error',worst)
