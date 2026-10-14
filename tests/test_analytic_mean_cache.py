#!/usr/bin/env python3
"""Lazy authoritative rows and grid-independent refinement of real facets."""
from pathlib import Path
import subprocess
import tempfile
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='mbs_mean_rows_')as directory:
    work=Path(directory);probe=work/'probe'
    subprocess.run(['g++','-std=gnu++11','-O3','-fopenmp','-I'+str(ROOT/'src'),
        '-I'+str(ROOT/'src/math'),'-I'+str(ROOT/'src/handler'),
        str(ROOT/'tests/analytic_facet_probe.cpp'),str(ROOT/'src/AnalyticFacetAverage.cpp'),
        str(ROOT/'src/AnalyticBackscatter.cpp'),str(ROOT/'src/math/Sobol.cpp'),'-o',str(probe)],check=True)
    reference=work/'theta.txt';np.savetxt(reference,np.linspace(0,180,4097),fmt='%.17g')
    cache=work/'means.cache'
    def run(angles,scale='1',ri='1.31'):
        return np.atleast_2d(np.loadtxt(subprocess.check_output([str(probe),'mean-cached',
            'rectangle',scale,ri,'0','circular',str(cache),str(reference),*map(str,angles)],text=True).splitlines()))
    first=run([0,90,180]);header=cache.read_text().splitlines()[0].split()
    assert header[0]=='MBS_PHYSICAL_FACET_MEAN_V3_COMPACT_ROWS' and int(header[2])==3
    run([45,135]);assert int(cache.read_text().splitlines()[0].split()[2])==5
    np.testing.assert_array_equal(run([0,90,180]),first)
    # The same angles must be bitwise identical when computed alone and in a
    # batch that includes a much harder, high-phase theta near backscatter.
    direct=np.loadtxt(subprocess.check_output([str(probe),'mean','rectangle','300',
        '1.3116','0','circular','0','90','180'],text=True).splitlines())
    cache.unlink();subset=run([0,90,180],scale='300',ri='1.3116')
    np.testing.assert_array_equal(subset[:,:3],direct[:,:3])
    for angle,row in zip([0,90,180],direct):
        cold=work/('cold'+str(angle)+'.cache')
        result=np.loadtxt(subprocess.check_output([str(probe),'mean-cached','rectangle','300',
            '1.3116','0','circular',str(cold),'none',str(angle)],text=True).splitlines())
        np.testing.assert_array_equal(result,row[:3])
    changed=run([0,90,180],scale='301',ri='1.3116')
    assert not np.array_equal(changed,subset)
    assert int(cache.read_text().splitlines()[0].split()[2])==3
    print('PASS: only3/4097 reference rows computed; subset append/reuse; bitwise per-row refinement; geometry invalidation')
