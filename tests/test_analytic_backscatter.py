#!/usr/bin/env python3
"""Independent phase quadrature plus full-solver residual/CLI regression.

Requires numpy/scipy for physical reference integration, not for the solver.
"""
import argparse
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn, sici

ROOT = Path(__file__).resolve().parents[1]


def amplitude(beta, n, B, L):
    s, c = math.sin(beta), math.cos(beta)
    gz, gx = math.sqrt(n*n-s*s), math.sqrt(n*n-c*c)
    def fresnel(mu, eta):
        g = complex(eta*eta-1+mu*mu)**.5
        return np.array([(mu-g)/(mu+g), (eta*eta*mu-g)/(eta*eta*mu+g)])
    ab = (1-fresnel(c,n)**2)*fresnel(gz/n,1/n)*fresnel(s/n,1/n)
    ax = (1-fresnel(s,n)**2)*fresnel(gx/n,1/n)*fresnel(c/n,1/n)
    def window(d, w):
        return max(0, min(d,w)-max(0,d-w))
    pref_b, pref_x = c*window(2*L*s/gz,B), s*window(2*B*c/gx,L)
    return pref_b*ab, pref_x*ax, L*(gz-c)-B*(gx-s)


def self_reference(beta, n, B, L, p, q):
    s,c=math.sin(beta),math.cos(beta)
    g=math.sqrt(n*n-s*s)
    def fresnel(mu,eta):
        root=complex(eta*eta-1+mu*mu)**.5
        return np.array([(mu-root)/(mu+root),(eta*eta*mu-root)/(eta*eta*mu+root)])
    dx=2*p*L*s/g-2*(q-1)*B
    width=max(0,min(dx,B)-max(0,dx-B))
    coeff=(1-fresnel(c,n)**2)*fresnel(g/n,1/n)**(2*p-1)*fresnel(s/n,1/n)**(2*q-1)
    return float(c*c*width*width*np.sum(abs(coeff)**2)/2)


def direct_reference(n, B, L, H, wave, order=1):
    k = 2*math.pi/wave
    points = [0., math.pi/2]
    for v in [B/(2*L), B/L]:
        t = n*v/math.sqrt(1+v*v)
        if t < 1: points.append(math.asin(t))
    for v in [L/(2*B), L/B]:
        t = n*v/math.sqrt(1+v*v)
        if t < 1: points.append(math.acos(t))
    if n < math.sqrt(2):
        t = math.asin(math.sqrt(n*n-1))
        points.extend([t, math.pi/2-t])
    for p in range(1,order+1):
        for q in range(1,order+2-p):
            for v in [2*q-2,2*q-1,2*q]:
                for width,depth,side in [(B,L,False),(L,B,True)]:
                    ratio=v*width/(2*p*depth)
                    t=n*ratio/math.sqrt(1+ratio*ratio)
                    if t<1:points.append(math.acos(t) if side else math.asin(t))
    points = sorted(set(points))
    def intensity(beta):
        a, b, d = amplitude(beta,n,B,L)
        self_value=0.
        for p in range(1,order+1):
            for q in range(1,order+2-p):
                if (p,q)==(1,1):continue
                self_value+=self_reference(beta,n,B,L,p,q)+self_reference(math.pi/2-beta,n,L,B,p,q)
        return float(np.sum(abs(a*np.exp(2j*k*d)+b)**2)/2)+self_value
    integral = sum(quad(intensity,a,b,epsabs=1e-8,epsrel=1e-9,limit=2000)[0]
                   for a,b in zip(points,points[1:]))
    leading = H/(8*math.pi*wave)*integral
    az = k*H
    finite_z = 2/az*(sici(2*az)[0]+(math.cos(2*az)-1)/(2*az))
    return leading*finite_z/(math.pi/az), leading


def data(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--binary', type=Path, default=ROOT/'cpu/bin/mbs_po_mpi')
    args = parser.parse_args()
    binary = str(args.binary.resolve())
    with tempfile.TemporaryDirectory(prefix='mbs_analytic_backscatter_') as folder:
        work = Path(folder)
        probe = work/'probe'
        subprocess.run(['g++','-std=c++11','-O2','-Wall','-Wextra','-pedantic',
                        '-I'+str(ROOT/'src'),str(ROOT/'tests/analytic_backscatter_probe.cpp'),
                        str(ROOT/'src/AnalyticBackscatter.cpp'),'-o',str(probe)], check=True)
        rows = subprocess.check_output([str(probe),'special'],text=True).splitlines()
        for row in rows:
            items = row.split(); x = float(items[1])
            if items[0] == 'sinc':
                expected = 2. if x == 0 else 2/x*(sici(2*x)[0]+(math.cos(2*x)-1)/(2*x))
                if x < .001:  # independent smooth angular quadrature, avoid cancellation
                    expected = quad(lambda t:float(np.sinc(x*t/math.pi)**2),-1,1)[0]
                assert abs(float(items[2])-expected) < 2e-12*max(expected,1e-20), row
            else:
                ell = int(items[2]); expected = spherical_jn(ell,x)
                assert abs(float(items[3])-expected) < 2e-12*max(abs(expected),1e-280), row
        print('PASS: independent Si/sinc and spherical Bessel references')
        for n in [1.31,1.53]:
            for scale in [1.,30.,100.]:
                B,L,H = scale,1.2*scale,.8*scale
                output = subprocess.check_output([str(probe),'strip',str(n),str(B),str(L),str(H),'.532','.6','.8','0'],text=True)
                mean,leading,point = map(float,output.split())
                expected_mean,expected_leading = direct_reference(n,B,L,H,.532)
                assert abs(mean/expected_mean-1) < 1e-6, (n,scale,mean,expected_mean)
                assert abs(leading/expected_leading-1) < 1e-6
                a,b,d = amplitude(math.atan2(.6,.8),n,B,L)
                expected_point = H*H/.532**2*np.sum(abs(a*np.exp(4j*math.pi/.532*d)+b)**2)/2
                assert abs(point-expected_point) < 1e-12*max(expected_point,1e-20)
        print('PASS: physical angular Fresnel and coherent phase means against direct integration')
        for n in [1.31,1.53]:
            for order in [2,3]:
                B,L,H=30.,36.,24.
                output=subprocess.check_output([str(probe),'strip',str(n),str(B),str(L),str(H),'.532','.6','.8','0',str(order)],text=True)
                mean,leading,point=map(float,output.split())
                expected_mean,expected_leading=direct_reference(n,B,L,H,.532,order)
                assert abs(mean/expected_mean-1)<1e-6,(n,order,mean,expected_mean)
                assert abs(leading/expected_leading-1)<1e-6
        print('PASS: longer image-family self means against independent integration')
        env = {k:v for k,v in os.environ.items() if not k.startswith('MBS_')}
        base = [binary,'--method','po','--backend','cpu','--particle','1','4','4',
                '--refractive-index','1.31','0','--wavelength-um','.532','--max-reflections','4',
                '--symmetry','1','1','--scattering-grid','178','180','2','2',
                '--cutoff-profile','off','--no-shadow-output','--close']
        for rule in [['--sobol-seed','64','42'],['--euler-quadrature','4','6'],
                     ['--lattice','32'],['--hammersley','32']]:
            raw_grid = None
            for label,flags,threads in [('raw',[],'1'),('hybrid1',['--analytic-backscatter'],'1'),
                                        ('hybrid4',['--analytic-backscatter'],'4')]:
                target = work/(rule[0][2:]+'_'+label)
                command = base+rule+['--threads',threads,'--output',str(target)]+flags
                with target.with_suffix('.log').open('w') as log:
                    subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
                grid = np.loadtxt(target/(target.name+'.dat'),skiprows=1)
                no_shadow = np.loadtxt(target/(target.name+'_noshadow.dat'),skiprows=1)
                if label == 'raw':
                    raw_grid,raw_ns = grid,no_shadow
                    continue
                summary = data(target/(target.name+'_analytic_backscatter.tsv'))[0]
                samples = data(target/(target.name+'_analytic_backscatter_samples.tsv'))
                w = np.array([float(r['weight']) for r in samples])
                y = np.array([float(r['raw_M11']) for r in samples])
                control = np.array([float(r['control_M11']) for r in samples])
                known = float(summary['analytic_control_mean'])
                assert int(summary['strips']) == 12
                expected_known = 12*direct_reference(1.31,math.sqrt(3)*2,4.,2.,.532)[0]
                assert abs(known/expected_known-1) < 1e-6
                # Validate detected facet normals/edge axes and Euler conventions
                # against the regular hexagon's explicit twelve physical strips.
                for sample in samples[:4]:
                    beta,gamma = math.radians(float(sample['beta_deg'])),math.radians(float(sample['gamma_deg']))
                    u = np.array([-math.sin(beta)*math.cos(gamma),math.sin(beta)*math.sin(gamma),math.cos(beta)])
                    expected_control=0.
                    for sign in [-1.,1.]:
                        for angle in np.arange(6)*math.pi/3+math.pi/6:
                            nx=np.array([math.cos(angle),math.sin(angle),0.]);nz=np.array([0.,0.,sign])
                            sx,sz,edge = u@nx,u@nz,u@np.cross(nx,nz)
                            if sx <= 0 or sz <= 0: continue
                            a,b,d=amplitude(math.atan2(sx,sz),1.31,math.sqrt(3)*2,4.)
                            value=2.**2/.532**2*np.sum(abs(a*np.exp(4j*math.pi/.532*d)+b)**2)/2
                            expected_control+=value*np.sinc(4/.532*edge)**2
                    assert abs(float(sample['control_M11'])-expected_control)<1e-10*max(expected_control,1e-12)
                assert abs(w.sum()-1) < 1e-12
                assert abs(w@y-raw_grid[-1,2]) < 1e-11*raw_grid[-1,2]
                assert abs(known+w@(y-control)-grid[-1,2]) < 1e-11*max(abs(grid[-1,2]),1.)
                correction = known-w@control
                assert abs(no_shadow[-1,2]-raw_ns[-1,2]-correction) < 1e-10
                # All other directions/elements retain the complete MBS calculation.
                difference = grid-raw_grid; difference[-1,2] = 0
                assert np.max(abs(difference)) < 1e-11*np.max(abs(raw_grid[:,2:]))
                if label == 'hybrid1': previous = grid
                else: assert np.max(abs(previous-grid)) < 1e-11*np.max(abs(grid[:,2:]))
        print('PASS: complete coherent residual, unequal weights, no-shadow, thread equivalence, four orientation rules')
        reduced = base.copy(); index=reduced.index('--symmetry');del reduced[index:index+3]
        target=work/'reduced_domain'
        completed=subprocess.run(reduced+['--sobol-seed','32','42','--analytic-backscatter',
                                          '--threads','2','--output',str(target)],
                                 env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
        assert completed.returncode == 0, completed.stdout
        assert 'beta_sym=90 deg, gamma_sym=60 deg' in completed.stdout
        reduced_known=float(data(target/(target.name+'_analytic_backscatter.tsv'))[0]['analytic_control_mean'])
        assert abs(reduced_known/known-1)<1e-12
        print('PASS: verified native hexagonal symmetry retains the same full-sphere control mean')
        # A closed rectangular prism also validates the physical longer-family
        # interpretation directly against ray tracing, rather than only the
        # control identity. Loader facet IDs are deterministic for this fixture.
        box=work/'box.particle'
        faces=[[[0,0,1.2],[1,0,1.2],[1,.8,1.2],[0,.8,1.2]],
               [[1,0,0],[1,.8,0],[1,.8,1.2],[1,0,1.2]],
               [[0,0,0],[0,.8,0],[1,.8,0],[1,0,0]],
               [[0,0,0],[0,0,1.2],[0,.8,1.2],[0,.8,0]],
               [[0,.8,0],[0,.8,1.2],[1,.8,1.2],[1,.8,0]],
               [[0,0,0],[1,0,0],[1,0,1.2],[0,0,1.2]]]
        box.write_text('0\n0\n180 360\n\n'+'\n\n'.join('\n'.join(' '.join(map(str,v)) for v in f) for f in faces)+'\n')
        for beta,p,q in [(25.,2,1),(40.,2,2),(25.,3,2)]:
            target=work/f'audit_{beta:g}_{p}_{q}'
            audit=work/(target.name+'.tsv')
            audit_env=dict(env,MBS_COHERENCE_AUDIT=str(audit),MBS_FORCE_TRACK_IDS='1')
            command=[binary,'--method','po','--backend','cpu','--particle-file',str(box),
                     '--refractive-index','1.31','0','--wavelength-um','.532','--max-reflections','12',
                     '--fixed-orientation',str(beta),'180','--scattering-grid','179.999','180','1','1',
                     '--cutoff-profile','off','--threads','1','--allow-experimental-environment',
                     '--output',str(target),'--close']
            subprocess.run(command,env=audit_env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True,timeout=30)
            from io import StringIO
            records=list(csv.DictReader(StringIO(audit.read_text().replace('\\t','\t').lstrip('# ')),delimiter='\t'))
            matrix=np.zeros((2,2),complex); rays=0
            for row in records:
                path=row['path'].split(',')
                if path[0]!='0' or path[-1]!='0':continue
                internal=path[1:-1]
                if internal.count('5')!=p or internal.count('0')!=p-1 or internal.count('2')!=q or internal.count('3')!=q-1:continue
                if len(internal)!=2*(p+q)-2:continue
                for a in range(2):
                    for b in range(2):matrix[a,b]+=complex(float(row[f'j{a}{b}_re']),float(row[f'j{a}{b}_im']))
                rays+=1
            native=float(np.sum(abs(matrix)**2)/2)
            expected=float(subprocess.check_output([str(probe),'self',str(math.radians(beta)),'1.31','1','1.2','.8','.532',str(p),str(q)],text=True))
            assert rays>0 and abs(native/expected-1)<1e-7,(beta,p,q,rays,native,expected)
            print(f'  native family beta={beta:g}, p={p}, q={q}: rays={rays}, relative S11 error={abs(native/expected-1):.3g}')
        print('PASS: longer return-family coherent self intensities match actual MBS ray groups')
        target=work/'long_return_control'
        long_command=base+['--sobol-seed','32','42','--analytic-backscatter',
                           '--analytic-return-order','2','--threads','2','--output',str(target)]
        index=long_command.index('--max-reflections');long_command[index+1]='12'
        subprocess.run(long_command,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True,timeout=30)
        summary=data(target/(target.name+'_analytic_backscatter.tsv'))[0]
        expected=12*direct_reference(1.31,math.sqrt(3)*2,4.,2.,.532,2)[0]
        assert int(summary['return_order'])==2 and abs(float(summary['analytic_control_mean'])/expected-1)<1e-6
        assert abs(float(summary['hybrid_M11'])-float(summary['analytic_control_mean'])-float(summary['residual_mean']))<1e-10
        print('PASS: native extended control preserves full residual and known spherical mean')
        point_command=base.copy();index=point_command.index('--scattering-grid')
        point_command[index:index+5]=['--scattering-grid','0','2','2']
        target=work/'point_cone'
        subprocess.run(point_command+['--sobol-seed','32','42','--analytic-backscatter',
                                     '--threads','2','--output',str(target)],
                       env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True,timeout=30)
        point_grid=np.loadtxt(target/(target.name+'.dat'),skiprows=1)
        summary=data(target/(target.name+'_analytic_backscatter.tsv'))[0]
        assert np.max(abs(point_grid[:,0]-180))<1e-12
        assert np.max(abs(point_grid[:,2]-float(summary['hybrid_M11'])))<1e-10
        print('PASS: zero-radius cone corrects every repeated exact-backward row')
        for flags,text in [(['--scattering-grid','178','179','2','1'],'exact 180-degree'),
                           (['--particle','2','4','4'],'no right-angle edges'),
                           (['--symmetry','3','6'],'cannot verify')]:
            # Replace the respective existing option rather than pass it twice.
            command = base.copy(); key=flags[0]; index=command.index(key)
            count=4 if key == '--scattering-grid' else (2 if key == '--symmetry' else 3)
            command[index:index+count+1]=flags
            target=work/('reject_'+key[2:])
            completed=subprocess.run(command+['--sobol','8','--analytic-backscatter','--output',str(target)],
                                     env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=30)
            assert completed.returncode != 0 and text in completed.stdout, completed.stdout
        print('PASS: unsupported scattering grid and geometry rejected')
    print('Analytic backscatter regression passed')


if __name__ == '__main__':
    main()
