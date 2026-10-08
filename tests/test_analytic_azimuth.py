#!/usr/bin/env python3
"""Exact conditional azimuth identity, native residual and backend agreement."""
import argparse,csv,os,subprocess,tempfile
from io import StringIO
from pathlib import Path
import numpy as np
from scipy.special import i0e
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[1]

def table(path):
    with path.open() as file:return list(csv.DictReader(file,delimiter='\t'))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--binary',type=Path,default=ROOT/'cpu/bin/mbs_po_mpi');parser.add_argument('--probe',type=Path);args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='mbs_azimuth_') as directory:
        work=Path(directory);probe=args.probe
        if probe is None:
            probe=work/'probe'
            subprocess.run(['g++','-std=gnu++11','-O2','-fopenmp','-I'+str(ROOT/'src'),
                str(ROOT/'tests/analytic_azimuth_probe.cpp'),str(ROOT/'src/AnalyticAzimuthGaussian.cpp'),'-o',str(probe)],check=True)
        points=np.r_[0,1e-10,np.geomspace(1e-6,1e8,300),49.999,50,50.001]
        values=np.loadtxt(StringIO(subprocess.check_output([str(probe),'i0'],input='\n'.join(map(str,points))+'\n',text=True)))
        relative=np.max(abs(values/i0e(points)-1));assert relative<5e-14,relative
        rng=np.random.default_rng(29871);rows=[]
        for kappa in [0,.01,1.,10.,50.,1e3,1e6]:
            for beam_theta in [0,.2,1.4,np.pi]:
                for delta in [0,.7/np.sqrt(max(kappa,1.)),.3]:
                    t=np.clip(beam_theta+delta,0,np.pi);phi=rng.uniform(-np.pi,np.pi)
                    d=(np.sin(beam_theta)*np.cos(phi),np.sin(beam_theta)*np.sin(phi),-np.cos(beam_theta))
                    rows.append([*d,kappa*.532**2/(2*np.pi),1.7,.532,t,.73])
        stream=''.join(' '.join(map(str,row))+'\n' for row in rows)
        result=np.atleast_2d(np.loadtxt(StringIO(subprocess.check_output([str(probe)],input=stream,text=True))))
        for row,actual in zip(rows,result):
            d=np.array(row[:3]);theta=row[-2];peak,kappa,point,mean=actual[:4]
            v=np.array([np.sin(theta)*np.cos(.73),np.sin(theta)*np.sin(.73),-np.cos(theta)])
            expected=peak*np.exp(-kappa*np.sum((d-v)**2)/2)
            assert abs(point-expected)<1e-12*max(peak,1.)
            B=kappa*np.sin(theta)*np.hypot(d[0],d[1]);scale=max(1.,np.sqrt(B))
            # Resolve narrow azimuth peaks explicitly rather than allowing an
            # adaptive quadrature to miss them. This integrates the original
            # exp(B*cos(phi)) angular factor, independently of the I0 series.
            integral=quad(lambda y:np.exp(-2*B*np.sin(y/(2*scale))**2),0,min(np.pi*scale,16),epsabs=1e-13)[0]/(np.pi*scale)
            beam_theta=np.arctan2(np.hypot(d[0],d[1]),-d[2])
            expected=peak*np.exp(-2*kappa*np.sin((theta-beam_theta)/2)**2)*integral
            assert abs(mean-expected)<2e-12*max(peak,1.),(row,actual,expected)
            if len(actual)>4:
                assert abs(actual[4]-actual[5])<2e-12*max(peak,1.)
                assert abs(mean-actual[6])<2e-12*max(peak,1.)
        print(f'PASS: scaled I0 relative error {relative:.3g}; exact conditional means for narrow/finite beams and poles')
        if args.probe:return
        env={k:v for k,v in os.environ.items() if not k.startswith('MBS_')};env['OMP_NUM_THREADS']='2'
        base=[str(args.binary.resolve()),'--method','po','--backend','cpu','--particle','1','4','4','--wavelength-um','.532',
            '--max-reflections','4','--sobol-seed','64','17','--haar-alpha','--analytic-facet-average',
            '--scattering-grid','0','180','2','12','--threads','2','--cutoff-profile','off','--no-shadow-output','--close']
        for index in [('1.31','0'),('1.53','.0018'),('1','0')]:
            reference=None
            for label,flags in [('base',[]),('azimuth',['--analytic-azimuth-gaussian','--analytic-facet-samples'])]:
                target=work/(index[0]+label)
                with target.with_suffix('.log').open('w') as log:
                    subprocess.run(base+['--refractive-index',*index,'--output',str(target)]+flags,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
                grid=np.loadtxt(target/(target.name+'.dat'),skiprows=1);ns=np.loadtxt(target/(target.name+'_noshadow.dat'),skiprows=1)
                if reference is None:reference=grid;reference_ns=ns;continue
                summary=table(target/(target.name+'_analytic_azimuth.tsv'));samples=table(target/(target.name+'_analytic_azimuth_samples.tsv'))
                for t,row in enumerate(summary):
                    selected=[s for s in samples if abs(float(s['theta_deg'])-float(row['theta_deg']))<1e-9]
                    correction=sum(float(s['weight'])*(float(s['conditional_mean'])-float(s['point'])) for s in selected)
                    scale=max(abs(reference[t,2]),1.)
                    assert abs(grid[t,2]-reference[t,2]-correction)<1e-10*scale
                    assert abs(ns[t,2]-reference_ns[t,2]-correction)<1e-10*scale
                difference=grid-reference;difference[:,2]=0;assert np.max(abs(difference))<1e-10*max(np.max(abs(grid)),1.)
                assert np.max(abs(grid[[0,-1],2]-reference[[0,-1],2]))<1e-10*max(np.max(abs(grid)),1.)
                if index[0]=='1':
                    assert all(float(row['sampled_point'])==0 and float(row['conditional_mean'])==0 for row in summary)
        print('PASS: retained internal-beam control, additive coherent/no-shadow residual, vacuum and unchanged remaining Mueller entries')

if __name__=='__main__':main()
