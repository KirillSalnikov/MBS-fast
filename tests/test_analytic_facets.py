#!/usr/bin/env python3
"""Independent Haar/phase references and native all-angle residual regression."""
import argparse,csv,math,os,subprocess,tempfile
from pathlib import Path
import numpy as np
from scipy.special import eval_legendre,roots_legendre
from scipy.stats import qmc

ROOT=Path(__file__).resolve().parents[1]

def polygon(shape,scale):
    if shape=='rectangle':p=np.array([[0,0],[1,0],[1,.6],[0,.6]])
    elif shape=='triangle':p=np.array([[0,0],[1,0],[.2,.7]])
    else:p=5.0510479784*np.column_stack((np.cos(np.arange(6)*np.pi/3),np.sin(np.arange(6)*np.pi/3)))
    return scale*p

def transform(p,qx,qy):
    # Direct barycentric area integral of each triangle, independent of the
    # production boundary formula. exprel(i*x)=exp(i*x/2)*sinc(x/2).
    result=np.zeros(qx.shape,complex)
    for j in range(1,len(p)-1):
        v=p[[0,j,j+1]];a=qx*(v[1,0]-v[0,0])+qy*(v[1,1]-v[0,1]);b=qx*(v[2,0]-v[0,0])+qy*(v[2,1]-v[0,1])
        phase=qx*v[0,0]+qy*v[0,1]
        jac=float((v[1,0]-v[0,0])*(v[2,1]-v[0,1])-(v[1,1]-v[0,1])*(v[2,0]-v[0,0]))
        swap=abs(a)>abs(b);aa=np.where(swap,b,a);bb=np.where(swap,a,b)
        small=np.maximum(abs(aa),abs(bb))<1e-3
        value=np.empty(qx.shape,complex)
        value[small]=.5+1j*(aa[small]+bb[small])/6-(aa[small]**2+bb[small]**2+aa[small]*bb[small])/24
        a0,b0=aa[~small],bb[~small]
        rel=lambda z:np.exp(.5j*z)*np.sinc(z/(2*np.pi))
        value[~small]=(np.exp(1j*b0)*rel(a0-b0)-rel(a0))/(1j*b0)
        result+=jac*np.exp(1j*phase)*value
    return result

def R(mu,index):
    g=np.sqrt(index*index-1+mu*mu+0j)
    return .5*(abs((mu-g)/(mu+g))**2+abs((index*index*mu-g)/(index*index*mu+g))**2)

def native_R(mu,index):
    # MBS preserves a real refracted-ray direction even for complex index,
    # then uses n*cosB in its external Jones coefficients. This convention
    # differs slightly from the ideal complex sqrt Fresnel law used by C.
    sine2=1-mu*mu;delta=(index*index).real-sine2
    q=.5*(delta+math.hypot(delta,2*index.real*index.imag))
    cosB=math.sqrt(q/(sine2+q))
    return .5*(abs((mu-index*cosB)/(mu+index*cosB))**2+abs((index*mu-cosB)/(index*mu+cosB))**2)

def direct_mean(p,index,theta,nmu=256,nphi=1024):
    h,c=math.sin(theta/2),math.cos(theta/2)
    if abs(theta-math.pi)<1e-13:h,c=1.,0.
    intervals=[(-1.,1.)] if theta==0 else [(-1.,-c),(-c,c),(c,1.)]
    x,w=roots_legendre(nmu);all_z=[];all_w=[]
    for a,b in intervals:
        if b-a<1e-13:continue
        all_z.extend((a+b)/2+(b-a)/2*x);all_w.extend(w*(b-a)/2)
    z=np.array(all_z);wz=np.array(all_w)
    nx,nw=roots_legendre(128);A=h*z;B=c*np.sqrt(np.maximum(0,1-z*z))
    limit=np.where(A>=B,np.pi,np.where(A<=-B,0,np.arccos(np.clip(-A/np.maximum(B,1e-300),-1,1))))
    mu=A[:,None]+B[:,None]*np.cos((nx[None,:]+1)*limit[:,None]/2)
    nu=2*A[:,None]-mu
    angular=h*h/np.pi*limit/2*np.sum(nw[None,:]*np.maximum(0,c*c+mu*nu)*R(np.maximum(mu,0),index),axis=1)
    phi=2*np.pi*(np.arange(nphi)+.5)/nphi
    reflected=0.;shadow=0.;k=2*np.pi/.532
    for start in range(0,len(z),32):
        zz=z[start:start+32];rr=np.sqrt(1-zz*zz)[:,None]
        for Q,weight,kind in [(2*k*h,angular[start:start+32],'ref'),(k*math.sin(theta),c**4*(1-zz*zz)/4,'shadow')]:
            amplitude=transform(p,Q*rr*np.cos(phi),Q*rr*np.sin(phi))
            value=float(np.sum(wz[start:start+32]*weight*np.mean(abs(amplitude)**2,axis=1)))/2/.532**2
            if kind=='ref':reflected+=value
            else:shadow+=value
    return np.array([reflected,shadow])

def table(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--binary',type=Path,default=ROOT/'cpu/bin/mbs_po_mpi');args=parser.parse_args()
    binary=str(args.binary.resolve())
    with tempfile.TemporaryDirectory(prefix='mbs_general_analytic_') as directory:
        work=Path(directory);probe=work/'probe'
        subprocess.run(['g++','-std=gnu++11','-O2','-march=native','-Wall','-Wextra','-fopenmp',
          '-I'+str(ROOT/'src'),'-I'+str(ROOT/'src/math'),'-I'+str(ROOT/'src/handler'),
          str(ROOT/'tests/analytic_facet_probe.cpp'),str(ROOT/'src/AnalyticFacetAverage.cpp'),
          str(ROOT/'src/AnalyticBackscatter.cpp'),str(ROOT/'src/math/Sobol.cpp'),'-o',str(probe)],check=True)
        points=np.loadtxt(subprocess.check_output([str(probe),'sobol','42','0','31'],text=True).splitlines())
        assert np.array_equal(points[:,:3],qmc.Sobol(3,scramble=False).random_base2(5)[1:])
        scrambled=np.loadtxt(subprocess.check_output([str(probe),'sobol','42','1','127'],text=True).splitlines())
        assert np.array_equal(scrambled[:,:2],scrambled[:,3:])
        print('PASS: third Sobol direction numbers match SciPy; old scrambled beta/gamma samples preserved')
        rng=np.random.default_rng(3791);mu=rng.uniform(1e-5,1,500);phi=rng.uniform(0,2*np.pi,500);theta=rng.uniform(0,np.pi,500);az=rng.uniform(0,2*np.pi,500)
        stream=''.join(' '.join(map(str,r))+'\n' for r in zip(mu,phi,theta,az))
        actual=np.array([list(map(float,r.split())) for r in subprocess.check_output([str(probe),'projection'],input=stream,text=True).splitlines()])
        nu=np.sqrt(1-mu*mu)*np.cos(phi-az)*np.sin(theta)-mu*np.cos(theta);h=np.sin(theta/2)
        assert np.max(abs(actual-h[:,None]**2*(1-h[:,None]**2+mu[:,None]*nu[:,None])))<1e-12
        print('PASS: arbitrary-angle physical polarization weights agree with native Jones projection')
        x,w=roots_legendre(512);m=(x+1)/2;w=w/2;spin=2*np.pi*(np.arange(256)+.5)/256
        for index in [1.31+0j,1.53+.0018j]:
            for degrees in [0,30,90,170,180]:
                th=math.radians(degrees);h,c=math.sin(th/2),math.cos(th/2);t=2*h*h-1
                coef=np.array(list(map(float,subprocess.check_output([str(probe),'moments',str(index.real),str(index.imag),str(degrees),'64'],text=True).splitlines())))
                z=h*m[:,None]+c*np.sqrt(1-m*m)[:,None]*np.cos(spin)
                nu=t*m[:,None]+2*h*c*np.sqrt(1-m*m)[:,None]*np.cos(spin)
                for l in [0,2,4,16,64]:
                    ref=(2*l+1)/2*np.sum(w*np.mean(h*h*(c*c+m[:,None]*nu)*R(m[:,None],index)*eval_legendre(l,z),axis=1))
                    assert abs(coef[l]-ref)<1e-12,(index,degrees,l,coef[l],ref)
        print('PASS: two-direction Legendre moments agree with independent angular-spin integration')
        for shape,scale,index in [('rectangle',1.,1.31+0j),('triangle',3.,1.53+.0018j),('hex',1.,1.31+0j)]:
            degrees=[0,30,90,170,180]
            native=np.loadtxt(subprocess.check_output([str(probe),'mean',shape,str(scale),str(index.real),str(index.imag),'facets',*map(str,degrees)],text=True).splitlines())
            for row,deg in zip(native,degrees):
                reference=direct_mean(polygon(shape,scale),index,math.radians(deg))
                fine=direct_mean(polygon(shape,scale),index,math.radians(deg),384,1536)
                assert np.max(abs(reference-fine)/np.maximum(abs(fine),1e-100))<2e-4,(shape,deg,reference,fine)
                error=np.max(abs(row[1:3]-fine)/np.maximum(abs(fine),1e-100))
                assert error<2e-4,(shape,deg,row,fine,error)
        print('PASS: reflected and shadow means at0/30/90/170/180 versus original Haar phase integrals')
        for index in [1.31+0j,1.53+.0018j]:
            for degrees in [30,90,170,180]:
                grouped=np.loadtxt(subprocess.check_output([str(probe),'groups',str(index.real),str(index.imag),str(degrees)],text=True).splitlines())
                assert np.max(abs(grouped[:2]-grouped[2:]))<1e-12*max(np.max(abs(grouped)),1.),(index,degrees,grouped)
        print('PASS: congruent-facet grouping equals independent self means under cyclic order, translation, tilt and reversal')
        # Cold means on enough rows to enter the OpenMP path, including an
        # oscillatory large facet and an absorbing non-rectangular facet.
        # Equality covers the per-row sums and the convergence reduction.
        for shape,scale,index,shadow,angles in [
            ('rectangle',300.,1.31+0j,'off',np.linspace(0,25,129)),
            ('triangle',4.,1.53+.0018j,'facets',np.linspace(0,180,33)),
            ('hex',1.,1.31+0j,'circular',np.linspace(0,180,33))]:
            command=[str(probe),'mean',shape,str(scale),str(index.real),str(index.imag),shadow,*map(str,angles)]
            outputs=[subprocess.check_output(command,env=dict(os.environ,OMP_NUM_THREADS=str(n))) for n in [1,4]]
            assert outputs[0]==outputs[1],(shape,'parallel analytic means changed')
        print('PASS: serial/four-thread analytic means are bitwise equal on large, absorbing and circular-shadow cases')
        env={k:v for k,v in os.environ.items() if not k.startswith('MBS_')}
        base=[binary,'--method','po','--backend','cpu','--particle','1','4','4','--refractive-index','1.31','0',
              '--wavelength-um','.532','--max-reflections','4','--scattering-grid','0','180','2','6',
              '--no-shadow-output','--cutoff-profile','off','--close']
        for rule in [['--sobol-seed','64','42','--haar-alpha'],['--so3-full-quaternion','64']]:
            raw_grid=None
            for name,flags,threads in [('raw',[],'1'),('facet1',['--analytic-facet-average','--analytic-facet-samples'],'1'),
                                       ('facet4',['--analytic-facet-average','--analytic-facet-samples'],'4')]:
                target=work/(rule[0][2:]+'_'+name)
                if flags:flags=flags+['--analytic-mean-cache',str(work/(rule[0][2:]+'_means.cache'))]
                with target.with_suffix('.log').open('w') as log:
                    subprocess.run(base+rule+flags+['--threads',threads,'--output',str(target)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=120)
                grid=np.loadtxt(target/(target.name+'.dat'),skiprows=1);ns=np.loadtxt(target/(target.name+'_noshadow.dat'),skiprows=1)
                if name=='raw':raw_grid=grid;raw_ns=ns;continue
                summary=table(target/(target.name+'_analytic_facets.tsv'));samples=table(target/(target.name+'_analytic_facet_samples.tsv'))
                for i,row in enumerate(summary):
                    selected=[s for s in samples if abs(float(s['theta_deg'])-float(row['theta_deg']))<1e-9]
                    weights=np.array([float(s['weight']) for s in selected]);y=np.array([float(s['raw_M11']) for s in selected])
                    ref=np.array([float(s['reflection_control']) for s in selected]);shadow=np.array([float(s['shadow_control']) for s in selected])
                    known=float(row['reflection_mean'])+float(row['shadow_mean'])
                    assert abs(weights@y-raw_grid[i,2])<1e-10*max(abs(raw_grid[i,2]),1e-20)
                    assert abs(known+weights@(y-ref-shadow)-grid[i,2])<1e-10*max(abs(grid[i,2]),1.)
                    assert abs(float(row['reflection_mean'])-weights@ref+raw_ns[i,2]-ns[i,2])<1e-10*max(abs(ns[i,2]),1.)
                difference=grid-raw_grid;difference[:,2]=0
                assert np.max(abs(difference))<1e-10*np.max(abs(raw_grid[:,2:]))
                if name=='facet1':first=grid
                else:assert np.max(abs(first-grid))<1e-10*np.max(abs(first[:,2:]))
        print('PASS: all-angle complete coherent residual, no-shadow, full quaternion and one/four threads')
        # Confirm the optional circular model and the two distinct analytic
        # controls compose on the SAME fields rather than overwriting output.
        target=work/'combined'
        subprocess.run(base+['--sobol-seed','32','42','--haar-alpha','--analytic-facet-average','--analytic-shadow-control','circular',
          '--analytic-backscatter','--analytic-facet-samples','--threads','2','--output',str(target)],env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True,timeout=120)
        facets=table(target/(target.name+'_analytic_facets.tsv'));returns=table(target/(target.name+'_analytic_backscatter.tsv'))[0]
        grid=np.loadtxt(target/(target.name+'.dat'),skiprows=1)
        correction=float(returns['analytic_control_mean'])-float(returns['sampled_control_mean'])
        assert abs(float(facets[-1]['hybrid_M11'])+correction-grid[-1,2])<1e-10
        assert np.max(abs(grid[:-1,2]-np.array([float(r['hybrid_M11']) for r in facets[:-1]])))<1e-10
        print('PASS: circular shadow option and old analytic backscatter combine additively')
        # Actual first reflected ray apertures of the closed unit cube. This
        # checks the phase vector, Fresnel law and native physical projection
        # together; no analytic-control estimator is involved in these runs.
        beta,gamma=math.radians(37),math.radians(19)
        rotation=np.array([[math.cos(beta)*math.cos(gamma),-math.cos(beta)*math.sin(gamma),math.sin(beta)],
                           [math.sin(gamma),math.cos(gamma),0],
                           [-math.sin(beta)*math.cos(gamma),math.sin(beta)*math.sin(gamma),math.cos(beta)]])
        u=rotation.T@np.array([0.,0.,1.])
        from io import StringIO
        for index in [1.31+0j,1.53+.0018j]:
            for degrees in [1e-6,30.,90.,170.,180.]:
                target=work/f'ray_{index.real}_{degrees:g}';audit=target.with_suffix('.tsv')
                audit_env=dict(env,MBS_COHERENCE_AUDIT=str(audit),MBS_FORCE_TRACK_IDS='1')
                command=[binary,'--method','po','--backend','cpu','--particle-file',str(ROOT/'examples/cube.particle'),
                    '--refractive-index',str(index.real),str(index.imag),'--wavelength-um','.532','--max-reflections','3',
                    '--fixed-orientation','37','19','--scattering-grid',str(max(0,degrees-.001)),str(degrees),'1','1',
                    '--threads','1','--cutoff-profile','off','--allow-experimental-environment','--output',str(target),'--close']
                subprocess.run(command,env=audit_env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True,timeout=30)
                rows=list(csv.DictReader(StringIO(audit.read_text().replace('\\t','\t').lstrip('# ')),delimiter='\t'))
                native=0.;shadow=0.;rays=0
                for row in rows:
                    intensity=sum(float(row[f'j{a}{b}_{part}'])**2 for a in range(2) for b in range(2) for part in ['re','im'])/2
                    if row['external']=='1':shadow+=intensity
                    elif ',' not in row['path'] and row['path']!='':native+=intensity;rays+=1
                theta=math.radians(degrees);v=rotation.T@np.array([math.sin(theta),0.,-math.cos(theta)])
                q=2*np.pi/.532*(u+v);h,c=math.sin(theta/2),math.cos(theta/2);expected=0.
                for axis in range(3):
                    for sign in [-1,1]:
                        mu=sign*u[axis]
                        if mu<=0:continue
                        nu=sign*v[axis];tangent=[i for i in range(3) if i!=axis]
                        fourier=np.prod(np.sinc(q[tangent]/(2*np.pi)))
                        expected+=h*h*max(0,c*c+mu*nu)*native_R(mu,index)*fourier*fourier/.532**2
                assert rays==3 and abs(native-expected)<1e-7*max(expected,1e-10),(index,degrees,native,expected,rays)
                print(f'  native first reflections: index={index}, theta={degrees:g}, relative error={abs(native-expected)/max(expected,1e-10):.3g}')
                if degrees<1e-5:
                    expected_shadow=np.sum(abs(u))**2/.532**2
                    assert abs(shadow/expected_shadow-1)<1e-9,(shadow,expected_shadow)
        print('PASS: real MBS reflected ray self intensities at five angles for real/absorbing indices; coherent forward shadow')
        for index in [1.53+.0018j,1.+0j]:
            absorption_base=base.copy();position=absorption_base.index('--refractive-index')
            absorption_base[position+1:position+3]=[str(index.real),str(index.imag)]
            reference=None
            for label,options in [('raw',[]),('analytic',['--analytic-facet-average'])]:
                target=work/f'index_{index.real}_{label}'
                subprocess.run(absorption_base+['--sobol-seed','32','7','--haar-alpha',
                    '--threads','1','--output',str(target)]+options,env=env,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,check=True,timeout=120)
                grid=np.loadtxt(target/(target.name+'.dat'),skiprows=1)
                if label=='raw':reference=grid;continue
                rows=table(target/(target.name+'_analytic_facets.tsv'))
                expected=np.array([float(r['hybrid_M11']) for r in rows])
                assert np.max(abs(grid[:,2]-expected))<1e-10*max(np.max(abs(grid[:,2])),1.)
                if index==1+0j:assert np.array_equal(grid,reference)
        print('PASS: absorbing full solver and exact vacuum null-control regression')
    print('General analytic facet regression passed')

if __name__=='__main__':main()
