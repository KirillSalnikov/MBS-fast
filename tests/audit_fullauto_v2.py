"""Independent reconstruction of local-budget native fullauto from worker files."""
import json
import math
from pathlib import Path
import subprocess
import numpy as np
from scipy.interpolate import CubicHermiteSpline

T=2.3646242510102993

def fnv(data):
    value=14695981039346656037
    for byte in data:value=((value^byte)*1099511628211)&((1<<64)-1)
    return format(value,'x')

def grid(path):return np.loadtxt(path,ndmin=1)
def matrix(path):return np.loadtxt(path,skiprows=1,ndmin=2)
def nearest_rows(angles,query):
    right=np.minimum(np.searchsorted(angles,query),len(angles)-1);left=np.maximum(right-1,0)
    return np.where(abs(angles[left]-query)<abs(angles[right]-query),left,right)

def stored_path(path,root):
    path=Path(path)
    if path.exists():return path
    name=str(path)
    if '/jobs/' in name:return root/'jobs'/name.split('/jobs/',1)[1]
    return path

def validate(folder,seen,root):
    folder=stored_path(folder,root)
    if folder in seen:return
    seen.add(folder);meta=(folder/'complete.txt').read_text().split()
    assert fnv((folder/'result/result.dat').read_bytes())==meta[2],folder
    if meta[3]!='none':assert fnv((folder/'result/result_analytic_facets.tsv').read_bytes())==meta[3]
    sources=folder/'angular_sources.tsv'
    if sources.exists():
        assert fnv(sources.read_bytes())==meta[4]
        for line in sources.read_text().splitlines():
            name,checksum=line.split('\t');assert fnv((stored_path(name,root)/'complete.txt').read_bytes())==checksum
            validate(name,seen,root)
    paired=folder/'paired_depth.tsv'
    if paired.exists():
        assert fnv(paired.read_bytes())==meta[4]
        for line in paired.read_text().splitlines():
            name,checksum=line.split();assert fnv(stored_path(name,root).read_bytes())==checksum

def intervals(samples):
    mean=samples.mean(0);variance=samples.var(0,ddof=1);scale=mean[:,0]
    assert np.all(scale>0)
    confidence=T*np.sqrt(variance/8)/scale[:,None]
    ratios=mean/scale[:,None];ratio_ci=np.zeros_like(mean);q=T*T/8
    for j in range(1,16):
        covariance=((samples[:,:,0]-scale)*(samples[:,:,j]-mean[:,j])).sum(0)/7
        aa=scale*scale-q*variance[:,0]
        transformed=(samples[:,:,j]-ratios[:,j]*samples[:,:,0]).var(0,ddof=1)
        disc=np.maximum(0,q*scale*scale*transformed-q*q*np.maximum(0,variance[:,0]*variance[:,j]-covariance*covariance))
        centre=(scale*mean[:,j]-q*covariance)/aa
        ratio_ci[:,j]=np.where(aa>0,abs(centre-ratios[:,j])+np.sqrt(disc)/aa,1e100)
    return mean,confidence,ratio_ci

def interpolation_bounds(predicted,direct,ratios):
    delta=predicted-direct;pm=predicted.mean(0);dm=direct.mean(0)
    bound=(abs(delta.mean(0))+T*delta.std(0,ddof=1)/math.sqrt(8))/dm[:,0,None]
    if ratios:
        p=pm/pm[:,0,None];d=dm/dm[:,0,None]
        influence=(predicted-p[None]*predicted[:,:,0,None])/pm[None,:,0,None]-(direct-d[None]*direct[:,:,0,None])/dm[None,:,0,None]
        ratio_bound=abs(p-d)+T*influence.std(0,ddof=1)/math.sqrt(8)
        bound=np.maximum(bound,ratio_bound)
    return bound

def local_spline(x,y):
    derivative=np.zeros_like(y);width=min(4,len(x))
    for i in range(len(x)):
        first=min(max(i-1,0),len(x)-width)
        for j in range(first,first+width):
            if j==i:coefficient=sum(1/(x[i]-x[k]) for k in range(first,first+width) if k!=i)
            else:
                coefficient=1/(x[j]-x[i])
                for k in range(first,first+width):
                    if k!=i and k!=j:coefficient*=(x[i]-x[k])/(x[j]-x[k])
            derivative[:,i]+=coefficient*y[:,j]
    return CubicHermiteSpline(x,y,derivative,axis=1)

def reconstruct(root,state,angles,budgets,source,seeds,records=None):
    levels=budgets.shape[1];data=np.zeros((levels,8,len(angles),16));filled=np.zeros((levels,2,8,len(angles)),bool)
    coefficients=matrix(root/'polarization_weights.tsv')
    seen=set();jobs=[]
    for line in source.read_text().splitlines() if records is None else records:
        level,sign,count,seed,name=line.split();level,sign,count,seed=map(int,[level,sign,count,seed]);slot=seeds.index(seed)
        folder=stored_path(name,root);validate(folder,seen,root);raw=matrix(folder/'result/result.dat');values=raw[:,2:].copy()
        if state['analytic_controls']:
            control=matrix(folder/'result/result_analytic_facets.tsv')
            for j in range(1,16):
                wr=np.interp(raw[:,0],coefficients[:,0],coefficients[:,1+2*j]);ws=np.interp(raw[:,0],coefficients[:,0],coefficients[:,2+2*j])
                values[:,j]+=wr*(control[:,2]-control[:,4])+ws*(control[:,3]-control[:,5])
        target=np.flatnonzero(budgets[:,level]==count);rows=nearest_rows(raw[:,0],angles[target]);ok=rows<len(raw)
        ok &= np.isclose(raw[np.minimum(rows,len(raw)-1),0],angles[target],rtol=0,atol=1e-10)
        target=target[ok];rows=rows[ok];sign_slot=0 if sign==1 else 1
        assert not filled[level,sign_slot,slot,target].any(),(level,sign,seed)
        filled[level,sign_slot,slot,target]=True;data[level,slot,target]+=sign*values[rows]
        jobs.append(dict(level=level,sign=sign,count=count,seed=seed,folder=folder,args=(folder/'command.args').read_text().splitlines()))
    assert filled[:,0].all() and filled[1:,1].all()
    return data,seen,jobs

def common_bounds(root,state,label,evaluated,support,seeds,actual_budget):
    budget=np.loadtxt(root/('interpolation_'+label+'_budgets.csv'),delimiter=',',skiprows=1,ndmin=2)
    lines=(root/('interpolation_'+label+'_sources.tsv')).read_text().splitlines();bound=np.zeros((len(evaluated),16));seen=set();width=min(4,len(support))
    for group in np.unique(budget[:,0]).astype(int):
        rows=budget[budget[:,0]==group];angles=rows[:,1];records=[' '.join(line.split()[1:]) for line in lines if int(line.split()[0])==group]
        levels,checked,_=reconstruct(root,state,angles,rows[:,2:].astype(int),None,seeds,records);seen.update(checked);data=levels.sum(0)
        for interval in range(len(support)-1):
            left=min(max(interval-1,0),len(support)-width);right=min(interval,len(support)-width)
            stencil=support[left:right+width];target=np.flatnonzero((evaluated>=support[interval])&(evaluated<=support[interval+1]))
            all_nodes=np.concatenate([stencil,evaluated[target]]);ns=actual_budget[nearest_rows(evaluated,all_nodes)].max(0)
            if not np.array_equal(ns,rows[0,2:]):continue
            positions=nearest_rows(angles,stencil);assert np.allclose(angles[positions],stencil,rtol=0,atol=1e-10)
            predicted=local_spline(stencil,data[:,positions])(evaluated[target]);direct=data[:,nearest_rows(angles,evaluated[target])]
            bound[target]=np.maximum(bound[target],interpolation_bounds(predicted,direct,state['ratio_accuracy']))
    return bound,seen

def run(root,reference,dense=False,dense_workers=1):
    state=json.loads((root/'fullauto_status.json').read_text());assert state['status']=='converged',state
    requested=grid(root/'requested_theta.csv');evaluated=grid(root/'theta_evaluated.csv');support=grid(root/'theta_support.csv')
    budget=np.loadtxt(root/'orientation_budget.csv',delimiter=',',skiprows=1,ndmin=2);assert np.allclose(budget[:,0],evaluated,rtol=0,atol=1e-10)
    levels,seen,jobs=reconstruct(root,state,evaluated,budget[:,1:].astype(int),root/'final_stage_sources.tsv',state['seeds']);samples=levels.sum(0)
    mean,confidence,ratio_ci=intervals(samples);indices=np.searchsorted(evaluated,support);assert np.allclose(evaluated[indices],support,rtol=0,atol=1e-10)
    spline=local_spline(support,samples[:,indices]);output=spline(requested)
    output_mean,output_ci,output_ratio=intervals(output)
    difference=float(np.max(abs(output_mean-matrix(root/'mueller_fullauto.dat')[:,2:])/output_mean[:,0,None]))
    direct_difference=float(np.max(abs(mean-matrix(root/'mueller_evaluated.dat')[:,2:])/mean[:,0,None]))
    saved_ci=np.loadtxt(root/'all_mueller_confidence.csv',delimiter=',',skiprows=1)[:,1:]
    saved_ratio=np.loadtxt(root/'normalized_mueller_confidence.csv',delimiter=',',skiprows=1)[:,1:]
    assert difference<1e-11 and direct_difference<1e-11
    assert np.max(abs(output_ci-saved_ci))<1e-11
    assert np.max(abs(output_ratio-saved_ratio))<1e-10
    residual=spline(evaluated)-samples
    guard=interpolation_bounds(spline(evaluated),samples,state['ratio_accuracy'])
    if json.loads((root/'interpolation_production_mode.json').read_text())['common_prefix']:
        guard,checked=common_bounds(root,state,'production',evaluated,support,state['seeds'],budget[:,1:].astype(int));seen.update(checked)
    assert np.max(abs(guard.max(1)-np.loadtxt(root/'theta_validation.csv',delimiter=',',skiprows=1)[:,3]))<1e-11
    assert guard.max()<=state['interpolation_budget'] and max(confidence.max(),output_ci.max())<=state['statistical_budget']
    if state['ratio_accuracy']:assert max(ratio_ci.max(),output_ratio.max())<=state['statistical_budget']
    validation=json.loads((root/'independent_final_validation.json').read_text())
    if validation['enabled']:
        assert validation['passed'] and set(validation['seeds']).isdisjoint(state['seeds'])
        angles=grid(root/'independent_theta.csv');vb=np.empty((len(angles),budget.shape[1]-1),int)
        for i,angle in enumerate(angles):
            right=min(np.searchsorted(evaluated,angle),len(evaluated)-1);left=max(0,right-1)
            if abs(evaluated[right]-angle)<1e-10:vb[i]=budget[right,1:]
            else:vb[i]=np.maximum(budget[left,1:],budget[right,1:])
        if (root/'independent_orientation_budget.csv').exists():
            frozen=np.loadtxt(root/'independent_orientation_budget.csv',delimiter=',',skiprows=1,ndmin=2)
            assert np.allclose(frozen[:,0],angles,rtol=0,atol=1e-10)
            vb=frozen[:,1:].astype(int)
        vlevels,vseen,_=reconstruct(root,state,angles,vb,root/'independent_final_sources.tsv',validation['seeds']);seen.update(vseen)
        values=vlevels.sum(0);validation_values=values.copy();vm,vci,vri=intervals(values);vid=np.searchsorted(angles,support)
        vs=local_spline(support,values[:,vid]);delta=vs(angles)-values
        bound=interpolation_bounds(vs(angles),values,state['ratio_accuracy'])
        if json.loads((root/'interpolation_validation_mode.json').read_text())['common_prefix']:
            bound,checked=common_bounds(root,state,'validation',angles,support,validation['seeds'],vb);seen.update(checked)
        assert bound.max()<=state['interpolation_budget']
        assert max(vci.max(),intervals(vs(requested))[1].max())<=state['statistical_budget']
        if state['ratio_accuracy']:assert max(vri.max(),intervals(vs(requested))[2].max())<=state['statistical_budget']
    assert state.get('max_estimated_analytic_mean_error_scaled_M11',0)<=state.get('analytic_mean_budget',0)
    report=dict(passed=True,schema=2,epsilon=state['statistical_budget']+state['interpolation_budget']+state.get('mirror_budget',0)+state.get('analytic_mean_budget',0)+state.get('rotation_budget',0)+state.get('cutoff_budget',0)+state.get('execution_budget',0),
                requested_theta=len(requested),evaluated_theta=len(evaluated),support_theta=len(support),final_counts_per_seed=state['counts_per_seed'],
                minimum_counts_per_seed=budget[:,1:].min(0).astype(int).tolist(),phi=state['phi'],complete_jobs_audited=len(seen),
                matrix_difference_scaled_M11=difference,evaluated_matrix_difference_scaled_M11=direct_difference,
                statistical95_scaled_M11=float(max(confidence.max(),output_ci.max())),ratio95=float(max(ratio_ci.max(),output_ratio.max())),
                guard_residual95_scaled_M11=float(guard.max()),independent_final_validation=validation,
                limits='pointwise statistical estimates; frozen target depth/cutoffs; mean refinement and physical sensitivity reported separately')
    if dense:
        # Recalculate each budget group at all requested nodes. Preserve
        # paired scramble prefixes; interpolate only the explicit budget map.
        dense_budgets=np.empty((len(requested),budget.shape[1]-1),int)
        for i,angle in enumerate(requested):
            right=min(np.searchsorted(evaluated,angle),len(evaluated)-1);left=max(0,right-1)
            dense_budgets[i]=budget[right,1:] if abs(evaluated[right]-angle)<1e-10 else np.maximum(budget[left,1:],budget[right,1:])
        full=np.zeros((8,len(requested),16));coefficients=matrix(root/'polarization_weights.tsv');destination=root/'independent_dense';destination.mkdir(exist_ok=True)
        for level in range(dense_budgets.shape[1]):
            for count in np.unique(dense_budgets[:,level]):
                target=np.flatnonzero(dense_budgets[:,level]==count);angles=requested[target]
                if len(angles)==1:angles=np.unique(np.append(angles,requested[0] if angles[0]!=requested[0] else requested[-1]))
                group_grid=destination/('grid_{}_{}.txt'.format(level,count));np.savetxt(group_grid,angles,fmt='%.17g')
                weights=destination/('weights_{}_{}.tsv'.format(level,count));np.savetxt(weights,np.column_stack([angles,np.interp(angles,coefficients[:,0],coefficients[:,1]),np.interp(angles,coefficients[:,0],coefficients[:,2])]),fmt='%.17g',header='theta_deg reflection_weight shadow_weight',comments='')
                for seed in state['seeds']:
                    for sign in [1,-1] if level else [1]:
                        template=next(j for j in jobs if j['level']==level and j['sign']==sign and j['seed']==seed)
                        args=template['args'].copy()
                        def setarg(flag,value):
                            if flag in args:args[args.index(flag)+1]=str(value)
                        # A native paired source is replayed as a standalone
                        # requested depth here, independently of the cache.
                        if '--paired-depth' in args:i=args.index('--paired-depth');del args[i:i+2]
                        if '--sobol-offset' in args:i=args.index('--sobol-offset');del args[i:i+2]
                        i=args.index('--sobol-seed');args[i+1]=str(count)
                        setarg('--theta-grid-file',group_grid);setarg('--analytic-control-weights',weights)
                        folder=destination/('l{}_n{}_s{}_sign{}'.format(level,count,seed,sign));setarg('--output',folder/'result')
                        if not (folder/'result/result.dat').exists():
                            folder.mkdir(exist_ok=True);(folder/'command.args').write_text('\n'.join(args)+'\n')
                            with (folder/'stdout.log').open('w') as stream:subprocess.run(args,stdout=stream,stderr=subprocess.STDOUT,check=True)
                        raw=matrix(folder/'result/result.dat');values=raw[:,2:].copy()
                        if state['analytic_controls']:
                            c=matrix(folder/'result/result_analytic_facets.tsv')
                            for j in range(1,16):values[:,j]+=np.interp(raw[:,0],coefficients[:,0],coefficients[:,1+2*j])*(c[:,2]-c[:,4])+np.interp(raw[:,0],coefficients[:,0],coefficients[:,2+2*j])*(c[:,3]-c[:,5])
                        rows=nearest_rows(raw[:,0],requested[target]);full[state['seeds'].index(seed),target]+=sign*values[rows]
        dm=full.mean(0);delta=output-full;bound=interpolation_bounds(output,full,state['ratio_accuracy'])
        if json.loads((root/'interpolation_production_mode.json').read_text())['common_prefix']:
            # The actual heterogeneous estimator's sampling difference is
            # already covered by its statistical budget. The interpolation
            # bias is audited on common prefixes over each local stencil.
            validation=json.loads((root/'independent_final_validation.json').read_text())
            angles=grid(root/'independent_theta.csv')
            assert len(angles)==len(requested) and np.allclose(angles,requested,rtol=0,atol=1e-10),'Use all-node independent validation for a dense common-prefix audit'
            vb=np.empty_like(dense_budgets)
            for i,angle in enumerate(angles):
                right=min(np.searchsorted(evaluated,angle),len(evaluated)-1);left=max(0,right-1)
                vb[i]=budget[right,1:] if abs(evaluated[right]-angle)<1e-10 else np.maximum(budget[left,1:],budget[right,1:])
            if (root/'independent_orientation_budget.csv').exists():
                frozen=np.loadtxt(root/'independent_orientation_budget.csv',delimiter=',',skiprows=1,ndmin=2)
                assert np.allclose(frozen[:,0],angles,rtol=0,atol=1e-10)
                vb=frozen[:,1:].astype(int)
            if json.loads((root/'interpolation_validation_mode.json').read_text())['common_prefix']:
                bound,checked=common_bounds(root,state,'validation',angles,support,validation['seeds'],vb);seen.update(checked)
            else:
                vid=np.searchsorted(angles,support)
                predicted=local_spline(support,validation_values[:,vid])(angles)
                bound=interpolation_bounds(predicted,validation_values,state['ratio_accuracy'])
        report['independent_dense']=dict(matrix_error_scaled_M11=float(np.max(abs(output_mean-dm)/dm[:,0,None])),residual_plus95_scaled_M11=float(bound.max()),all_requested_nodes_pass=bool(np.all(bound<=state['interpolation_budget'])),scope='interpolation bias on common prefixes, sampling variability checked separately')
        assert report['independent_dense']['all_requested_nodes_pass'],report['independent_dense']
    (root/'independent_fullauto_audit.json').write_text(json.dumps(report,indent=2)+'\n');return report
