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

def interpolation_components(predicted,direct,ratios):
    delta=predicted-direct;pm=predicted.mean(0);dm=direct.mean(0)
    bound=(abs(delta.mean(0))+T*delta.std(0,ddof=1)/math.sqrt(8))/dm[:,0,None]
    ratio_bound=np.zeros_like(bound)
    if ratios:
        p=pm/pm[:,0,None];d=dm/dm[:,0,None]
        influence=(predicted-p[None]*predicted[:,:,0,None])/pm[None,:,0,None]-(direct-d[None]*direct[:,:,0,None])/dm[None,:,0,None]
        ratio_bound=abs(p-d)+T*influence.std(0,ddof=1)/math.sqrt(8)
        ratio_bound[:,0]=0
    return bound,ratio_bound

def interpolation_bounds(predicted,direct,ratios):
    raw,normalized=interpolation_components(predicted,direct,ratios)
    return np.maximum(raw,normalized)

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

def reconstruct_angular(root,state,angles,seeds):
    requested=grid(root/'requested_theta.csv');support=grid(root/'angular_audit_support.csv')
    frozen=np.loadtxt(root/'angular_audit_budget.csv',delimiter=',',skiprows=1,ndmin=2)
    assert np.allclose(frozen[:,0],angles,rtol=0,atol=1e-10)
    anchor_budget=frozen[:,2::2].astype(int)
    records=[line.split() for line in (root/'angular_audit_sources.tsv').read_text().splitlines()]
    anchor_records=[' '.join(line[2:]) for line in records if line[0]=='anchor']
    anchor,seen,_=reconstruct(root,state,support,anchor_budget[nearest_rows(angles,support)],None,seeds,anchor_records)
    result=np.array([local_spline(support,level)(angles) for level in anchor]);assigned=np.isin(angles,support)
    width=min(4,len(support))
    for group in json.loads((root/'angular_audit_groups.json').read_text())['groups']:
        rows=requested[np.array(group['theta_indices'],int)]
        source=[' '.join(line[2:]) for line in records if line[0]=='residual' and int(line[1])==group['id']]
        if not source:continue
        budget=np.tile(np.array(group['counts'],int),(len(rows),1))
        data,checked,_=reconstruct(root,state,rows,budget,None,seeds,source);seen.update(checked)
        for interval in group['intervals']:
            left=min(max(interval-1,0),len(support)-width);right=min(interval,len(support)-width)
            knots=support[left:right+width];positions=nearest_rows(rows,knots)
            assert np.allclose(rows[positions],knots,rtol=0,atol=1e-10)
            target=np.flatnonzero((angles>support[interval])&(angles<support[interval+1]))
            assert not assigned[target].any();assigned[target]=True
            raw=nearest_rows(rows,angles[target]);assert np.allclose(rows[raw],angles[target],rtol=0,atol=1e-10)
            for level in range(len(result)):
                predicted=local_spline(knots,data[level][:,positions])(angles[target])
                result[level][:,target]+=data[level][:,raw]-predicted
    assert assigned.all(),'angular audit skipped checked nodes'
    return result,seen,[]

def common_bounds(root,state,label,evaluated,support,seeds,actual_budget,separate=False):
    budget=np.loadtxt(root/('interpolation_'+label+'_budgets.csv'),delimiter=',',skiprows=1,ndmin=2)
    lines=(root/('interpolation_'+label+'_sources.tsv')).read_text().splitlines();bound=np.zeros((len(evaluated),16));ratio_bound=np.zeros_like(bound);seen=set();width=min(4,len(support))
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
            raw,normalized=interpolation_components(predicted,direct,state['ratio_accuracy'])
            bound[target]=np.maximum(bound[target],raw);ratio_bound[target]=np.maximum(ratio_bound[target],normalized)
    return ((bound,ratio_bound) if separate else np.maximum(bound,ratio_bound)),seen

def check_error_budget(state,contract,requested,checked,ci,ri,oci,ori,raw,ratios):
    """Independent acceptance reconstruction; never interpolate error estimates."""
    assert np.isfinite(raw).all() and np.isfinite(ratios).all()
    assert max(raw.max(),ratios.max())<=state['interpolation_budget']
    if contract.get('error_budget_policy')!='pointwise-stat-plus-interpolation-v1':
        assert max(ci.max(),oci.max())<=state['statistical_budget']
        if state['ratio_accuracy']:assert max(ri.max(),ori.max())<=state['statistical_budget']
        return None
    target=contract['fixed_error_reserve']+state['statistical_budget']+state['interpolation_budget']
    reserve=contract['fixed_error_reserve']
    output_raw=np.full_like(oci,state['interpolation_budget']);output_ratio=np.full_like(ori,state['interpolation_budget']);output_ratio[:,0]=0
    ix=nearest_rows(requested,checked);assert np.allclose(requested[ix],checked,rtol=0,atol=1e-10)
    output_raw[ix]=raw;output_ratio[ix]=ratios
    direct=ci+raw;output=oci+output_raw
    if state['ratio_accuracy']:
        direct=np.maximum(direct,ri+ratios);output=np.maximum(output,ori+output_ratio)
    assert np.isfinite(direct).all() and np.isfinite(output).all()
    assert max(direct.max(),output.max())<=state['statistical_budget']+state['interpolation_budget']
    return dict(direct=float(direct.max()+reserve),output=float(output.max()+reserve),target=target)

def run(root,reference,dense=False,dense_workers=1):
    state=json.loads((root/'fullauto_status.json').read_text());assert state['status']=='converged',state
    contract=json.loads((root/'precision_contract.json').read_text())
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
    raw_guard,ratio_guard=interpolation_components(spline(evaluated),samples,state['ratio_accuracy'])
    if json.loads((root/'interpolation_production_mode.json').read_text())['common_prefix']:
        (raw_guard,ratio_guard),checked=common_bounds(root,state,'production',evaluated,support,state['seeds'],budget[:,1:].astype(int),True);seen.update(checked)
    guard=np.maximum(raw_guard,ratio_guard)
    assert np.max(abs(guard.max(1)-np.loadtxt(root/'theta_validation.csv',delimiter=',',skiprows=1)[:,3]))<1e-11
    pointwise=check_error_budget(state,contract,requested,evaluated,confidence,ratio_ci,output_ci,output_ratio,raw_guard,ratio_guard)
    validation=json.loads((root/'independent_final_validation.json').read_text())
    angular=False
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
        angular=validation.get('estimator')=='angular_residual_v1'
        if angular:vlevels,vseen,_=reconstruct_angular(root,state,angles,validation['seeds'])
        else:vlevels,vseen,_=reconstruct(root,state,angles,vb,root/'independent_final_sources.tsv',validation['seeds'])
        seen.update(vseen)
        values=vlevels.sum(0);validation_values=values.copy();vm,vci,vri=intervals(values);vid=np.searchsorted(angles,support)
        vs=local_spline(support,values[:,vid]);delta=vs(angles)-values
        raw_bound,ratio_bound=interpolation_components(vs(angles),values,state['ratio_accuracy'])
        if not angular and json.loads((root/'interpolation_validation_mode.json').read_text())['common_prefix']:
            (raw_bound,ratio_bound),checked=common_bounds(root,state,'validation',angles,support,validation['seeds'],vb,True);seen.update(checked)
        _,voci,vori=intervals(vs(requested))
        independent_pointwise=check_error_budget(state,contract,requested,angles,vci,vri,voci,vori,raw_bound,ratio_bound)
        if independent_pointwise is not None:
            assert abs(independent_pointwise['direct']-validation['max_combined_direct95'])<1e-10
            assert abs(independent_pointwise['output']-validation['max_combined_output95'])<1e-10
    assert state.get('max_estimated_analytic_mean_error_scaled_M11',0)<=state.get('analytic_mean_budget',0)
    report=dict(passed=True,schema=3,pointwise_budget=pointwise,epsilon=state['statistical_budget']+state['interpolation_budget']+state.get('mirror_budget',0)+state.get('analytic_mean_budget',0)+state.get('rotation_budget',0)+state.get('cutoff_budget',0)+state.get('execution_budget',0)+state.get('depth_budget',0),
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
            if not angular and json.loads((root/'interpolation_validation_mode.json').read_text())['common_prefix']:
                bound,checked=common_bounds(root,state,'validation',angles,support,validation['seeds'],vb);seen.update(checked)
            else:
                vid=np.searchsorted(angles,support)
                predicted=local_spline(support,validation_values[:,vid])(angles)
                bound=interpolation_bounds(predicted,validation_values,state['ratio_accuracy'])
        report['independent_dense']=dict(matrix_error_scaled_M11=float(np.max(abs(output_mean-dm)/dm[:,0,None])),residual_plus95_scaled_M11=float(bound.max()),all_requested_nodes_pass=bool(np.all(bound<=state['interpolation_budget'])),scope='interpolation bias on common prefixes, sampling variability checked separately')
        assert report['independent_dense']['all_requested_nodes_pass'],report['independent_dense']
    (root/'independent_fullauto_audit.json').write_text(json.dumps(report,indent=2)+'\n');return report
