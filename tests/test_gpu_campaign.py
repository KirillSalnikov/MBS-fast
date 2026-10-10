#!/usr/bin/env python3
"""Independent regression tests for control calibration and telescoping."""
import sys
from pathlib import Path
import unittest
import tempfile
import threading
import time

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from gpu_campaign import JobRunner, allocation, confidence, fit_control_weights, reweight
from run_gpu_optimized import combine, difference, validate_level_order


class EstimatorTests(unittest.TestCase):
    def test_dynamic_queue_retains_seed_order_and_one_job_per_device(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);binary=root/'binary';binary.write_bytes(b'test')
            theta=root/'theta.csv';theta.write_text('170\n180\n')
            runner=JobRunner(binary,root/'run',theta_file=theta,gpus=[0,1,2],controls='off',scheduling='queue')
            active=set();lock=threading.Lock();assigned=[]
            def job(gpu,n,seed,depth,cutoff,phi,**options):
                with lock:
                    self.assertNotIn(gpu,active);active.add(gpu);assigned.append((gpu,seed))
                time.sleep(.003 if gpu==0 else .02)
                matrix=np.zeros((2,18));matrix[:,0]=[170,180];matrix[:,1]=1;matrix[:,2]=seed+100
                path=root/f'{seed}.dat';np.savetxt(path,matrix,comments='',header='columns')
                with lock:active.remove(gpu)
                return dict(data=str(path),seconds=.003 if gpu==0 else .02,gpu=gpu)
            runner.job=job
            data,_,_,jobs=runner.stage(128,seeds=list(range(8)))
            np.testing.assert_array_equal(data[:,0,2],np.arange(8)+100)
            self.assertEqual(sorted(seed for _,seed in assigned),list(range(8)))
            self.assertGreater(sum(gpu==0 for gpu,_ in assigned),3)

    def test_equal_depth_cutoff_refinement_has_a_distinct_order(self):
        validate_level_order([8,12,18,18],['.001','.001','.001','off'])
        validate_level_order([8,12,18],['.001','off','off'])
        for depths,cutoffs in [([18,18],['off','off']),([18,18],['off','.001']),
                               ([18,18],['.001','.01']),([18,12],['.001','off'])]:
            with self.assertRaises(ValueError):validate_level_order(depths,cutoffs)

    def test_test3_preset_includes_exact_backscatter(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);binary=root/'binary';binary.write_bytes(b'test')
            runner=JobRunner(binary,root/'run',case='test3',gpus=[1,2,3],controls='off')
            self.assertEqual((runner.height,runner.diameter),(316.2,123.8))
            self.assertEqual(len(runner.theta),227)
            self.assertEqual((runner.theta[0],runner.theta[-1]),(170.,180.))
            self.assertEqual(runner.physical['analytic_controls'],'off')

    def test_gpu_scheduler_covers_eight_seeds_with_one_three_or_four_devices(self):
        # Exercise real batching with synthetic completed jobs. Check every
        # seed is used exactly once and no physical device runs two jobs.
        for devices in ([3],[1,2,3],[0,1,2,3]):
            with tempfile.TemporaryDirectory() as temporary:
                root=Path(temporary);binary=root/'binary';binary.write_bytes(b'test')
                theta=root/'theta.txt';theta.write_text('170\n180\n')
                runner=JobRunner(binary,root/'run',theta_file=theta,gpus=devices,controls='off')
                seen=[];active=set();lock=threading.Lock()
                def job(gpu,n,seed,depth,cutoff,phi,**options):
                    with lock:
                        self.assertNotIn(gpu,active);active.add(gpu);seen.append(seed)
                    matrix=np.zeros((2,18));matrix[:,0]=[170,180];matrix[:,1]=1;matrix[:,2]=seed+100
                    path=root/f'{seed}.dat';np.savetxt(path,matrix,header='columns',comments='')
                    with lock:active.remove(gpu)
                    return dict(data=str(path),seconds=1.)
                runner.job=job
                samples,summaries,wall,jobs=runner.stage(128,seeds=list(range(8)))
                self.assertEqual(sorted(seen),list(range(8)))
                self.assertEqual(len(jobs),8)
                self.assertEqual(wall,(8+len(devices)-1)//len(devices))
                np.testing.assert_array_equal(samples[:,0,2],np.arange(8)+100)
                for sample,summary in zip(samples,summaries):
                    np.testing.assert_array_equal(summary['raw_M11'],sample[:,2])
                    self.assertTrue(np.all(summary['sampled_reflection']==0))

    def test_joint_seed_covariance_is_retained(self):
        base=np.zeros((8,2,18));base[:,:,0]=[0,180];base[:,:,1]=[1,1]
        x=np.arange(8)-3.5
        base[:,:,2]=100+x[:,None]
        delta=base.copy();delta[:,:,2:]=0;delta[:,:,2]=-x[:,None]
        final=combine([base,delta])
        metrics,mean,_=confidence(final)
        self.assertEqual(metrics['max_M11_pointwise95'],0.)
        np.testing.assert_array_equal(mean[:,2],[100,100])

    def test_full_mueller_telescoping_not_amplitude_averaging(self):
        rng=np.random.default_rng(73)
        values=[]
        for _ in range(3):
            d=np.zeros((8,3,18));d[:,:,0]=[0,15,180];d[:,:,1]=.01
            d[:,:,2:]=rng.normal(size=(8,3,16));d[:,:,2]+=100
            values.append(d)
        result=combine([values[0],difference(values[1],values[0]),difference(values[2],values[1])])
        np.testing.assert_allclose(result,values[2],rtol=0,atol=4e-14)

    def test_independent_control_fit_on_known_polynomial_integral(self):
        # Uniform x in[-1,1]: E[x]=0, E[x²]=1/3. The physical signal
        # 100+3*x²+2*x has known exact mean101 and positive intensity.
        def sample(seed):
            rng=np.random.default_rng(seed);x=rng.uniform(-1,1,(8,128))
            c1=(x*x).mean(axis=1);c2=x.mean(axis=1);y=100+3*c1+2*c2
            dtype=[(n,float) for n in ['raw_M11','reflection_mean','shadow_mean','sampled_reflection','sampled_shadow']]
            summaries=[];data=np.zeros((8,1,18));data[:,:,0]=30;data[:,:,1]=1
            for i in range(8):
                s=np.zeros(1,dtype=dtype);s['raw_M11']=y[i];s['reflection_mean']=1/3
                s['sampled_reflection']=c1[i];s['sampled_shadow']=c2[i];summaries.append(s)
            data[:,:,2]=np.array(y-c1+1/3-c2)[:,None]
            return data,summaries
        train,st=sample(111);valid,sv=sample(222)
        weights,_=fit_control_weights(train,valid,st,sv)
        production,sp=sample(333)
        corrected=reweight(production,sp,weights)
        self.assertLess(corrected[:,:,2].var(),production[:,:,2].var()/5)
        self.assertLess(abs(corrected[:,:,2].mean()-101),.05)
        # Exact controls with known coefficients remove the complete sampled
        # dependence and recover the analytic integral on a third sample.
        exact=reweight(production,sp,np.array([[3.,2.]]))
        np.testing.assert_allclose(exact[:,:,2],101,rtol=0,atol=3e-14)

    def test_allocation_prefers_variance_and_cheap_levels(self):
        d=np.zeros((8,2,18));d[:,:,0]=[0,180];d[:,:,1]=1
        d[:,:,2]=100+(np.arange(8)-3.5)[:,None]
        delta=d.copy();delta[:,:,2:]=0;delta[:,:,2]=(np.arange(8)-3.5)[:,None]/10
        counts=allocation([d,delta],[1,10],[128,128],np.array([100,100]),.03)
        self.assertGreater(counts[0],counts[1])
        self.assertTrue(all(n>=1 and n&(n-1)==0 for n in counts))


if __name__=='__main__':unittest.main()
