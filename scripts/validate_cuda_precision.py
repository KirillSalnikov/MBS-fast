#!/usr/bin/env python3
from pathlib import Path
import importlib.util,json,subprocess,os,argparse
import sys
from types import SimpleNamespace
parser=argparse.ArgumentParser(description='Validate explicit reduced trigonometry against precise FP32 on full Mueller curves; not a physical-convergence or timing certificate')
parser.add_argument('--prefix',type=Path,default=Path.home()/'.local/opt/mbs-fast')
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--shapes',help='Optional comma-separated shape names')
parser.add_argument('--sizes',default='.01,1,20,100',help='Comma-separated k_eq values')
options=parser.parse_args()
root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('b',root/'scripts/benchmark_hardware.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
out=options.output.expanduser().absolute();out.mkdir(parents=True,exist_ok=False)
a=SimpleNamespace(prefix=options.prefix,profile='gpu-fp32',binary=None,python=None,threads=8,device=0,output=out,label='RTX5080 reduced profile validation')
suite=json.loads((root/'benchmarks/hardware_v1/suite.json').read_text());build=b.build_info(a,suite);candidate=Path(build['native_binary']).parent/'mbs_po_gpu_float_reduced';env=b.environment(a,build)
shapes=[('column',['1','1.8','1.4']),('plate',['1','.15','1.4']),('bullet',['2','1.8','1.4']),('bullet-rosette',['3','1.8','1.4']),('droxtal',['4','1']),('concave-column',['10','1.8','1.4','20']),('aggregate2',['12','1.8','1.4','2']),('certain-aggregate',['999','1'])]
report=dict(passed=False,precise_build=build,candidate_binary=str(candidate),candidate_sha256=b.sha(candidate),candidate_version=b.capture([candidate,'--version'],env),angles=181,phi_points=16,orientations=64,cases=[])
selected=set(options.shapes.split(',')) if options.shapes else {name for name,particle in shapes}
if selected-{name for name,particle in shapes}:parser.error('Unknown shape name')
sizes=[float(v) for v in options.sizes.split(',')]
if not sizes or any(not (0<v<float('inf')) for v in sizes):parser.error('Positive finite k_eq values required')
b.validate_mbs_version('gpu-fp32-reduced',report['candidate_version'])
report['driver_sha256']=b.sha(Path(__file__))
for name,particle in shapes:
 if name not in selected:continue
 for size in sizes:
  c=dict(id=name+'-'+str(size),kind='mbs',angles=181,work=dict(orientations=64),arguments=['--method','po','--particle']+particle+['--k-eq',str(size),'--ri','1.3116','0','--wavelength-um','.532','--sobol-seed','64','314159','--haar-alpha','--scattering-grid','0','180','16','180','--max-reflections','8','--cutoff-profile','off','--orientation-chunk','64','--threads','8','--profile-phases'])
  matrices=[];runs=[]
  for label,binary in [('precise',build['native_binary']),('reduced',str(candidate))]:
   directory=out/c['id']/label;command=b.command_for(c,build,'gpu-fp32',8,directory);command[0]=binary
   r=b.timed_process(command,env,directory,0,180);matrix,metrics,work=b.read_result(c,directory,'gpu-fp32');matrices.append(matrix);runs.append(dict(profile=label,wall_s=r['wall_s'],command=command,work=work,metrics=metrics))
  error=b.matrix_error(matrices[1],matrices[0]);passed=error['max_scaled_M11']<2e-4 and error['relative_L2']<5e-5
  record=dict(shape=name,k_eq=size,validation=error,passed=passed,runs=runs);report['cases'].append(record);b.write(out/'results.json',report);print(name,size,error,'PASS' if passed else 'FAIL',flush=True)
  if not passed:raise RuntimeError('Reduced phase matrix gate failed')
# Existing absorbing irregular-particle cancellation regression, both coherent
# and incoherent, with its independent 3e-5 matrix/ratio gate.
p=subprocess.run([sys.executable,str(root/'tests/test_gpu_reduced_profiles.py'),'--reference',build['native_binary'],'--candidate',str(candidate)],env=env,text=True,capture_output=True,check=True)
report['absorbing_dust_regression']=p.stdout;print(p.stdout,flush=True)
# The range reducer is also tested against long-double libm, independent of
# either scattering executable.
host=Path(build['installation']['compiler']);native_source=Path(build['source_root']);cpu=out/'phase_host_test';gpu=out/'phase_device_test'
subprocess.run([str(host),'-O2','-std=c++11',str(native_source/'tests/test_gpu_phase_reduction.cpp'),'-o',str(cpu)],env=env,check=True)
report['phase_host_test']=subprocess.check_output([str(cpu)],env=env,text=True)
cuda=Path(build['installation']['cuda_home']);subprocess.run([str(cuda/'bin/nvcc'),'-O3','-std=c++11','-ccbin',str(host),'-arch=sm_120','-U_GNU_SOURCE','-D_DEFAULT_SOURCE','-D_XOPEN_SOURCE=700','-include',str(native_source/'src/cuda/CudaHostCompat.h'),str(native_source/'tests/test_gpu_phase_reduction.cu'),'-o',str(gpu)],env=env,check=True)
report['phase_device_test']=subprocess.check_output([str(gpu)],env=env,text=True)
report['passed']=True;b.write(out/'results.json',report);print('ALL',len(report['cases']),'SHAPE/SIZE AND INDEPENDENT PHASE CHECKS PASSED',flush=True)
