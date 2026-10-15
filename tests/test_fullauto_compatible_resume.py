#!/usr/bin/env python3
"""Check the explicit 8ae9e85 memory-fix compatibility bridge on real data."""
import argparse, json, os, subprocess
from pathlib import Path
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--old-binary', required=True, type=Path)
parser.add_argument('--binary', required=True, type=Path)
parser.add_argument('--output', required=True, type=Path)
a = parser.parse_args()
root = a.output.resolve(); root.mkdir(parents=True, exist_ok=True)
old, new = str(a.old_binary.resolve()), str(a.binary.resolve())
env = dict(os.environ)
# A short complete real calculation exercises exact resume and all guards.
target=root/'calculation'
base=['--method','po','--backend','cpu','--particle','1','.2','.1','--ri','1.3116','0','--wavelength-um','1','--max-reflections','4','--threads','1','--fullauto','.1','--fullauto-resources','off','--fullauto-cutoffs','off','--cutoff-profile','off','--fullauto-mirror','off','--fullauto-rotation-symmetry','off','--fullauto-analytic-controls','off','--phi-points','1','--scattering-grid','0','5','1','16','--fullauto-pilot','32','--fullauto-initial','64','--fullauto-min-correction','32','--fullauto-theta-start','5','--max-orientations','16384','--convergence-passes','1','--output',str(target),'--close']
def run(name,exe,arguments):
 with (root/(name+'.log')).open('w') as f:
  result=subprocess.run([exe]+arguments,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=600)
 return result.returncode,(root/(name+'.log')).read_text()
code,log=run('resume_old',old,base);assert code==0,log[-2000:]
files={str(p):p.read_bytes() for p in target.glob('jobs/*/complete.txt')}
original=(target/'mueller_fullauto.dat').read_bytes()
code,log=run('resume_wrong_binary',new,base);assert code!=0 and 'resume inputs/binary/environment changed' in log,log[-2000:]
compat=base+['--fullauto-resume-compatible',old]
code,log=run('resume_compatible',new,compat);assert code==0 and 'no calculation repeated' in log,log[-2000:]
assert all(Path(p).read_bytes()==v for p,v in files.items())
assert (target/'mueller_fullauto.dat').read_bytes()==original
assert (target/'controller_resume_command.json').exists()
changed=compat[:];changed[changed.index('--fullauto')+1]='.2'
code,log=run('resume_changed_accuracy',new,changed);assert code!=0 and 'resume inputs/binary/environment changed' in log,log[-2000:]
changed=compat[:];changed[changed.index('--fullauto-resume-compatible')+1]=new
code,log=run('resume_wrong_origin',new,changed);assert code!=0 and 'only from releases 8ae9e85 and 89f83a3' in log,log[-2000:]
changed=compat[:];changed[changed.index('--output')+1]=str(root/'no_checkpoint')
code,log=run('resume_missing_checkpoint',new,changed);assert code!=0 and 'requires an existing fullauto result' in log,log[-2000:]
# A version-report fixture must be rejected before any checkpoint is changed.
wrong_profile = root/'wrong_profile.py'
version = subprocess.check_output([old, '--version'], env=env).decode().replace('backend-build:', 'backend-build: incompatible')
wrong_profile.write_text('#!/usr/bin/env python3\nprint('+repr(version)+', end="")\n')
wrong_profile.chmod(0o755)
changed=compat[:];changed[changed.index('--fullauto-resume-compatible')+1]=str(wrong_profile)
code,log=run('resume_wrong_profile',new,changed)
assert code!=0 and 'identical CUDA architecture' in log,log[-2000:]
# The saved resume recipe must use the fixed executable on subsequent runs.
replay = Path(__file__).resolve().parents[1]/'scripts/replay_fullauto.py'
with (root/'resume_saved_recipe.log').open('w') as f:
 subprocess.run(['python3',str(replay),str(target)],env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=600)
assert 'no calculation repeated' in (root/'resume_saved_recipe.log').read_text()
assert (target/'mueller_fullauto.dat').read_bytes()==original
(root/'resume_validation.json').write_text(json.dumps({'completed_job_markers_preserved':len(files),'matrix_unchanged':True,'ordinary_binary_mismatch_rejected':True,'changed_accuracy_rejected':True,'wrong_origin_rejected':True,'missing_checkpoint_rejected':True,'wrong_profile_rejected':True,'saved_recipe_reuses_fixed_binary':True},indent=2))
print('resume tests passed',flush=True)
