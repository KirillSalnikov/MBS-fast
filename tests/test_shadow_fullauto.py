#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser(description="Verify shadow-model propagation and reject mixed-model FullAuto resume")
parser.add_argument('binary', type=Path)
parser.add_argument('--out', required=True, type=Path)
options = parser.parse_args()
root = options.out.resolve()
root.mkdir(parents=True, exist_ok=True)
binary = str(options.binary.resolve())
base = ['--method','po','--backend','cpu','--particle','1','.2','.1','--ri','1.3116','0',
        '--wavelength-um','1','--max-reflections','4','--threads','1','--fullauto','.1',
        '--fullauto-resources','off','--fullauto-cutoffs','off','--cutoff-profile','off',
        '--fullauto-mirror','off','--fullauto-rotation-symmetry','off','--fullauto-analytic-controls','off',
        '--phi-points','1','--scattering-grid','0','5','1','16','--fullauto-pilot','32',
        '--fullauto-initial','64','--fullauto-min-correction','32','--fullauto-theta-start','5',
        '--max-orientations','16384','--convergence-passes','1','--close']
reports = {}
for model in ['projected', 'surface']:
    target = root / model
    cmd = [binary] + base + ['--output', str(target)]
    if model == 'surface':
        cmd += ['--po-shadow-model', model]
    with (root / (model+'.log')).open('w') as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True,
                       env={**os.environ, 'OMP_NUM_THREADS':'1'}, timeout=600)
    # Inspect real worker recipes, not only the controller's requested options.
    commands = []
    for file in target.glob('jobs/j*/command.args'):
        args = file.read_text().splitlines()
        assert '--po-shadow-model' in args, file
        assert args[args.index('--po-shadow-model')+1] == model, file
        commands.append(str(file))
    assert commands, (model, 'worker command did not record the shadow model')
    output = target / 'mueller_fullauto.dat'
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    changed = cmd[:]
    if '--po-shadow-model' in changed:
        changed[changed.index('--po-shadow-model')+1] = 'projected'
    else:
        changed += ['--po-shadow-model', 'surface']
    with (root / (model+'_mismatch.log')).open('w') as log:
        result = subprocess.run(changed, stdout=log, stderr=subprocess.STDOUT,
                                env={**os.environ, 'OMP_NUM_THREADS':'1'}, timeout=60)
    message = (root / (model+'_mismatch.log')).read_text()
    assert result.returncode != 0 and 'resume inputs/binary/environment changed' in message, message[-3000:]
    assert hashlib.sha256(output.read_bytes()).hexdigest() == digest
    reports[model] = {'worker_recipe_examples': commands[:3], 'changed_model_resume_rejected': True,
                      'completed_output_preserved': True}
(root/'verification.json').write_text(json.dumps(reports, indent=2)+'\n')
print(json.dumps(reports, indent=2))
