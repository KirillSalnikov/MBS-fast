#!/usr/bin/env python3
"""Resume a fullauto result with its exact saved executable, CLI and environment."""
import argparse
import json
import os
from pathlib import Path
import shlex

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('result', type=Path, help='Fullauto output directory')
parser.add_argument('--show', action='store_true', help='Print the saved command without running')
parser.add_argument('--compatible-executable', type=Path,
                    help='Use the analytic CUDA memory fix with an existing 8ae9e85 result; the binary verifies compatibility')
args = parser.parse_args()
resume = args.result/'controller_resume_command.json'
recipe = json.loads((resume if resume.exists() and not args.compatible_executable
                     else args.result/'controller_command.json').read_text())
if recipe['schema'] != 1:
    parser.error('Unsupported saved command schema')
command = [recipe['executable']] + recipe['argv']
if args.compatible_executable:
    command += ['--fullauto-resume-compatible', recipe['executable']]
    recipe['executable'] = str(args.compatible_executable.resolve())
    command[0] = recipe['executable']
if args.show:
    print('cd '+shlex.quote(recipe['cwd']))
    print(' '.join(shlex.quote(argument) for argument in command))
else:
    environment = os.environ.copy()
    for key in list(environment):
        if key.startswith('MBS_') or key in ('CUDA_VISIBLE_DEVICES', 'LD_LIBRARY_PATH', 'OMP_NUM_THREADS'):
            environment.pop(key)
    environment.update(recipe['environment'])
    os.chdir(recipe['cwd'])
    os.execve(recipe['executable'], command, environment)
