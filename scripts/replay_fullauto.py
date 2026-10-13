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
args = parser.parse_args()
recipe = json.loads((args.result/'controller_command.json').read_text())
if recipe['schema'] != 1:
    parser.error('Unsupported saved command schema')
command = [recipe['executable']] + recipe['argv']
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
