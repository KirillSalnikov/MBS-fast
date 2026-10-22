#!/usr/bin/env python3
"""Test isolated angular continuation using frozen real server inventories.

The inventory directory must contain gpu1/dima/epyc1 requested_theta.csv,
adaptive_checkpoint.tsv, execution_tuning.tsv, automatic_depth.tsv,
depth_selection.json, configuration.txt, fullauto_parameters.json,
polarization_weights.tsv, and means_shared.cache. The input is never modified.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inventory', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=True)
objects = sorted(p for p in (root / 'cpu/build/obj').rglob('*.o') if p.name not in ['main.o', 'FullAuto.o'])
assert objects, 'Build the isolated cpu target first'
probe = output / 'fullauto_window_probe'
command = ['/usr/bin/mpicxx', '-std=c++11', '-O2', '-fopenmp', '-DUSE_MPI', '-I' + str(root / 'src'),
           '-I' + str(root / 'src/cuda'), str(root / 'tests/fullauto_window_probe.cpp'),
           *map(str, objects), '-lm', '-lgomp', '-o', str(probe)]
subprocess.run(command, check=True)
summary = {}
for host in ['gpu1', 'dima', 'epyc1']:
    parent = output / host / 'parent'
    assert not parent.exists(), 'Use a fresh output directory to retain test evidence'
    parent.mkdir(parents=True)
    for file in (args.inventory / host).iterdir():
        if file.is_file():
            shutil.copy2(file, parent / file.name)
    (parent / 'jobs').mkdir()
    params = json.loads((parent / 'fullauto_parameters.json').read_text())
    before = {str(f.relative_to(parent)): hashlib.sha256(f.read_bytes()).hexdigest()
              for f in parent.rglob('*') if f.is_file()}
    result = subprocess.run([str(probe), str(parent), str(params['initial_count'])],
                            capture_output=True, text=True)
    after = {name: hashlib.sha256((parent / name).read_bytes()).hexdigest() for name in before}
    summary[host] = {'returncode': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr,
                     'original_parent_files_checked': len(before), 'all_original_files_byte_identical': before == after}
    if result.returncode == 0:
        summary[host]['result'] = json.loads(result.stdout.strip().splitlines()[-1])
        summary[host]['readonly_cache'] = json.loads((parent.parent / 'readonly_cache/result.json').read_text())
    print(host, json.dumps(summary[host]), flush=True)
    (output / 'results.json').write_text(json.dumps(summary, indent=2) + '\n')
assert all(row['returncode'] == 0 and row['all_original_files_byte_identical'] for row in summary.values()), summary
print('PASS: three real checkpoint imports/resume; exact N/phi/coefficients; 5%/1.25%; parent immutable; eight-seed partial cache reuse')
