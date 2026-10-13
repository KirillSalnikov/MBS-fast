#!/usr/bin/env python3
"""GO fullauto: independent native audit, exact prefixes and durable recovery.

Keeps all commands, logs and matrices. Requires NumPy; never touches GPU jobs.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
T95 = 2.3646242510102993


def read_matrix(path):
    return np.loadtxt(path, skiprows=1, ndmin=2)[:, 2:]


def run(command, log, expected=0):
    Path(str(log)+'.args').write_text('\n'.join(map(str, command))+'\n')
    started = time.monotonic()
    with Path(log).open('w') as output:
        result = subprocess.run(list(map(str, command)), stdout=output,
                                stderr=subprocess.STDOUT, timeout=600)
    assert result.returncode == expected, (result.returncode, Path(log).read_text()[-3000:])
    return time.monotonic()-started


def history(folder):
    return [json.loads(line) for line in (folder/'fullauto_history.jsonl').read_text().splitlines()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=ROOT/'bin/mbs_po')
    parser.add_argument('--old-binary', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    binary, output = args.binary.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    core = [binary, '--method', 'go', '--backend', 'cpu', '--particle', '1', '.2', '.1',
            '--ri', '1.31', '0', '--wavelength-um', '1', '--max-reflections', '2',
            '--threads', '1', '--close']
    native = core+['--scattering-grid', '0', '180', '1', '6']
    report = {}
    # The first two dimensions must retain legacy seeded GO exactly.
    run(native+['--sobol-seed', '256', '11', '--output', output/'native'], output/'native.log')
    if args.old_binary:
        run([args.old_binary.resolve()]+native[1:]+['--sobol-seed', '256', '11',
            '--output', output/'legacy'], output/'legacy.log')
        for ending in ('_all.dat', '_forward.dat', '_back.dat'):
            a = np.loadtxt(next((output/'native').glob('*'+ending)), skiprows=1)
            b = np.loadtxt(next((output/'legacy').glob('*'+ending)), skiprows=1)
            np.testing.assert_array_equal(a, b)
        report['legacy_GO_unchanged'] = True
    # Offset windows must assemble the full high-precision Mueller product.
    whole = output/'whole'
    run(native+['--symmetry', '1', '1', '--go-fullauto-worker', '--sobol-seed', '256', '11',
        '--output', whole], output/'whole.log')
    parts = []
    for offset in range(0, 256, 64):
        part = output/('offset'+str(offset))
        run(native+['--symmetry', '1', '1', '--go-fullauto-worker', '--sobol-seed', '64', '11',
            '--sobol-offset', str(offset), '--output', part], output/('offset'+str(offset)+'.log'))
        parts.append(read_matrix(part/(part.name+'.dat')))
    expected = read_matrix(whole/'whole.dat')
    difference = np.max(np.abs(np.mean(parts, axis=0)-expected))/np.max(expected[:, 0])
    assert difference < 1e-12, difference
    report['exact_windows_max_all16_scaled_global_M11'] = float(difference)
    # Recover a genuine 1% calculation after its first matrix snapshot.
    folder = output/'one_percent_recovery'
    command = core+['--fullauto', '.01', '--fullauto-initial', '256',
        '--max-orientations', '262144', '--fullauto-max-rounds', '16',
        '--scattering-grid', '0', '180', '1', '6', '--fullauto-cpu-workers', '2',
        '--fullauto-work-queue', 'on', '--output', folder]
    Path(output/'recovery.args').write_text('\n'.join(map(str, command))+'\n')
    with (output/'interrupted.log').open('w') as log:
        proc = subprocess.Popen(list(map(str, command)), stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True)
        deadline = time.monotonic()+60
        try:
            while not (folder/'intermediate_current.json').exists():
                assert proc.poll() is None, (output/'interrupted.log').read_text()
                assert time.monotonic() < deadline
                time.sleep(.005)
            ready = {p: (p.read_bytes(), p.stat().st_mtime_ns)
                     for p in (folder/'jobs').glob('*/complete.txt')}
            assert ready
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=10)
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait(timeout=10)
    # Reap-time metadata distinguishes zombies, and old workers are stopped.
    time.sleep(.1)
    seconds = run(command, output/'resumed.log')
    for path, (contents, stamp) in ready.items():
        assert path.read_bytes() == contents and path.stat().st_mtime_ns == stamp
    status = json.loads((folder/'fullauto_status.json').read_text())
    assert status['status'] == 'converged' and status['target_depth'] == 2
    assert status['phi'] == 1 and status['interpolation_budget'] == 0
    assert not status['mirror_gamma']
    count = status['counts_per_seed'][0]
    independent = json.loads((folder/'independent_final_validation.json').read_text())
    assert independent['passed'] and independent['conservative_bound95'] <= .01
    # Fresh standalone workers use whole prefixes, rather than controller windows.
    samples = []
    for seed in status['seeds']:
        part = output/('audit_seed'+str(seed))
        run(native+['--symmetry', '1', '1', '--go-fullauto-worker', '--sobol-seed', str(count),
            str(seed), '--output', part], output/('audit_seed'+str(seed)+'.log'))
        samples.append(read_matrix(part/(part.name+'.dat')))
    samples = np.array(samples)
    mean = samples.mean(axis=0)
    actual = np.vstack([read_matrix(folder/'mueller_fullauto.dat'),
                       read_matrix(folder/'mueller_forward.dat'),
                       read_matrix(folder/'mueller_backward.dat')])
    error = np.max(np.abs(actual-mean))/np.max(mean[:, 0])
    assert error < 2e-12, error
    accuracy = np.genfromtxt(folder/'parameters_accuracy.csv', delimiter=',', names=True,
                             dtype=None, encoding='utf-8')
    for row in range(len(mean)):
        if mean[row, 0] > 0:
            ci = T95*np.std(samples[:, row], axis=0, ddof=1)/np.sqrt(8)/mean[row, 0]
            stored = np.array([accuracy[row]['estimated95_M%d%d_scaled_M11' % (j//4+1, j%4+1)]
                               for j in range(16)])
            np.testing.assert_allclose(stored, ci, rtol=1e-8, atol=1e-12)
            assert np.max(stored) <= .01
    events = history(folder)
    starts = [x['data']['folder'] for x in events if x['event'] == 'worker_start']
    assert all(starts.count(str(p.parent)) == 1 for p in ready), 'Completed work repeated'
    before = len(starts)
    run(command, output/'completed_resume.log')
    assert sum(x['event']=='worker_start' for x in history(folder)) == before
    # Corrupt a hashed final contract: verification must precede writing it.
    contract = folder/'precision_contract.json'
    saved = contract.read_bytes()
    contract.write_bytes(saved+b' ')
    run(command, output/'corrupt_final.log', expected=2)
    assert 'completed output changed' in (output/'corrupt_final.log').read_text()
    contract.write_bytes(saved)
    run(command, output/'restored_final.log')
    assert json.loads((folder/'fullauto_status.json').read_text())['status']=='converged'
    report['one_percent'] = dict(N_per_seed=count, total_orientations=8*count,
        resumed_seconds=seconds, preserved_jobs=len(ready), independent=independent,
        full_prefix_matrix_difference=float(error), completed_resume_no_new_workers=True,
        final_corruption_rejected=True)
    # No-scattering bins must not be reported as a converged 1% result.
    emptycore = list(core)
    emptycore[emptycore.index('--ri')+1] = '1'
    empty = output/'empty_bins'
    run(emptycore+['--fullauto', '.01', '--fullauto-initial', '256',
        '--max-orientations', '256', '--scattering-grid', '0', '180', '1', '6',
        '--output', empty], output/'empty_bins.log', expected=3)
    assert json.loads((empty/'fullauto_status.json').read_text())['status']=='requires_further_sampling'
    report['empty_bins_not_certified'] = True
    # Other physical geometries and absorption exercise native GO, with no controls.
    spec = importlib.util.spec_from_file_location('shapes', ROOT/'tests/test_fullauto_arbitrary.py')
    shapes = importlib.util.module_from_spec(spec); spec.loader.exec_module(shapes)
    tetra = output/'asymmetric_tetra.particle';shapes.asymmetric_tetrahedron(tetra)
    cases = [('asymmetric', tetra, '1.31', '0'),
             ('concave', ROOT/'examples/particles/concave_hexagonal.particle', '1.31', '0'),
             ('absorbing', tetra, '1.45', '.02')]
    report['geometries'] = []
    for name, shape, real, imag in cases:
        result = output/name
        cmd = [binary, '--method', 'go', '--backend', 'cpu', '--particle-file', shape,
            '--ri', real, imag, '--wavelength-um', '1', '--max-reflections', '4',
            '--fullauto', '.1', '--fullauto-initial', '256', '--max-orientations', '65536',
            '--scattering-grid', '0', '180', '1', '6', '--fullauto-cpu-workers', '2',
            '--threads', '1', '--output', result, '--close']
        duration = run(cmd, output/(name+'.log'))
        state = json.loads((result/'fullauto_status.json').read_text())
        assert state['status']=='converged' and not state['mirror_gamma']
        check = json.loads((result/'independent_final_validation.json').read_text())
        assert check['passed']
        report['geometries'].append(dict(case=name, seconds=duration,
            N_per_seed=state['counts_per_seed'][0], audit=check))
    # A deliberately under-resolved reflection depth must fail rather than mutate n.
    physical = output/'depth_audit_failure'
    run(core+['--fullauto', '.1', '--fullauto-initial', '256', '--max-orientations', '65536',
        '--scattering-grid', '0', '180', '1', '6', '--fullauto-physics-audit', 'on',
        '--fullauto-cpu-workers', '2', '--trace-cutoff', '.001', '--output', physical],
        output/'depth_audit_failure.log', expected=4)
    audit = json.loads((physical/'physics_validation.json').read_text())
    assert not audit['passed'] and 'cutoffs_off_residual95_scaled_M11' in audit
    assert json.loads((physical/'fullauto_status.json').read_text())['target_depth']==2
    report['physics_audit_preserves_depth'] = audit
    (output/'test_report.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
