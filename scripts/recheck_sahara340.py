#!/usr/bin/env python3
"""Recheck saved Sahara340 statistics into a separate output directory."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from sahara340_statistics import fieller


def recheck(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output == source or source in output.parents:
        raise ValueError('Output must be separate from the source ensemble tree')
    reports = []
    for status_file in sorted(source.glob('*/status.json')):
        data_file = status_file.with_name('replicas.npz')
        if not data_file.exists():
            continue
        old = json.loads(status_file.read_text())
        with np.load(data_file) as archive:
            data = archive['replicas']
        if data.shape != (8, 181, 16) or not np.isfinite(data).all():
            raise ValueError('Invalid ensemble replicas: '+str(data_file))
        mean = data.mean(axis=0)
        if (mean[:, 0] <= 0).any():
            raise ValueError('Nonpositive ensemble M11: '+str(data_file))
        raw = float(np.max(2.364624251*data.std(axis=0, ddof=1)/np.sqrt(8)/mean[:, :1]))
        ratio = float(np.max(fieller(data[:, :, 1:], data[:, :, :1])))
        ldr = float(fieller(data[:, -1, 0]-data[:, -1, 5], data[:, -1, 0]+data[:, -1, 5]))
        policy = old.get('stopping_policy', {})
        eps = policy.get('epsilon')
        change = old.get('refinement_change')
        initial_ok = change is None and policy.get('allow_initial_stop', False)
        passed = (eps is not None and max(raw, ratio, ldr) <= eps and
                  (initial_ok or change is not None and change <= eps))
        report = dict(job=status_file.parent.name, original_status=old['status'],
                      max_pointwise95_scaled_M11=raw, max_pointwise95_ratio=ratio if np.isfinite(ratio) else None,
                      ldr180_halfwidth95=ldr if np.isfinite(ldr) else None,
                      refinement_change=change, tolerance=eps, numerical_criteria_passed=passed,
                      source_sha256=hashlib.sha256(data_file.read_bytes()).hexdigest(),
                      controller_status_changed=False, solver_work_repeated=False)
        reports.append(report)
    output.mkdir(parents=True, exist_ok=True)
    (output/'statistics_recheck.json').write_text(json.dumps(reports, indent=2, allow_nan=False)+'\n')
    return reports


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(recheck(args.source, args.output), indent=2))
