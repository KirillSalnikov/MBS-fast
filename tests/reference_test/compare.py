#!/usr/bin/env python3
"""Compare MBS-fast output with reference Test_Ii_PO data."""
import numpy as np
import sys, os

def load_ref(path):
    """Load reference file: Name L theta CSA M11..M44"""
    theta, m11 = [], []
    with open(path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().replace('"','').split()
            if len(parts) >= 5:
                theta.append(float(parts[2]))
                m11.append(float(parts[4]))
    return np.array(theta), np.array(m11)

def load_mbs(path):
    """Load MBS-fast output: theta 2pi*dcos M11..M44"""
    d = np.loadtxt(path, skiprows=1)
    return d[:,0], d[:,2]

if len(sys.argv) < 4:
    print("Usage: python compare.py ref_file.txt mbs_file.dat label")
    sys.exit(1)

ref_file, mbs_file, label = sys.argv[1], sys.argv[2], sys.argv[3]
rt, rm = load_ref(ref_file)
mt, mm = load_mbs(mbs_file)

print(f"=== {label} ===")
print(f"Reference: {len(rt)} points, MBS: {len(mt)} points")

# Interpolate MBS onto reference theta grid
from numpy import interp
mm_interp = interp(rt, mt, mm)

# Compare
mask = rm > 0
rdiff = np.abs(mm_interp[mask] - rm[mask]) / rm[mask] * 100

print(f"Relative difference: mean={np.mean(rdiff):.2f}%, max={np.max(rdiff):.2f}%, median={np.median(rdiff):.2f}%")
print(f"M11 at key angles:")
for angle in [0, 1, 5, 10, 22, 160, 170, 175, 180]:
    ir = np.argmin(np.abs(rt - angle))
    im = np.argmin(np.abs(mt - angle))
    if abs(rt[ir] - angle) < 1 and abs(mt[im] - angle) < 1:
        rd = abs(mm[im] - rm[ir]) / max(abs(rm[ir]), 1e-30) * 100
        print(f"  theta={angle:>3}: ref={rm[ir]:.4e}  mbs={mm[im]:.4e}  diff={rd:.2f}%")
