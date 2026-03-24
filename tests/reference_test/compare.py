#!/usr/bin/env python3
"""Compare MBS output with reference Test_Ii_PO data or old vs new.

Usage:
  python3 compare.py <ref_file.txt> <mbs_file.dat>            # MBS vs reference
  python3 compare.py old_vs_new <old.dat> <new.dat>            # Old vs New
"""
import numpy as np
import sys

def load_ref(path):
    """Load reference: Name L theta CSA M11..M44 (16 Mueller)"""
    data = []
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().replace('"','').split()
            if len(parts) >= 20:
                theta = float(parts[2])
                mueller = [float(x) for x in parts[4:20]]
                data.append([theta] + mueller)
    return np.array(data)

def load_mbs(path):
    """Load MBS phi-averaged: ScAngle 2pi*dcos M11..M44"""
    data = []
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 18:
                theta = float(parts[0])
                mueller = [float(x) for x in parts[2:18]]
                data.append([theta] + mueller)
    return np.array(data)

def compare(d1, d2, l1="D1", l2="D2"):
    t1, t2 = d1[:,0], d2[:,0]

    # Interpolate d2 onto d1 grid if needed
    if len(t1) != len(t2) or not np.allclose(t1, t2, atol=0.01):
        from scipy.interpolate import interp1d
        cols = []
        for c in range(1, d2.shape[1]):
            f = interp1d(t2, d2[:,c], kind='linear', bounds_error=False, fill_value=0)
            cols.append(f(t1))
        d2 = np.column_stack([t1] + cols)

    m11a, m11b = d1[:,1], d2[:,1]
    mask = np.abs(m11a) > 1e-6

    rel = np.zeros_like(m11a)
    rel[mask] = np.abs(m11b[mask] - m11a[mask]) / np.abs(m11a[mask])

    print(f"\n{'='*70}")
    print(f"{l1} vs {l2}")
    print(f"{'='*70}")
    print(f"Points: {len(t1)}, range [{t1[0]:.4f}, {t1[-1]:.4f}]°")
    print()

    # M11 table
    step = max(1, len(t1) // 30)
    print(f"{'theta':>8s}  {'M11_'+l1:>14s}  {'M11_'+l2:>14s}  {'ratio':>10s}  {'rel_err':>10s}")
    for i in range(0, len(t1), step):
        r = m11b[i]/m11a[i] if abs(m11a[i]) > 1 else 0
        print(f"{t1[i]:8.3f}  {m11a[i]:14.4f}  {m11b[i]:14.4f}  {r:10.6f}  {rel[i]:10.2e}")

    print(f"\nM11 stats (theta > 0.1°):")
    far = t1 > 0.1
    m = mask & far
    if np.any(m):
        print(f"  Max rel diff:  {np.max(rel[m]):.4e} ({np.max(rel[m])*100:.3f}%)")
        print(f"  Mean rel diff: {np.mean(rel[m]):.4e}")
        print(f"  Median:        {np.median(rel[m]):.4e}")

    # Key angle comparison
    print(f"\nKey angles:")
    for angle in [0, 0.5, 1, 2, 5, 10, 15, 20, 22, 25]:
        idx = np.argmin(np.abs(t1 - angle))
        if abs(t1[idx] - angle) < 1.5:
            r = m11b[idx]/m11a[idx] if abs(m11a[idx]) > 1 else 0
            print(f"  θ={angle:5.1f}°: {l1}={m11a[idx]:.6e}  {l2}={m11b[idx]:.6e}  ratio={r:.6f}")

    # Normalized elements
    print(f"\nNormalized |Mij/M11| max abs diff (θ>1°):")
    names = ['M12','M13','M14','M21','M22','M23','M24','M31','M32','M33','M34','M41','M42','M43','M44']
    far2 = t1 > 1.0
    for c, name in enumerate(names, start=2):
        n1 = d1[far2,c] / np.maximum(np.abs(d1[far2,1]), 1e-20)
        n2 = d2[far2,c] / np.maximum(np.abs(d2[far2,1]), 1e-20)
        m = np.abs(n1) > 1e-6
        if np.any(m):
            maxd = np.max(np.abs(n2[m] - n1[m]))
            if maxd > 1e-4:
                print(f"  {name}/M11: {maxd:.4e}")

    return rel

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    if sys.argv[1] == "old_vs_new":
        d1 = load_mbs(sys.argv[2])
        d2 = load_mbs(sys.argv[3])
        compare(d1, d2, "Old", "New")
    else:
        ref = load_ref(sys.argv[1])
        mbs = load_mbs(sys.argv[2])
        compare(ref, mbs, "Ref", "MBS")
