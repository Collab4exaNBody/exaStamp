#!/usr/bin/env python3
"""Compare exastamp_snap_global_small.txt (compute_descriptor_snap_global) against
lammps_snap_global_small.txt (LAMMPS compute snap, dumped via explicit c_ID[i][j] print
statements -- see in.snap_global_small) for the same Ta_small configuration.

LAMMPS's real compute snap carries a trailing reference-label column (energy/force/virial) that
compute_descriptor_snap_global deliberately does not -- dropped from the LAMMPS side before
comparing (with pair_style zero, those labels are 0 anyway, but the shapes must match regardless).

Unlike compare_descriptors.py / compare_derivative.py, no position-matching is needed here: rows
are already keyed identically in both outputs (row 0 = summed descriptor, row 1+3*id+xyz = that
atom's gradient row, last 6 rows = virial), so this is a direct flat numeric diff.

Usage: python3 compare_global.py [exastamp_file] [lammps_file]
"""
import sys
import numpy as np

def read_matrix(path):
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.split()
            if not vals:
                continue
            rows.append([float(v) for v in vals])
    return np.array(rows)

def main():
    exastamp_file = sys.argv[1] if len(sys.argv) > 1 else "exastamp_snap_global_small.txt"
    lammps_file   = sys.argv[2] if len(sys.argv) > 2 else "lammps_snap_global_small.txt"

    a = read_matrix(exastamp_file)
    b_full = read_matrix(lammps_file)
    b = b_full[:, :-1]  # drop LAMMPS's trailing reference-label column

    if a.shape != b.shape:
        print(f"FAIL: shape mismatch: exaStamp={a.shape} LAMMPS(minus label col)={b.shape}")
        sys.exit(1)

    diff = np.abs(a - b)
    max_abs = diff.max()
    denom = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-12)
    max_rel = (diff / denom).max()

    print(f"matrix shape           = {a.shape}")
    print(f"max abs component err  = {max_abs:.3e}")
    print(f"max rel component err  = {max_rel:.3e}")

    if max_abs > 1e-6:
        print("FAIL: snap_global arrays do not match")
        sys.exit(1)
    print("PASS: exaStamp compute_descriptor_snap_global matches LAMMPS compute snap")

if __name__ == "__main__":
    main()
