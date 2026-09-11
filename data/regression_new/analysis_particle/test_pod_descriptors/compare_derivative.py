#!/usr/bin/env python3
"""Cross-validate exaStamp's compute_descriptor_pod compute_derivative aggregate (dumped by
write_descriptor_pod) against real LAMMPS compute podd/atom, on a small test system.

LAMMPS's podd/atom is NOT a compact per-atom aggregate like SNAP's snad/atom -- it's a dense
per-local-atom-row matrix, size_peratom_cols = 3*natoms*Mdesc*nClusters (see
ML-POD/compute_podd_atom.cpp), with row i's nonzero blocks at column-atom i ("+= bdd") and
column-atom = each of i's neighbors ("-= bdd"). exaStamp's PodDescriptorOp instead reduces this
straight into a compact per-atom slot: agg[a] = sum over every row i of pod[i][:, a]. This script
performs that same reduction on LAMMPS's dense dump and compares the result to exaStamp's output.

Only practical on small systems (dense matrix scales as natoms^2) -- see
data/regression_new/analysis_particle/test_pod_descriptors/Ta_small.xyz/.lmp and
in.pod_derivative_small.

Usage: python3 compare_derivative.py [exastamp_file] [lammps_file]
"""
import sys
import numpy as np

def read_exastamp(path):
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.split()
            if not vals:
                continue
            pos = tuple(float(v) for v in vals[1:4])
            agg = np.array([float(v) for v in vals[4:]])  # m*3+xyz order
            rows.append((pos, agg))
    return rows

def read_lammps_dense(path):
    """Returns (positions[natoms,3], dense[natoms(rows), 3*natoms*Mdesc(cols)])."""
    with open(path) as f:
        lines = f.readlines()
    i = 0
    rows = []
    positions = []
    while i < len(lines):
        if lines[i].startswith("ITEM: ATOMS"):
            i += 1
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                vals = lines[i].split()
                if vals:
                    positions.append([float(v) for v in vals[1:4]])
                    rows.append([float(v) for v in vals[4:]])
                i += 1
        else:
            i += 1
    return np.array(positions), np.array(rows)

def main():
    exastamp_file = sys.argv[1] if len(sys.argv) > 1 else "exastamp_pod_derivative_small.txt"
    lammps_file   = sys.argv[2] if len(sys.argv) > 2 else "lammps_pod_derivative_small.txt"

    a = read_exastamp(exastamp_file)
    positions, dense = read_lammps_dense(lammps_file)
    natoms = positions.shape[0]
    ncols = dense.shape[1]
    if ncols % natoms != 0:
        print(f"FAIL: LAMMPS dense matrix column count {ncols} not divisible by natoms {natoms}")
        sys.exit(1)
    block = ncols // natoms   # 3*Mdesc
    if block % 3 != 0:
        print(f"FAIL: per-atom block width {block} not divisible by 3")
        sys.exit(1)
    mdesc = block // 3

    if len(a) != natoms:
        print(f"FAIL: atom count mismatch: exaStamp={len(a)} LAMMPS={natoms}")
        sys.exit(1)
    nc = len(a[0][1])
    if nc != mdesc*3:
        print(f"FAIL: derivative length mismatch: exaStamp={nc} LAMMPS(3*Mdesc)={mdesc*3}")
        sys.exit(1)

    # reduce LAMMPS's dense (row=atom_i, col=atom_a*3*Mdesc + xyz*Mdesc + m) matrix: sum down the
    # rows (over every atom i) to get, per column-atom a, a (xyz,m) block -- exactly what
    # exaStamp's PodDescriptorOp scatter accumulates directly into each atom's own slot.
    reshaped = dense.reshape(natoms, natoms, 3, mdesc)   # [row_i, col_a, xyz, m]
    agg_lammps = reshaped.sum(axis=0)                    # [col_a, xyz, m]
    agg_lammps = agg_lammps.transpose(0, 2, 1).reshape(natoms, mdesc*3)  # [atom, m*3+xyz]

    # match rows by rounded position (robust to any ordering/id-offset difference)
    lut = {}
    for idx in range(natoms):
        key = tuple(round(c, 5) for c in positions[idx])
        lut[key] = agg_lammps[idx]

    max_abs = 0.0
    max_rel = 0.0
    unmatched = 0
    for pos, agg_a in a:
        key = tuple(round(c, 5) for c in pos)
        agg_b = lut.get(key)
        if agg_b is None:
            unmatched += 1
            continue
        diff = np.abs(agg_a - agg_b)
        max_abs = max(max_abs, diff.max())
        denom = np.maximum(np.maximum(np.abs(agg_a), np.abs(agg_b)), 1e-12)
        max_rel = max(max_rel, (diff / denom).max())

    print(f"atoms compared        = {len(a) - unmatched} / {len(a)}")
    print(f"unmatched positions   = {unmatched}")
    print(f"max abs component err = {max_abs:.3e}")
    print(f"max rel component err = {max_rel:.3e}")

    if unmatched > 0 or max_abs > 1e-6:
        print("FAIL: derivative aggregate does not match LAMMPS compute podd/atom")
        sys.exit(1)
    print("PASS: exaStamp compute_descriptor_pod compute_derivative matches LAMMPS compute podd/atom")

if __name__ == "__main__":
    main()
