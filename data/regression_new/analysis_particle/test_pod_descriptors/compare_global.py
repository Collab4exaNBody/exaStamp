#!/usr/bin/env python3
"""Compare exastamp_pod_global_small.txt (compute_descriptor_pod_global) against
lammps_pod_global_small.txt (LAMMPS compute pod/global, dumped via explicit c_ID[i][j] print
statements -- see in.pod_global_small) for the same Ta_small configuration.

exaStamp's array now has 6 extra virial rows LAMMPS's own pod/global doesn't have (POD fitting there
just doesn't use stress, not a structural limitation -- added to match compute_descriptor_snap_global's
shape). Rows 0..3*natoms-1 are compared directly against LAMMPS as before. The virial rows have no
LAMMPS ground truth, so instead they get a self-consistency check: recomputed independently in this
script from Ta_small.xyz's own known atom positions and exaStamp's own (LAMMPS-verified) gradient
rows, using the same Voigt formula [xx,yy,zz,yz,xz,xy] = sum_atom pos . gradient_row(atom).

Usage: python3 compare_global.py [exastamp_file] [lammps_file]
"""
import sys
import numpy as np

NATOMS = 8

def read_matrix(path):
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.split()
            if not vals:
                continue
            rows.append([float(v) for v in vals])
    return np.array(rows)

def read_positions(xyz_path):
    with open(xyz_path) as f:
        lines = f.readlines()
    n = int(lines[0])
    pos = np.zeros((n, 3))
    for i in range(n):
        toks = lines[2 + i].split()
        pos[i] = [float(toks[1]), float(toks[2]), float(toks[3])]
    return pos

def main():
    exastamp_file = sys.argv[1] if len(sys.argv) > 1 else "exastamp_pod_global_small.txt"
    lammps_file   = sys.argv[2] if len(sys.argv) > 2 else "lammps_pod_global_small.txt"
    xyz_file      = sys.argv[3] if len(sys.argv) > 3 else "Ta_small.xyz"

    a_full = read_matrix(exastamp_file)
    b = read_matrix(lammps_file)

    n_grad_rows = 1 + 3 * NATOMS
    a = a_full[:n_grad_rows]

    if a.shape != b.shape:
        print(f"FAIL: shape mismatch (rows 0..{n_grad_rows-1}): exaStamp={a.shape} LAMMPS={b.shape}")
        sys.exit(1)

    diff = np.abs(a - b)
    max_abs = diff.max()
    denom = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-12)
    max_rel = (diff / denom).max()

    print(f"descriptor+gradient rows shape = {a.shape}")
    print(f"max abs component err          = {max_abs:.3e}")
    print(f"max rel component err          = {max_rel:.3e}")

    ok = True
    if max_abs > 1e-6:
        print("FAIL: pod_global descriptor+gradient rows do not match LAMMPS")
        ok = False
    else:
        print("PASS: descriptor+gradient rows match LAMMPS compute pod/global")

    # Virial rows: self-consistency check against exaStamp's own (LAMMPS-verified) gradient rows.
    virial_rows = a_full[n_grad_rows:n_grad_rows + 6]
    if virial_rows.shape[0] != 6:
        print(f"FAIL: expected 6 virial rows, got {virial_rows.shape[0]}")
        sys.exit(1)

    pos = read_positions(xyz_file)
    ncols = a_full.shape[1]
    expected_virial = np.zeros((6, ncols))
    for atom in range(NATOMS):
        grad = a_full[1 + 3*atom : 1 + 3*atom + 3]  # dx,dy,dz rows for this atom
        rx, ry, rz = pos[atom]
        expected_virial[0] += grad[0] * rx  # xx
        expected_virial[1] += grad[1] * ry  # yy
        expected_virial[2] += grad[2] * rz  # zz
        expected_virial[3] += grad[2] * ry  # yz
        expected_virial[4] += grad[2] * rx  # xz
        expected_virial[5] += grad[1] * rx  # xy

    vdiff = np.abs(virial_rows - expected_virial)
    vmax_abs = vdiff.max()
    print(f"virial rows max abs err (self-consistency) = {vmax_abs:.3e}")
    if vmax_abs > 1e-9:
        print("FAIL: virial rows inconsistent with exaStamp's own gradient rows")
        ok = False
    else:
        print("PASS: virial rows match independent recomputation from exaStamp's own gradient rows")

    if not ok:
        sys.exit(1)
    print("PASS: exaStamp compute_descriptor_pod_global matches LAMMPS compute pod/global (+ self-consistent virial rows)")

if __name__ == "__main__":
    main()
