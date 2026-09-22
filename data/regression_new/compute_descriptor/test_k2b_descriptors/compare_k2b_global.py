#!/usr/bin/env python3
"""Verify exaStamp's compute_descriptor_k2b_global global design matrix.

k2b has no external (e.g. LAMMPS) reference implementation -- it's a paper-specific potential
(Dezaphie et al. 2025), not a LAMMPS pair style. Instead, this builds an independent numpy
reimplementation (k2b_descriptor_and_deriv: D_i[k] = sum_j g_k(r_ij), minimum-image PBC, same
closed-form Gaussian-RBF math as potential.h/k2b_descriptor_op.h -- originally shared with
compare_k2b.py's own per-atom check before that script was removed, inlined here directly now that
this is the only surviving k2b regression check) and builds the expected global matrix from it:

  row 0             = sum_i D_i[k]
  rows 1..3*natoms  = deriv_agg[atom_id, k*3+xyz]  (deriv_agg is already the per-atom aggregate
                       compute_k2b's own compute_derivative computes, by construction)
  rows 3N+1..3N+6   = sum_atom pos[atom] . deriv_agg[atom], Voigt order [xx,yy,zz,yz,xz,xy]

No position-matching needed: field::id is assigned in file order by read_xyz_file_with_xform
(0-indexed, contiguous), so Ta_small.xyz's row i IS atom id i directly.

Usage: python3 compare_k2b_global.py
"""
import sys
import numpy as np

N_RBF, R_MIN, R_CUT, SIGMA = 8, 0.5, 6.0, 0.4
BOX = 16.0
GLOBAL_TXT = "exastamp_k2b_global_small.txt"

def read_positions_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0])
    pos = []
    for line in lines[2:2+n]:
        parts = line.split()
        pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(pos)

def k2b_descriptor_and_deriv(pos, box):
    """Reference numpy implementation, minimum-image PBC, same math as potential.h /
    k2b_descriptor_op.h. Returns D[natoms,K] and deriv_agg[natoms,K*3] (k*3+xyz order)."""
    n = pos.shape[0]
    s = R_MIN + np.arange(N_RBF) * ((R_CUT - R_MIN) / (N_RBF - 1))
    inv_2sig2 = 0.5 / SIGMA**2
    inv_sig2 = 1.0 / SIGMA**2
    D = np.zeros((n, N_RBF))
    deriv = np.zeros((n, N_RBF, 3))
    for i in range(n):
        for j in range(n):
            if i == j: continue
            dr = pos[j] - pos[i]                       # matches exaStamp's buf.drx convention
            dr -= box * np.round(dr / box)              # minimum image
            r = np.linalg.norm(dr)
            if r <= 0.0 or r > box/2: continue
            for k in range(N_RBF):
                d = r - s[k]
                arg = d*d*inv_2sig2
                if arg > 20.0: continue
                g = np.exp(-arg)
                D[i, k] += g
                dgdr = -d * inv_sig2 * g
                c = -dgdr / r
                v = c * dr
                deriv[i, k, :] += v
                deriv[j, k, :] -= v
    return D, deriv.reshape(n, N_RBF*3)

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
    pos = read_positions_xyz("Ta_small.xyz")
    natoms = pos.shape[0]

    D, deriv_agg = k2b_descriptor_and_deriv(pos, BOX)  # D[natoms,K], deriv_agg[natoms,K*3] (k*3+xyz)

    expected = np.zeros((1 + 3*natoms + 6, N_RBF))
    expected[0] = D.sum(axis=0)
    for a in range(natoms):
        grad = deriv_agg[a].reshape(N_RBF, 3).T  # [3,K]: grad[0]=dx, grad[1]=dy, grad[2]=dz for atom a
        expected[1 + 3*a + 0] = grad[0]
        expected[1 + 3*a + 1] = grad[1]
        expected[1 + 3*a + 2] = grad[2]

    virial_row0 = 1 + 3*natoms
    for a in range(natoms):
        dx, dy, dz = expected[1+3*a+0], expected[1+3*a+1], expected[1+3*a+2]
        rx, ry, rz = pos[a]
        expected[virial_row0+0] += dx*rx  # xx
        expected[virial_row0+1] += dy*ry  # yy
        expected[virial_row0+2] += dz*rz  # zz
        expected[virial_row0+3] += dz*ry  # yz
        expected[virial_row0+4] += dz*rx  # xz
        expected[virial_row0+5] += dy*rx  # xy

    actual = read_matrix(GLOBAL_TXT)

    if actual.shape != expected.shape:
        print(f"FAIL: shape mismatch: exaStamp={actual.shape} expected={expected.shape}")
        sys.exit(1)

    diff = np.abs(actual - expected)
    max_abs = diff.max()
    print(f"matrix shape           = {actual.shape}")
    print(f"max abs component err  = {max_abs:.3e}")

    if max_abs > 1e-9:
        print("FAIL: compute_descriptor_k2b_global does not match independent reimplementation")
        sys.exit(1)
    print("PASS: exaStamp compute_descriptor_k2b_global matches independent numpy reimplementation")

if __name__ == "__main__":
    main()
