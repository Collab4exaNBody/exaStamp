#!/usr/bin/env python3
"""Verify exaStamp's compute_descriptor_k2b descriptor + compact derivative aggregate.

k2b has no external (e.g. LAMMPS) reference implementation to cross-validate against -- it's a
paper-specific potential (Dezaphie et al. 2025), not a LAMMPS pair style. Instead:

  1. Independent re-implementation: recompute D_i[k] = sum_j g_k(r_ij) and the compact per-atom
     derivative aggregate directly in numpy from the same closed-form Gaussian-RBF math in
     potential.h, applied to the same periodic 8-atom system, and diff against exaStamp's
     write_descriptor_k2b output (PBC-aware, minimum-image convention, box read from Ta_small.xyz).
  2. Finite-difference cross-check: perturb one atom coordinate by +-eps in a fresh exaStamp run
     (no compute_derivative) and compare the central difference of D_i[k] against the analytic
     derivative aggregate component for that atom/DOF -- this exercises the actual exaStamp
     descriptor-only code path, independent of the numpy reimplementation.

Usage: python3 compare_k2b.py
"""
import subprocess
import sys
import numpy as np

N_RBF, R_MIN, R_CUT, SIGMA = 8, 0.5, 6.0, 0.4
BOX = 16.0
EXASTAMP = "/home/lafourcadep/local/exaStampGPU/bin/exaStamp"
DERIV_MSP = "compute_derivative_k2b_small.msp"
DERIV_TXT = "exastamp_k2b_derivative_small.txt"

def read_positions_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0])
    pos = []
    for line in lines[2:2+n]:
        parts = line.split()
        pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(pos)

def read_exastamp(path):
    ids, xyz, desc, deriv = [], [], [], []
    with open(path) as f:
        for line in f:
            vals = line.split()
            if not vals: continue
            ids.append(int(float(vals[0])))
            xyz.append([float(v) for v in vals[1:4]])
            desc.append([float(v) for v in vals[4:4+N_RBF]])
            deriv.append([float(v) for v in vals[4+N_RBF:4+N_RBF+N_RBF*3]])
    return np.array(ids), np.array(xyz), np.array(desc), np.array(deriv)

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

def main():
    pos_ref = read_positions_xyz("Ta_small.xyz")
    ids, xyz, desc_exa, deriv_exa = read_exastamp(DERIV_TXT)

    # match exaStamp rows to reference positions by rounded coordinate (id order isn't guaranteed)
    order = []
    for row_xyz in xyz:
        d2 = np.sum((pos_ref - row_xyz)**2, axis=1)
        order.append(np.argmin(d2))
    order = np.array(order)
    pos_ordered = pos_ref[order]
    assert np.allclose(pos_ordered, xyz, atol=1e-6), "position matching failed"

    D_ref, deriv_ref = k2b_descriptor_and_deriv(pos_ref, BOX)
    D_ref, deriv_ref = D_ref[order], deriv_ref[order]

    max_abs_d = np.max(np.abs(desc_exa - D_ref))
    max_abs_g = np.max(np.abs(deriv_exa - deriv_ref))
    print(f"[1] independent re-implementation vs exaStamp:")
    print(f"    descriptor  max abs err = {max_abs_d:.3e}")
    print(f"    derivative  max abs err = {max_abs_g:.3e}")
    ok1 = max_abs_d < 1e-9 and max_abs_g < 1e-9

    # --- [2] finite-difference cross-check on the real exaStamp descriptor-only code path ---
    eps = 1e-5
    atom_idx, axis = 1, 0   # perturb atom #1's x coordinate
    header = open("Ta_small.xyz").readlines()[:2]
    lines = open("Ta_small.xyz").readlines()[2:]

    def write_perturbed(path, delta):
        p = pos_ref.copy()
        p[atom_idx, axis] += delta
        with open(path, "w") as f:
            f.writelines(header)
            for k in range(len(lines)):
                f.write(f"Ta {p[k,0]:.15e} {p[k,1]:.15e} {p[k,2]:.15e}\n")

    write_perturbed("Ta_small_plus.xyz", eps)
    write_perturbed("Ta_small_minus.xyz", -eps)

    def run_descriptor_only(xyzfile, outtxt):
        msp = f"""species:
  - Ta: {{ mass: 180.95 Da, z: 73, charge: 0 e- }}

setup_system:
  - domain:
      cell_size: 8.0 ang
      periodic: [ true, true, true ]
      expandable: false
  - read_xyz_file_with_xform:
      filename: "{xyzfile}"

global:
  max_iteration: 1
  dt: 1.0e-3 ps
  rcut_max: 6.0 ang

simulation_epilog:
  - compute_descriptor_k2b:
      rcut: 6.0 ang
      parameters: {{ n_rbf: 8, r_min: 0.5, r_cut: 6.0, sigma: 0.4, delta: 1.5, w: [ -1.057, -2.0949, 0.90561, -2.56538, 0.21529, -0.80587, -2.65201, 0.04461 ] }}
  - write_descriptor_k2b:
      filename: "{outtxt}"
      fields: [ id, x, y, z, descriptor ]
"""
        mspfile = outtxt.replace(".txt", ".msp")
        with open(mspfile, "w") as f: f.write(msp)
        r = subprocess.run([EXASTAMP, mspfile], capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout[-3000:]); print(r.stderr[-3000:])
            sys.exit(1)

    run_descriptor_only("Ta_small_plus.xyz", "fd_plus.txt")
    run_descriptor_only("Ta_small_minus.xyz", "fd_minus.txt")

    def read_desc_only(path):
        ids, xyz, desc = [], [], []
        with open(path) as f:
            for line in f:
                vals = line.split()
                if not vals: continue
                ids.append(int(float(vals[0])))
                xyz.append([float(v) for v in vals[1:4]])
                desc.append([float(v) for v in vals[4:4+N_RBF]])
        return np.array(xyz), np.array(desc)

    xyz_p, D_p = read_desc_only("fd_plus.txt")
    xyz_m, D_m = read_desc_only("fd_minus.txt")

    # match rows to original unperturbed order via nearest to pos_ref (perturbation is tiny)
    def match(xyz_rows):
        idx = []
        for row in xyz_rows:
            d2 = np.sum((pos_ref - row)**2, axis=1)
            idx.append(np.argmin(d2))
        return np.array(idx)

    op = match(xyz_p); om = match(xyz_m)
    D_p, D_m = D_p[op], D_m[om]

    fd = (D_p - D_m) / (2*eps)   # [natoms, K], d(D_i[k]) / d(pos[atom_idx, axis])

    # analytic: deriv_exa[atom_idx's row in original ordering, k*3+axis] gives d(D_atom_idx[k])/d(pos[atom_idx,axis])...
    # but we need d(D_i[k])/d(pos[atom_idx,axis]) for EVERY atom i, which is exactly column
    # (k*3+axis) of deriv_exa read at row = atom_idx (since deriv_agg[a][k*3+xyz] is defined as
    # sum over i of d(D_i[k])/d(r_a) -- see compute_descriptor_k2b's doc string).
    orig_row_for_atom = np.where(order == atom_idx)[0][0]
    analytic = deriv_exa[orig_row_for_atom, axis::3]   # k*3+axis for k=0..K-1 -> stride 3 from offset axis

    # sum over all atoms i of fd[i,k] should equal analytic[k] (fd gives d(D_i[k])/dpos per atom,
    # deriv_agg gives the SAME quantity already summed over i by construction)
    fd_summed = fd.sum(axis=0)
    max_abs_fd = np.max(np.abs(fd_summed - analytic))
    print(f"[2] finite-difference (atom {atom_idx}, axis {axis}) vs analytic aggregate:")
    print(f"    max abs err = {max_abs_fd:.3e}  (eps={eps})")
    ok2 = max_abs_fd < 1e-5

    if ok1 and ok2:
        print("PASS")
    else:
        print("FAIL")
        sys.exit(1)

if __name__ == "__main__":
    main()
