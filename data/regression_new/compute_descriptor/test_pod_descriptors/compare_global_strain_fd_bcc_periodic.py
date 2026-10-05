#!/usr/bin/env python3
"""Periodic BCC Ta variant of compare_global_strain_fd.py: independent finite-difference verification
of compute_descriptor_pod_global's virial/stress descriptor rows against its own energy descriptor
(row 0), this time on a fully PERIODIC simulation cell instead of an isolated cluster.

Math (same as compare_global_strain_fd.py): for the affine transform x' = F x with
F = I + eps * e_b e_a^T (a,b in {x,y,z}), dE_desc/deps|eps=0 = -virial_row[a,b] given the row's
already-force-signed convention (F_atom = +coeff . row). For a periodic system the identity only holds
if the strain is applied to the cell vectors AND the positions together (homogeneous deformation of the
whole periodic crystal) -- which is exactly what this script does.

How the strained cell reaches exaStamp: read_xyz_file_with_xform takes the Lattice="..." rows as the
cell vectors (H, row-major), computes reduced coords s = H^-T r, stores |a|,|b|,|c|-scaled coords in
the grid and sets domain.xform() = H^T diag(|a|,|b|,|c|)^-1 (times whatever scaling compute_domain_bounds
needed), so the real positions seen by the potential are exactly H^T s = r and the periodic images sit
at the strained cell vectors. Straining both the Lattice rows and the positions with the same F
(v' = F v for every row of H and every atom) is therefore exactly the homogeneous deformation above,
including shear (non-uniform xform) components.

System: 4x4x4 conventional BCC Ta cells (128 atoms, a0 = 3.304 ang -> box 13.216 ang > 2*rcut = 10 ang
with rcut = 5.0 from Ta_param.pod), with small random displacements (fixed seed) so that none of the six
virial components vanishes by cubic symmetry -- on a perfect BCC lattice the shear rows would be ~0
and the diagonal rows identical, which would make the check much weaker.

Usage: python3 compare_global_strain_fd_bcc_periodic.py
"""
import subprocess
import sys
import numpy as np

EXASTAMP = "/local_home/lafourcadep/local/exaStamp/bin/exaStamp"
EPS = 1.0e-6
A0 = 3.304          # BCC Ta lattice constant (ang)
NCELL = 4           # conventional cells per direction
RATTLE = 0.05       # max random displacement per component (ang)
SEED = 12345

BOX = A0 * NCELL
H0 = np.diag([BOX, BOX, BOX])   # rows = cell vectors a, b, c


def make_bcc(a0, n):
    basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    cells = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)], dtype=float)
    return ((cells[:, None, :] + basis[None, :, :]).reshape(-1, 3)) * a0


rng = np.random.default_rng(SEED)
POS0 = make_bcc(A0, NCELL) + rng.uniform(-RATTLE, RATTLE, size=(2 * NCELL**3, 3))
NATOMS = len(POS0)

VOIGT = ["xx", "yy", "zz", "yz", "xz", "xy"]
VOIGT_AB = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]  # (a,b): x[:,b] += eps*x[:,a]

MSP_TEMPLATE = """species:
  - Ta: {{ mass: 180.95 Da, z: 73, charge: 0 e- }}

init_parameters:
  - species
  - pod_init:
      parameters:
        pod_file: "Ta_param.pod"
        coeff_file: "Ta_coefficients.pod"

setup_system:
  - domain:
      cell_size: {cell_size} ang
      periodic: [ true, true, true ]
      expandable: false
  - read_xyz_file_with_xform:
      filename: "{xyzfile}"
      verbose: {verbose}

global:
  max_iteration: 1
  dt: 1.0e-3 ps

simulation_epilog:
  - compute_descriptor_pod: {{ compute_derivative: true }}
  - compute_descriptor_pod_global
  - write_descriptor_pod_global: {{ filename: "{outfile}" }}
"""


def strain(m, a, b, eps):
    """Apply x' = (I + eps e_b e_a^T) x to every row of m (atom positions or cell vectors)."""
    out = m.copy()
    out[:, b] += eps * m[:, a]
    return out


def write_xyz(path, h, pos):
    with open(path, "w") as f:
        f.write(f"{len(pos)}\n")
        lat = " ".join(f"{v:.15e}" for v in h.reshape(-1))
        f.write(f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for p in pos:
            f.write(f"Ta {p[0]:.15e} {p[1]:.15e} {p[2]:.15e}\n")


def run_pod_global(h, pos, tag, verbose=False):
    xyzfile = f"strain_bcc_{tag}.xyz"
    mspfile = f"strain_bcc_{tag}.msp"
    outfile = f"strain_bcc_{tag}_global.txt"
    write_xyz(xyzfile, h, pos)
    with open(mspfile, "w") as f:
        f.write(MSP_TEMPLATE.format(cell_size=2.0 * A0, xyzfile=xyzfile, outfile=outfile,
                                    verbose="true" if verbose else "false"))
    r = subprocess.run(["mpirun", "-np", "1", "bash", EXASTAMP, mspfile], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        sys.exit(1)
    if verbose:
        with open(f"strain_bcc_{tag}.log", "w") as f:
            f.write(r.stdout)
    rows = []
    with open(outfile) as f:
        for line in f:
            vals = line.split()
            if vals:
                rows.append([float(v) for v in vals])
    return np.array(rows)


def main():
    baseline = run_pod_global(H0, POS0, "baseline", verbose=True)
    ncols = baseline.shape[1]
    grad_rows = 1 + 3 * NATOMS
    virial_row0 = grad_rows
    virial_rows = baseline[virial_row0:virial_row0 + 6]
    if virial_rows.shape != (6, ncols):
        print(f"FAIL: expected 6 virial rows x {ncols} cols, got {virial_rows.shape}")
        sys.exit(1)

    # Descriptor sums grow with N, so scale the absolute tolerance by the largest virial entry.
    tol = 1e-4 * max(1.0, np.abs(virial_rows).max())
    print(f"{NATOMS} atoms, box {BOX:.4f} ang, {ncols} descriptors, tol = {tol:.3e}")

    ok = True
    for k, ((a, b), name) in enumerate(zip(VOIGT_AB, VOIGT)):
        e_plus = run_pod_global(strain(H0, a, b, EPS), strain(POS0, a, b, EPS), f"plus_{name}")[0]
        e_minus = run_pod_global(strain(H0, a, b, -EPS), strain(POS0, a, b, -EPS), f"minus_{name}")[0]

        fd = (e_plus - e_minus) / (2.0 * EPS)
        expected = -virial_rows[k]   # row is already force-signed, see compare_global_strain_fd.py
        diff = np.abs(fd - expected)
        denom = np.maximum(np.maximum(np.abs(fd), np.abs(expected)), 1e-10)
        max_abs = diff.max()
        max_rel = (diff / denom).max()
        status = "PASS" if max_abs < tol else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"{name}: max|virial| = {np.abs(expected).max():.3e}, max abs err = {max_abs:.3e}, "
              f"max rel err = {max_rel:.3e}  [{status}]")

    if ok:
        print("PASS: compute_descriptor_pod_global's virial rows match finite-difference of the "
              "energy descriptor (row 0) under a homogeneous strain of a periodic BCC Ta cell")
    else:
        print("FAIL: virial rows do not match the finite-difference derivative of the energy descriptor")
        sys.exit(1)


if __name__ == "__main__":
    main()
