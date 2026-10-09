#!/usr/bin/env python3
"""Finite-difference check of the gradient and virial rows of compute_descriptor_<family>_global.

The energy descriptor (row 0) is central-differenced and compared with the unstrained run:
  - moving one atom: dE_desc/dr_i = SIGN * row[1+3*i+xyz]
  - homogeneous strain x' = (I + eps e_b e_a^T) x of the cell vectors AND the positions of a periodic
    cell: dE_desc/deps = SIGN * virial_row[a,b]
SIGN = -1: the rows of every family are force-signed (F_atom = +coeff . row). The strained cell goes through
read_xyz_file_with_xform, so the domain xform is not the identity (shear included): this checks that
the virial rows use real-frame positions.

System: 4x4x4 conventional BCC cells (128 atoms, a0 = 3.304 ang, box 13.216 ang > 2*rcut for every
family below), with small random displacements so that no virial component vanishes by symmetry.

Usage: EXASTAMP=/path/to/bin/exaStamp python3 check_global_virial_fd.py {pod,snap,mtp,k2b} [...]
       (default: all four families)
"""
import os
import subprocess
import sys
import tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXASTAMP = os.environ.get("EXASTAMP", "exaStamp")
EPS = 1.0e-6
A0 = 3.304
NCELL = 4
RATTLE = 0.05
SEED = 12345
BOX = A0 * NCELL
H0 = np.diag([BOX, BOX, BOX])   # rows = cell vectors a, b, c

FAMILIES = {
    "pod": ("Ta", "Ta: {{ mass: 180.95 Da, z: 73, charge: 0 e- }}",
            """  - pod_init:
      parameters:
        pod_file: "{d}/test_pod_descriptors/Ta_param.pod"
        coeff_file: "{d}/test_pod_descriptors/Ta_coefficients.pod"
"""),
    "snap": ("Ta", "Ta: {{ mass: 180.95 Da, z: 73, charge: 0 e- }}",
             """  - snap_init:
      parameters:
        param: "{d}/test_snap_descriptors/param.txt"
        coef: "{d}/test_snap_descriptors/coeff.txt"
"""),
    "mtp": ("Cu", "Cu: {{ mass: 63.546 Da, z: 29, charge: 0 e- }}",
            """  - mtp_init:
      parameters:
        mtp_file: "{d}/test_mtp_descriptors/pot.almtp"
"""),
    "k2b": ("Ta", "Ta: {{ mass: 180.95 Da, z: 73, charge: 0 e- }}",
            """  - k2b_init:
      rcut: 6.0 ang
      parameters: {{ n_rbf: 8, r_min: 0.5, r_cut: 6.0, sigma: 0.4, delta: 1.5, w: [ -1.057, -2.0949, 0.90561, -2.56538, 0.21529, -0.80587, -2.65201, 0.04461 ] }}
"""),
}

MSP_TEMPLATE = """species:
  - {species}

init_parameters:
  - species
{init}
setup_system:
  - domain:
      cell_size: {cell_size} ang
      periodic: [ true, true, true ]
      expandable: false
  - read_xyz_file_with_xform:
      filename: "{xyzfile}"
      verbose: false

global:
  max_iteration: 1
  dt: 1.0e-3 ps

simulation_epilog:
  - compute_descriptor_{fam}: {{{{ compute_derivative: true }}}}
  - compute_descriptor_{fam}_global
  - write_descriptor_{fam}_global: {{{{ filename: "{outfile}" }}}}
"""

SIGN = {"pod": -1.0, "snap": -1.0, "mtp": -1.0, "k2b": -1.0}
GRAD_PROBES = [(0, 0), (5, 1), (77, 2)]   # (atom, component)

VOIGT = ["xx", "yy", "zz", "yz", "xz", "xy"]
VOIGT_AB = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]  # (a,b): x[:,b] += eps*x[:,a]


def make_bcc(a0, n):
    basis = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    cells = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)], dtype=float)
    return ((cells[:, None, :] + basis[None, :, :]).reshape(-1, 3)) * a0


POS0 = make_bcc(A0, NCELL) + np.random.default_rng(SEED).uniform(-RATTLE, RATTLE, size=(2 * NCELL**3, 3))
NATOMS = len(POS0)


def strain(m, a, b, eps):
    out = m.copy()
    out[:, b] += eps * m[:, a]
    return out


def run_global(workdir, fam, h, pos, tag):
    elem, species, init = FAMILIES[fam]
    xyzfile = os.path.join(workdir, f"{tag}.xyz")
    mspfile = os.path.join(workdir, f"{tag}.msp")
    outfile = os.path.join(workdir, f"{tag}_global.txt")
    with open(xyzfile, "w") as f:
        f.write(f"{len(pos)}\n")
        lat = " ".join(f"{v:.15e}" for v in h.reshape(-1))
        f.write(f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        for p in pos:
            f.write(f"{elem} {p[0]:.15e} {p[1]:.15e} {p[2]:.15e}\n")
    msp = MSP_TEMPLATE.format(species=species, init=init, cell_size=2.0 * A0, xyzfile=xyzfile,
                              fam=fam, outfile=outfile)
    with open(mspfile, "w") as f:
        f.write(msp.format(d=HERE))
    r = subprocess.run(["mpirun", "-np", "1", EXASTAMP, mspfile], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        sys.exit(1)
    return np.loadtxt(outfile, ndmin=2)


def check(fam):
    with tempfile.TemporaryDirectory(prefix=f"fd_{fam}_") as workdir:
        base = run_global(workdir, fam, H0, POS0, "baseline")
        sign = SIGN[fam]
        ok = True
        for atom, c in GRAD_PROBES:
            pp, pm = POS0.copy(), POS0.copy()
            pp[atom, c] += EPS
            pm[atom, c] -= EPS
            fd = (run_global(workdir, fam, H0, pp, "gplus")[0] - run_global(workdir, fam, H0, pm, "gminus")[0]) / (2.0 * EPS)
            row = base[1 + 3 * atom + c]
            tol = 1e-4 * max(1.0, np.abs(row).max())
            err = np.abs(fd - sign * row).max()
            ok = ok and err < tol
            print(f"  {fam} grad atom {atom} {'xyz'[c]}: max|row| = {np.abs(row).max():.3e}, max abs err = {err:.3e} "
                  f"[{'PASS' if err < tol else 'FAIL'}]")
        virial = base[1 + 3 * NATOMS:1 + 3 * NATOMS + 6]
        tol = 1e-4 * max(1.0, np.abs(virial).max())   # descriptor sums grow with N
        for k, ((a, b), name) in enumerate(zip(VOIGT_AB, VOIGT)):
            e_plus = run_global(workdir, fam, strain(H0, a, b, EPS), strain(POS0, a, b, EPS), f"plus_{name}")[0]
            e_minus = run_global(workdir, fam, strain(H0, a, b, -EPS), strain(POS0, a, b, -EPS), f"minus_{name}")[0]
            fd = (e_plus - e_minus) / (2.0 * EPS)
            err = np.abs(fd - sign * virial[k]).max()
            ok = ok and err < tol
            print(f"  {fam} {name}: max|virial| = {np.abs(virial[k]).max():.3e}, max abs err = {err:.3e} "
                  f"[{'PASS' if err < tol else 'FAIL'}]")
    return ok


def main():
    fams = sys.argv[1:] or list(FAMILIES)
    results = {fam: check(fam) for fam in fams}
    for fam, ok in results.items():
        print(f"{fam}: {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if all(results.values()) else 1)


if __name__ == "__main__":
    main()
