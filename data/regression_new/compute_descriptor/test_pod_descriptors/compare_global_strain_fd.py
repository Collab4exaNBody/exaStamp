#!/usr/bin/env python3
"""Independent finite-difference verification of compute_descriptor_pod_global's virial/stress
descriptor rows against its own energy descriptor (row 0) -- LAMMPS's compute pod/global has no
stress output to cross-validate against (see compare_global.py's own note), so this is the only
real ground truth available for those 6 rows; compare_global.py's existing "virial self-consistency"
check only recomputes the same formula from the same gradient rows, it is not independent.

Math: for E_desc(r_1..r_N) (any component of the row-0 descriptor vector) and the affine transform
r_i' = r_i + eps * r_i[a] * e_b (a,b in {x,y,z}, scaling only the eps=0 component's Voigt pairing),
the chain rule gives dE_desc/deps|eps=0 = sum_i (dE_desc/dr_i)[b] * r_i[a] = -virial_row[a,b] as
compute_descriptor_pod_global.cu codes it. That leading minus is a REAL FINDING, not a test bug: the
central+=/neighbor-= scatter in pod_descriptor_op.h (rij = r_neighbor - r_central) actually stores
-dE/dr (already force-signed) into the pda_* aggregate, not +dE/dr as the old docstrings claimed --
confirmed by an independent by-hand chain-rule derivation of the same scatter, matching this FD
result to ~1e-6/1e-7 (the eps-limited precision), not just approximately. Fixed in
compute_descriptor_pod_global.cu's (and snap/k2b/mtp's) documentation: F_atom = +coeff . row, not
-coeff . row. So this script checks the finite difference against -virial_row (dE_desc/deps ==
-virial_row is the CORRECT relation given the row's already-force-signed convention), and this is
EXACTLY the virial-row formula compute_descriptor_pod_global.cu codes directly (see its Voigt-order
comments: xx=(a=x,b=x), yy=(y,y), zz=(z,z), yz=(a=y,b=z), xz=(a=x,b=z), xy=(a=x,b=y)). This holds for
ANY fixed absolute coordinate origin (it's a chain-rule identity, not a statement about mechanical
equilibrium), so central-differencing E_desc(+eps)/E_desc(-eps) under this exact transform and
comparing to the baseline run's own virial row is a real, independent check of the code.

Uses a fresh, deliberately NON-PERIODIC 8-atom cluster (same relative geometry as Ta_small.xyz,
recentered) instead of the periodic Ta_small.xyz + compute_pod_global_small.msp system: that system's
actual simulated box is 8x8x8 ang (from the .msp's own domain cell_size, NOT the file's 16x16x16
Lattice header -- verified directly, Hbis reduces to identity there since domain_xform ends up
uniform-scale 1.0), with rcut=5.0 ang from Ta_param.pod, so atoms genuinely wrap and interact with
periodic images. Reusing it here would need the strain applied consistently to domain.xform() AND
positions together (read_xyz_file_with_xform's own uniform_scale/Hbis logic, non-trivial to get
right for an anisotropic/shear strain) -- going non-periodic sidesteps all of that: domain.xform()
stays identity always, so a pure position-only affine transform in Python is exactly correct, with
zero box/PBC ambiguity. Verified via a verbose dry run that the resulting domain is genuinely
identity-xform, bounds (0,0,0)-(30,30,30) -- no unexpected rescaling.

Usage: python3 compare_global_strain_fd.py
"""
import subprocess
import sys
import numpy as np

EXASTAMP = "/home/lafourcadep/local/exaStampGPU/bin/exaStamp"
EPS = 1.0e-6
BOX = 30.0

# Same 8 Ta atoms as Ta_small.xyz, recentered from their original centroid (8,8,8) to (15,15,15)
# so they sit comfortably inside a big non-periodic [0,BOX]^3 box with margin on every side.
_RAW = np.array([
    [8.1750456936, 8.0675146880, 7.9586966329],
    [6.5749544236, 6.6824013740, 9.4028288149],
    [6.5627856886, 6.7272519940, 6.2177479229],
    [9.8281017016, 6.4342601570, 9.3882153169],
    [9.7309075366, 6.5807854330, 6.1069658509],
    [6.6697696516, 10.007944756, 9.3914108119],
    [9.8428670096, 9.8600903200, 9.4579364889],
    [6.6155682946, 9.6397512780, 6.0761981609],
])
POS0 = _RAW - _RAW.mean(axis=0) + np.array([15.0, 15.0, 15.0])
NATOMS = len(POS0)

VOIGT = ["xx", "yy", "zz", "yz", "xz", "xy"]
VOIGT_AB = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]  # (a,b): perturb pos[:,b] += eps*pos[:,a]

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
      cell_size: {box} ang
      periodic: [ false, false, false ]
      expandable: false
  - read_xyz_file_with_xform:
      filename: "{xyzfile}"
      verbose: false

global:
  max_iteration: 1
  dt: 1.0e-3 ps

simulation_epilog:
  - compute_descriptor_pod: {{ compute_derivative: true }}
  - compute_descriptor_pod_global
  - write_descriptor_pod_global: {{ filename: "{outfile}" }}
"""


def write_xyz(path, pos):
    with open(path, "w") as f:
        f.write(f"{len(pos)}\n")
        f.write(f'Lattice="{BOX} 0.0 0.0 0.0 {BOX} 0.0 0.0 0.0 {BOX}"\n')
        for p in pos:
            f.write(f"Ta {p[0]:.15e} {p[1]:.15e} {p[2]:.15e}\n")


def run_pod_global(pos, tag):
    xyzfile = f"strain_{tag}.xyz"
    mspfile = f"strain_{tag}.msp"
    outfile = f"strain_{tag}_global.txt"
    write_xyz(xyzfile, pos)
    with open(mspfile, "w") as f:
        f.write(MSP_TEMPLATE.format(box=BOX, xyzfile=xyzfile, outfile=outfile))
    r = subprocess.run(["mpirun", "-np", "1", "bash", EXASTAMP, mspfile], capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:])
        print(r.stderr[-3000:])
        sys.exit(1)
    rows = []
    with open(outfile) as f:
        for line in f:
            vals = line.split()
            if vals:
                rows.append([float(v) for v in vals])
    return np.array(rows)


def main():
    baseline = run_pod_global(POS0, "baseline")
    ncols = baseline.shape[1]
    grad_rows = 1 + 3 * NATOMS
    virial_row0 = grad_rows
    virial_rows = baseline[virial_row0:virial_row0 + 6]
    if virial_rows.shape != (6, ncols):
        print(f"FAIL: expected 6 virial rows x {ncols} cols, got {virial_rows.shape}")
        sys.exit(1)

    ok = True
    fd_matrix = np.zeros((6, ncols))
    for k, ((a, b), name) in enumerate(zip(VOIGT_AB, VOIGT)):
        pos_plus = POS0.copy()
        pos_plus[:, b] += EPS * POS0[:, a]
        pos_minus = POS0.copy()
        pos_minus[:, b] -= EPS * POS0[:, a]

        e_plus = run_pod_global(pos_plus, f"plus_{name}")[0]   # row 0 = energy descriptor
        e_minus = run_pod_global(pos_minus, f"minus_{name}")[0]

        fd = (e_plus - e_minus) / (2.0 * EPS)
        fd_matrix[k] = fd

        # row is already force-signed (F_atom = +coeff.row, see compute_descriptor_pod_global.cu's
        # header comment) -- so dE_desc/deps == -virial_row is the correct relation, not +virial_row.
        expected = -virial_rows[k]
        diff = np.abs(fd - expected)
        denom = np.maximum(np.maximum(np.abs(fd), np.abs(expected)), 1e-10)
        max_abs = diff.max()
        max_rel = (diff / denom).max()
        status = "PASS" if max_abs < 1e-4 else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"{name}: max abs err = {max_abs:.3e}, max rel err = {max_rel:.3e}  [{status}]")

    if ok:
        print("PASS: compute_descriptor_pod_global's virial rows match finite-difference of the "
              "energy descriptor (row 0) under an affine strain -- independent of compare_global.py's "
              "self-consistency check")
    else:
        print("FAIL: virial rows do not match the finite-difference derivative of the energy descriptor")
        sys.exit(1)


if __name__ == "__main__":
    main()
