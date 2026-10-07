#!/usr/bin/env python3
"""Validation of the coulombic operators (wolf, dsf, ewald, pppm) against LAMMPS, on 12000 atoms of disturbed UO2.

Only inputs are stored here. First generate the reference (LAMMPS with KSPACE, metal units), from this folder :
    lmp -in in.coul -var case <case> -log log_<case>.lammps    # writes dump_<case>.{0,10}.txt
    case = wolf|dsf|ewald_auto|ewald_fixed|ewald_tri|ewald_tri_fixed|pppm|pppm_tri
exaStamp runs (from this folder, 1 thread is enough) :
    exaStamp exastamp_<case>.msp                 # e.g. exastamp_wolf.msp, exastamp_ewald_fixed+sym.msp
    mpirun -np 2 exaStamp exastamp_<case>.msp    # MPI check, same files
Comparison (needs numpy) :
    python compare.py <case> [<case> ...]       # e.g. wolf dsf ewald_auto ewald_fixed wolf+sym

  - thermo : potential energy (eV) and pressure tensor (bar) at every step 0..10
  - step 10: per-atom positions, forces (eV/ang) and energies (eV), matched by id (exaStamp id + 1 = LAMMPS id)
Variants <case>+<tag> are compared with the LAMMPS run of <case>.
Expected (2026-10) : |dPE| ~1e-5 eV (csv precision), dP ~9e-8 relative (LAMMPS nktv2p constant), dx ~3e-9 ang,
dF ~3e-8 eV/ang, dE_atom ~3e-8 eV : differences are at output precision.
Cases ewald_tri / ewald_tri_fixed / pppm_tri use a triclinic cell (UO2_tri.lmp / UO2_tri_ext.xyz, tilts xy=5 xz=3 yz=4 ang).
Cases pppm / pppm_tri (coulombic_pppm, LAMMPS kspace_style pppm 1e-5 order 5 diff ik) : g_ewald and mesh are chosen
automatically by both codes and agree (pppm 0.33415064, 72x72x72 ; pppm_tri 0.33026023, 72x75x75) ; same expected
differences as Ewald. Variant pppm_tri+species reads species charges.
Case pppm_ad (kspace_modify diff ad, orthogonal box only as in LAMMPS ; in.coul switches the box to ortho) : same mesh
(96x96x96), but the automatic g_ewald differs (0.32925083 vs LAMMPS 0.32923668) : LAMMPS's ad estimate is a difference
of nearly equal sums and g_ewald comes from a finite difference derivative of it, so a 1e-13 relative change of the sum
of q^2 already moves g_ewald by 3e-5 ; both meet the accuracy. Hence |dPE| ~1e-3 eV, dF ~9e-8 eV/ang. Variant pppm_ad+g
uses LAMMPS's g_ewald and gives the usual agreement (dF 3.2e-8 eV/ang, dE_atom 2.9e-8 eV, dP 9e-8). Variants pppm+gpu / pppm_tri+gpu run with
nogpu: false (GPU kernels + cuFFT on Cuda builds) and give the same agreement.
Cases pppm_slab / pppm_slab_auto / pppm_tri_slab (slab correction EW3DC, z non periodic, kspace_modify slab 3.0 /
slab auto / slab 3.0 with an xy tilt) use UO2_slab*.lmp / UO2_slab*_ext.xyz (made by mk_slab.py) : g_ewald, mesh,
auto volfactor (2.5130678) and estimated accuracy are those of LAMMPS ; dF ~1.5e-7 eV/ang for max|F| ~15 eV/ang
(surface atoms), dE_atom ~7e-7 eV, dP ~1.5e-7 ; step 0 energies agree to 3e-6 eV, the trajectories then drift apart
(dx ~1.4e-8 ang at step 10). Variant pppm_slab+ad (diff ad + slab) is compared with LAMMPS diff ik + slab : agreement
at the accuracy level only (dF 3.4e-4 x,y 7.6e-4 z ; LAMMPS's own ik/ad difference is ~3e-4). LAMMPS diff ad + slab
(case pppm_ad_slab) gives wrong z forces (fieldforce_ad uses nz/zprd instead of nz/zprd_slab, |dF_z| up to 13 eV/ang
vs ik), exaStamp uses the extended spacing : at fixed g_ewald, ik and ad converge to the same forces when the mesh is
refined (fine mesh difference 1.6e-5 eV/ang in all directions).
Variants +sym use symmetric pair computation (use_symmetry). Variants +fold (ewald_fixed, wolf, pppm) add
ghost_fold_back : pairs from owned cells only on half neighbor lists, ghost contributions folded back by
update_virial_force_energy_from_ghost (CPU, 1 thread, ewald real space 140 -> 89 ms/step vs LAMMPS 69 ms) ; same agreement.
Variants +pair use the pair potential template front-ends (coul_wolf_pair, coul_dsf) with species charges. For dsf+pair,
energies differ by 3.74e-2 eV in total (8e-6 eV per atom) : the template subtracts e(rcut) per pair, which is not exactly
0 for DSF (A&S erfc) ; forces are identical.
"""
import sys
import numpy as np

EV_INTERNAL = 1.602176634e-19 / (1.66053906892e-27 * 1e-20 / 1e-24)  # eV in exaStamp internal energy units
PA_PER_BAR = 1.0e5


def lammps_thermo(case):
    rows, on = [], False
    for line in open(f"log_{case}.lammps"):
        w = line.split()
        if w[:2] == ["Step", "PotEng"]:
            on = True
            continue
        if on:
            if not w or not w[0].isdigit():
                break
            rows.append([float(x) for x in w])
    a = np.array(rows)
    # step pe ke etotal press pxx pyy pzz pxy pxz pyz
    return {int(r[0]): (r[1], r[5:11]) for r in a}


def exastamp_thermo(case):
    out = {}
    for line in open(f"thermo_exastamp_{case}.csv"):
        if line.startswith("#"):
            continue
        w = line.split()
        n = float(w[2])
        out[int(w[0])] = (float(w[5]) * n, np.array([float(x) for x in w[7:13]]) / PA_PER_BAR)
    return out


def lammps_dump(fname):
    lines = open(fname).read().splitlines()
    i = lines.index(next(l for l in lines if l.startswith("ITEM: ATOMS")))
    cols = lines[i].split()[2:]
    a = np.array([[float(x) for x in l.split()] for l in lines[i + 1:]])
    a = a[np.argsort(a[:, cols.index("id")])]
    get = lambda *k: a[:, [cols.index(c) for c in k]]
    return get("x", "y", "z"), get("fx", "fy", "fz"), get("c_pe")[:, 0]


def exastamp_xyz(fname):
    lines = open(fname).read().splitlines()
    props = next(t for t in lines[1].split() if t.startswith("Properties=")).split("=")[1].split(":")
    # name:type:count triplets -> column offsets
    off, cols = 0, {}
    for k in range(0, len(props), 3):
        cols[props[k]] = (off, int(props[k + 2]))
        off += int(props[k + 2])
    rows = [l.split() for l in lines[2:]]
    num = lambda name: np.array([[float(r[cols[name][0] + j]) for j in range(cols[name][1])] for r in rows])
    ids = num("id")[:, 0].astype(int)
    o = np.argsort(ids)
    return ids[o], num("pos")[o], num("force")[o], num("ep")[o, 0], num("type")[o, 0].astype(int)


def compare(case):
    print(f"===== {case} =====")
    base = case.split("+")[0]  # variants (e.g. wolf+sym) compare with the LAMMPS run of their base case
    lt, et = lammps_thermo(base), exastamp_thermo(case)
    worst_e, worst_p = 0.0, 0.0
    for s in sorted(set(lt) & set(et)):
        de = abs(et[s][0] - lt[s][0])
        dp = np.max(np.abs(et[s][1] - lt[s][1]) / np.maximum(np.abs(lt[s][1]), 1.0))
        worst_e, worst_p = max(worst_e, de), max(worst_p, dp)
        if s in (0, 10):
            print(f"step {s:2d}: PE exaStamp={et[s][0]:.10e} LAMMPS={lt[s][0]:.10e} |dE|={de:.3e} eV ; "
                  f"max rel dP={dp:.3e} (Pxx {et[s][1][0]:.6e} vs {lt[s][1][0]:.6e} bar)")
    print(f"steps 0..10 : max |dPE| = {worst_e:.3e} eV , max rel dP = {worst_p:.3e}")

    lx, lf, le = lammps_dump(f"dump_{base}.10.txt")
    ids, ex, ef, ee, et = exastamp_xyz(f"exastamp_{case}_000000010.xyz")
    assert np.array_equal(ids + 1, np.arange(1, len(lx) + 1)), "id mismatch"
    # minimum image in the (possibly triclinic) cell read from the exaStamp snapshot (Lattice rows = cell vectors)
    lat = open(f"exastamp_{case}_000000010.xyz").read().splitlines()[1].split('"')[1].split()
    H = np.array([float(v) for v in lat]).reshape(3, 3).T
    sd = np.linalg.solve(H, (ex - lx).T).T
    sd -= np.round(sd)
    dx = (H @ sd.T).T
    # at snapshot time the force field holds accelerations (F/m, internal units) : F[eV/ang] = a * m / EV_INTERNAL
    mass = np.array([15.999, 238.02891])[et]  # species order O, U
    fscale = (mass / EV_INTERNAL)[:, None]
    escale = 1.0 if np.max(np.abs(ee)) < 100 * np.max(np.abs(le)) else 1.0 / EV_INTERNAL
    df = np.abs(ef * fscale - lf)
    print(f"step 10: max|dx| = {np.max(np.abs(dx)):.3e} ang , max|dF| = {np.max(df):.3e} eV/ang (max|F| {np.max(np.abs(lf)):.3e})"
          f" , max|dE_atom| = {np.max(np.abs(ee * escale - le)):.3e} eV"
          )


for c in sys.argv[1:]:
    compare(c)
