#!/usr/bin/env python3
"""Step-by-step diff of exaStamp thermodynamic_state.csv vs LAMMPS log.lammps (thermo_style of in.lammps).
usage: compare.py <case_dir> [exaStamp_csv]   -> prints max relative error per quantity"""
import sys, os
d = sys.argv[1]
lmp, cols = {}, None
for l in open(os.path.join(d, "log.lammps")):
    t = l.split()
    if t[:1] == ["Step"]: cols = t; continue
    if cols and len(t) == len(cols):
        try: lmp[int(t[0])] = dict(zip(cols, map(float, t)))
        except ValueError: cols = None
    elif cols: cols = None
exa, names = {}, None
KNOWN = ["Step", "Time (ps)", "Particles", "Temp. (K)", "Tot. E. (eV/part)", "Pot. E. (eV/part)", "Kin. E. (eV/part)",
         "Press. (Pa)", "Pxx (Pa)", "Pyy (Pa)", "Pzz (Pa)", "Pxy (Pa)", "Pxz (Pa)", "Pyz (Pa)", "sMises (Pa)",
         "Vol. (ang^3)", "Mass"]  # "Mv/Ext/Imb." has no value
for l in open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(d, "thermodynamic_state.csv")):
    if l.split()[:1] == ["Step"]:  # header names contain spaces and may be 1 space apart: locate known names
        names = [n for _, n in sorted((l.find(n), n) for n in KNOWN if n in l)]; continue
    t = l.split()
    if not t or not t[0].isdigit(): continue
    c = dict(zip(names, map(float, t)))
    r = dict(TotEng=c["Tot. E. (eV/part)"], KinEng=c["Kin. E. (eV/part)"], PotEng=c["Pot. E. (eV/part)"],
             Temp=c["Temp. (K)"])   # both codes use 3N-3 dof (COM removed)
    if "Pxx (Pa)" in c:
        for q in ["Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz"]: r[q] = c[q + " (Pa)"]
        r["Press"] = (r["Pxx"] + r["Pyy"] + r["Pzz"]) / 3
    elif "Press. (Pa)" in c: r["Press"] = c["Press. (Pa)"]
    exa[int(t[0])] = r
bar = 1e5  # LAMMPS metal pressure unit (bar) -> Pa
# LAMMPS 'units metal' hardcodes rounded constants (update.cpp); exaStamp uses CODATA 2018.
# Map LAMMPS thermo onto exact constants so only genuine discrepancies remain:
#   KE  scales with mvv2e, T with mvv2e/boltz, scalar P = nktv2p*(2*KE_tot + tr W)/(3V).
eV, amu, kB = 1.602176634e-19, 1.66053906892e-27, 1.380649e-23
mvv2e_ex, boltz_ex, nktv2p_ex = amu * 1e4 / eV, kB / eV, eV * 1e30 / 1e5
mvv2e_l, boltz_l, nktv2p_l = 1.0364269e-4, 8.617343e-5, 1.6021765e6
xyz = open(os.path.join(d, "init.xyz")).read().splitlines()
N = int(xyz[0]); lat = xyz[1].split('"')[1].split(); V = float(lat[0]) * float(lat[4]) * float(lat[8])
rk = mvv2e_ex / mvv2e_l
for r in lmp.values():
    trW = 3 * r["Press"] * V / nktv2p_l - 2 * N * r["KinEng"]      # eV, virial trace
    r["KinEng"] *= rk
    r["Temp"] *= rk * boltz_l / boltz_ex
    r["Press"] = nktv2p_ex * (2 * N * r["KinEng"] + trW) / (3 * V)
    for q in ["Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz"]:
        r[q] *= nktv2p_ex / nktv2p_l   # kinetic share of components still carries mvv2e (<=6.4e-8)
steps = sorted(set(lmp) & set(exa))
print(f"{d}: {len(steps)} common steps ({steps[0]}..{steps[-1]})")
for q in [q for q in ["PotEng", "KinEng", "TotEng", "Temp", "Press", "Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz"] if q in exa[steps[0]]]:
    s = bar if q[0] == "P" and q != "PotEng" else 1.0
    scale = max(abs(lmp[k][q] * s) for k in steps) or 1.0
    err, k = max((abs(exa[k][q] - lmp[k][q] * s) / scale, k) for k in steps)
    print(f"  {q:7s} max rel err {err:.3e} (step {k}: exa {exa[k][q]:.12g} lmp {lmp[k][q]*s:.12g})")
