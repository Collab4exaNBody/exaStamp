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
exa = {}
for l in open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(d, "thermodynamic_state.csv")):
    t = l.split()
    if not t or t[0].startswith("#"): continue
    v = list(map(float, t))  # step time n TotE KinE PotE T Pxx Pyy Pzz Pxy Pxz Pyz ...
    # exaStamp T uses 3N dof, LAMMPS 3N-3 (COM removed): rescale to LAMMPS convention
    exa[int(v[0])] = dict(TotEng=v[3], KinEng=v[4], PotEng=v[5], Temp=v[6] * 3 * v[2] / (3 * v[2] - 3),
                          Pxx=v[7], Pyy=v[8], Pzz=v[9], Pxy=v[10], Pxz=v[11], Pyz=v[12],
                          Press=(v[7] + v[8] + v[9]) / 3)
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
for q in ["PotEng", "KinEng", "TotEng", "Temp", "Press", "Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz"]:
    s = bar if q[0] == "P" and q != "PotEng" else 1.0
    scale = max(abs(lmp[k][q] * s) for k in steps) or 1.0
    err, k = max((abs(exa[k][q] - lmp[k][q] * s) / scale, k) for k in steps)
    print(f"  {q:7s} max rel err {err:.3e} (step {k}: exa {exa[k][q]:.12g} lmp {lmp[k][q]*s:.12g})")
