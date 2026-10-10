#!/usr/bin/env python3
"""Step-by-step diff of exaStamp thermodynamic_state.csv vs LAMMPS log.lammps (snap_ttm cases).
usage: compare.py <case_dir> [exaStamp_csv]   -> prints max relative error per quantity
Step 0 is skipped for constant_te: in LAMMPS, f_myte (fix ave/atom) is still 0 at setup, so the step-0 energy is the Te=0 one."""
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
KNOWN = ["Step", "Time (ps)", "Particles", "Temp. (K)", "Tot. E. (eV/part)", "Pot. E. (eV/part)", "Kin. E. (eV/part)",
         "Press. (Pa)", "sMises (Pa)", "Vol. (ang^3)", "Mass", "Elec. E. (eV)", "Ion Transf. E. (eV)"]  # "Mv/Ext/Imb." has no value
exa = {}
names = None
xyz = open(os.path.join(d, "init.xyz")).read().splitlines()
N = int(xyz[0]); lat = xyz[1].split('"')[1].split(); V = float(lat[0]) * float(lat[4]) * float(lat[8])
for l in open(sys.argv[2] if len(sys.argv) > 2 else os.path.join(d, "thermodynamic_state.csv")):
    if l.split()[:1] == ["Step"]:  # header names contain spaces and may be 1 space apart: locate known names
        names = [n for _, n in sorted((l.find(n), n) for n in KNOWN if n in l)]; continue
    t = l.split()
    if not t or not t[0].isdigit(): continue
    c = dict(zip(names, map(float, t)))
    r = dict(PotEng=c["Pot. E. (eV/part)"], KinEng=c["Kin. E. (eV/part)"],
             # exaStamp T uses 3N dof, LAMMPS 3N-3 (COM removed): rescale to LAMMPS convention
             Temp=c["Temp. (K)"] * 3 * N / (3 * N - 3))
    # exaStamp TotE includes the electronic energy when TTM is on, LAMMPS etotal does not
    r["TotEng"] = c["Tot. E. (eV/part)"] - c.get("Elec. E. (eV)", 0.0) / N
    if "Press. (Pa)" in c: r["Press"] = c["Press. (Pa)"]
    exa[int(t[0])] = r
bar = 1e5
# LAMMPS 'units metal' rounded constants -> CODATA 2018 (same mapping as ../snap/compare.py)
eV, amu, kB = 1.602176634e-19, 1.66053906892e-27, 1.380649e-23
mvv2e_ex, boltz_ex, nktv2p_ex = amu * 1e4 / eV, kB / eV, eV * 1e30 / 1e5
mvv2e_l, boltz_l, nktv2p_l = 1.0364269e-4, 8.617343e-5, 1.6021765e6
rk = mvv2e_ex / mvv2e_l
for r in lmp.values():
    trW = 3 * r["Press"] * V / nktv2p_l - 2 * N * r["KinEng"]
    r["KinEng"] *= rk
    r["Temp"] *= rk * boltz_l / boltz_ex
    r["Press"] = nktv2p_ex * (2 * N * r["KinEng"] + trW) / (3 * V)
first = 0 if "c_avgte" in next(iter(lmp.values())) else 1   # ttm_coupled: fix ttm pre_force sets Te at step 0 too
steps = sorted(k for k in set(lmp) & set(exa) if k >= first)
print(f"{d}: {len(steps)} common steps ({steps[0]}..{steps[-1]})")
for q in [q for q in ["PotEng", "KinEng", "TotEng", "Temp", "Press"] if q in exa[steps[0]]]:
    s = bar if q == "Press" else 1.0
    scale = max(abs(lmp[k][q] * s) for k in steps) or 1.0
    err, k = max((abs(exa[k][q] - lmp[k][q] * s) / scale, k) for k in steps)
    print(f"  {q:7s} max rel err {err:.3e} (step {k}: exa {exa[k][q]:.12g} lmp {lmp[k][q]*s:.12g})")
