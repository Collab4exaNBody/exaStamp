#!/usr/bin/env python3
"""Compare write_descriptor_snap's single, MPI-gathered exastamp_combined.txt (fields:
[id, x, y, z, descriptor, derivative]) against LAMMPS's per-atom sna/atom + snad/atom
references, regardless of how many MPI ranks produced the exaStamp file. This is the
real multi-rank validation -- unlike exastamp_descriptors.txt/exastamp_snad.txt (each
rank overwrites the same file, so those only mean anything on a single-rank run),
exastamp_combined.txt is already the full, globally-gathered result.

Usage: python3 compare_combined.py
"""
import sys

NCOEFF = 30

def read_combined(path):
    rows = []
    with open(path) as f:
        for line in f:
            v = line.split()
            if not v:
                continue
            pos = tuple(float(x) for x in v[1:4])
            desc = [float(x) for x in v[4:4+NCOEFF]]
            deriv = [float(x) for x in v[4+NCOEFF:]]  # our [k*3+xyz] layout
            rows.append((pos, desc, deriv))
    return rows

def read_lammps_dump(path):
    rows = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: ATOMS"):
            i += 1
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                v = lines[i].split()
                if v:
                    pos = tuple(float(x) for x in v[1:4])
                    vec = [float(x) for x in v[4:]]
                    rows.append((pos, vec))
                i += 1
        else:
            i += 1
    return rows

def main():
    combined = read_combined("exastamp_combined.txt")
    lammps_desc = { tuple(round(c,6) for c in pos): vec for pos, vec in read_lammps_dump("lammps_descriptors.txt") }
    lammps_snad = { tuple(round(c,6) for c in pos): vec for pos, vec in read_lammps_dump("lammps_snad.txt") }

    if len(combined) != len(lammps_desc):
        print(f"FAIL: atom count mismatch: exaStamp={len(combined)} LAMMPS={len(lammps_desc)}")
        sys.exit(1)

    max_abs_desc, max_abs_deriv = 0.0, 0.0
    unmatched = 0
    for pos, desc, deriv in combined:
        key = tuple(round(c,6) for c in pos)
        d_ref = lammps_desc.get(key)
        s_ref = lammps_snad.get(key)
        if d_ref is None or s_ref is None:
            unmatched += 1
            continue
        for a,b in zip(desc, d_ref):
            max_abs_desc = max(max_abs_desc, abs(a-b))
        for k in range(NCOEFF):
            for c in range(3):
                a = deriv[k*3+c]
                b = s_ref[c*NCOEFF+k]  # LAMMPS [xyz][k] layout
                max_abs_deriv = max(max_abs_deriv, abs(a-b))

    print(f"atoms compared        = {len(combined) - unmatched} / {len(combined)}")
    print(f"unmatched positions   = {unmatched}")
    print(f"max abs descriptor err = {max_abs_desc:.3e}")
    print(f"max abs derivative err = {max_abs_deriv:.3e}")

    if unmatched > 0 or max_abs_desc > 1e-6 or max_abs_deriv > 1e-4:
        print("FAIL: exastamp_combined.txt does not match LAMMPS")
        sys.exit(1)
    print("PASS: write_descriptor_snap's exastamp_combined.txt matches LAMMPS (sna/atom + snad/atom), any rank count")

if __name__ == "__main__":
    main()
