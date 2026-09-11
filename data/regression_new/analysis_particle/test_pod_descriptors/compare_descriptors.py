#!/usr/bin/env python3
"""Compare exastamp_pod_descriptors.txt (compute_descriptor_pod) against
lammps_pod_descriptors.txt (LAMMPS compute pod/atom) for the same Ta_disturbed
configuration. Matches rows by nearest position (robust to any id/ordering
convention difference between the two codes) rather than assuming a fixed
id offset.

Usage: python3 compare_descriptors.py [exastamp_file] [lammps_file]
"""
import sys

def read_exastamp(path):
    rows = []
    with open(path) as f:
        for line in f:
            vals = line.split()
            if not vals:
                continue
            pos = tuple(float(v) for v in vals[1:4])
            desc = [float(v) for v in vals[4:]]
            rows.append((pos, desc))
    return rows

def read_lammps(path):
    rows = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: ATOMS"):
            i += 1
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                vals = lines[i].split()
                if vals:
                    pos = tuple(float(v) for v in vals[1:4])
                    desc = [float(v) for v in vals[4:]]
                    rows.append((pos, desc))
                i += 1
        else:
            i += 1
    return rows

def main():
    exastamp_file = sys.argv[1] if len(sys.argv) > 1 else "exastamp_pod_descriptors.txt"
    lammps_file   = sys.argv[2] if len(sys.argv) > 2 else "lammps_pod_descriptors.txt"

    a = read_exastamp(exastamp_file)
    b = read_lammps(lammps_file)

    if len(a) != len(b):
        print(f"FAIL: atom count mismatch: exaStamp={len(a)} LAMMPS={len(b)}")
        sys.exit(1)
    ncoeff = len(a[0][1])
    if ncoeff != len(b[0][1]):
        print(f"FAIL: descriptor length mismatch: exaStamp={ncoeff} LAMMPS={len(b[0][1])}")
        sys.exit(1)

    # index LAMMPS rows by rounded position for fast nearest-position matching
    # (box is periodic and noise is small, so plain rounding is enough here --
    # positions come from the same source file, just re-serialized). 5 decimal
    # places (down from 6) avoids spurious mismatches when a coordinate lands
    # exactly on a rounding tie (e.g. 16.5635725) and the two codes' printf-style
    # float formatting break the tie in opposite directions -- atoms are always
    # far more than 1e-5 apart, so this loses no real disambiguating power.
    lut = {}
    for pos, desc in b:
        key = tuple(round(c, 5) for c in pos)
        lut[key] = desc

    max_abs = 0.0
    max_rel = 0.0
    unmatched = 0
    for pos, desc_a in a:
        key = tuple(round(c, 5) for c in pos)
        desc_b = lut.get(key)
        if desc_b is None:
            unmatched += 1
            continue
        for da, db in zip(desc_a, desc_b):
            diff = abs(da - db)
            max_abs = max(max_abs, diff)
            denom = max(abs(da), abs(db), 1e-12)
            max_rel = max(max_rel, diff / denom)

    print(f"atoms compared        = {len(a) - unmatched} / {len(a)}")
    print(f"unmatched positions   = {unmatched}")
    print(f"max abs component err = {max_abs:.3e}")
    print(f"max rel component err = {max_rel:.3e}")

    if unmatched > 0 or max_abs > 1e-6:
        print("FAIL: descriptors do not match")
        sys.exit(1)
    print("PASS: exaStamp compute_descriptor_pod matches LAMMPS compute pod/atom")

if __name__ == "__main__":
    main()
