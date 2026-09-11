#!/usr/bin/env python3
"""Direct comparison of exastamp_snad.txt (compute_descriptor_snap's derivative aggregate,
written via dump_descriptor_snap_aggregate) against lammps_snad.txt (LAMMPS compute snad/atom) --
no reshaping or aggregation needed, both are already in the same [xyz-block][coeff]
per-atom row layout. Rows are matched by position (exaStamp ids are 0-indexed, LAMMPS
ids are 1-indexed).

Usage: python3 compare_aggregate.py
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
            vec = [float(v) for v in vals[4:]]
            rows.append((pos, vec))
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
                    vec = [float(v) for v in vals[4:]]
                    rows.append((pos, vec))
                i += 1
        else:
            i += 1
    return rows

def main():
    a = read_exastamp("exastamp_snad.txt")
    b = read_lammps("lammps_snad.txt")

    if len(a) != len(b):
        print(f"FAIL: atom count mismatch: exaStamp={len(a)} LAMMPS={len(b)}")
        sys.exit(1)

    lut = { tuple(round(c, 6) for c in pos): vec for pos, vec in b }

    max_abs = 0.0
    max_rel = 0.0
    unmatched = 0
    for pos, vec_a in a:
        key = tuple(round(c, 6) for c in pos)
        vec_b = lut.get(key)
        if vec_b is None:
            unmatched += 1
            continue
        if len(vec_a) != len(vec_b):
            print(f"FAIL: row length mismatch: exaStamp={len(vec_a)} LAMMPS={len(vec_b)}")
            sys.exit(1)
        for da, db in zip(vec_a, vec_b):
            diff = abs(da - db)
            max_abs = max(max_abs, diff)
            denom = max(abs(da), abs(db), 1e-8)
            max_rel = max(max_rel, diff / denom)

    print(f"atoms compared        = {len(a) - unmatched} / {len(a)}")
    print(f"unmatched positions   = {unmatched}")
    print(f"max abs component err = {max_abs:.3e}")
    print(f"max rel component err = {max_rel:.3e}")

    if unmatched > 0 or max_abs > 1e-4:
        print("FAIL: exaStamp's derivative aggregate does not match LAMMPS compute snad/atom")
        sys.exit(1)
    print("PASS: exaStamp's derivative aggregate (dump_descriptor_snap_aggregate) matches LAMMPS compute snad/atom directly")

if __name__ == "__main__":
    main()
