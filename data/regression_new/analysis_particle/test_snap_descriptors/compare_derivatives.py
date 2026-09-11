#!/usr/bin/env python3
"""Compare compute_descriptor_snap's compute_derivative CSR output (per-neighbor-pair
bispectrum Jacobian, dump_descriptor_snap_derivative) against LAMMPS's compute snad/atom,
for the same tantalum_noisy configuration already validated by compare_descriptors.py.

compute_descriptor_snap now also computes the derivative aggregate directly (the same
aggregation this script does in Python, done in C++ via atomic scatter) -- see
dump_descriptor_snap_aggregate + compare_aggregate.py for the simpler, direct comparison.
This script instead validates the finer-grained CSR (full per-neighbor-pair Jacobian,
strictly more information than snad/atom's aggregate) by reducing it to snad/atom's shape.

compute_descriptor_snap's raw per-row values use exaStamp's internal relative-position
convention (buf.drx = r_center - r_neighbor, matching LAMMPS compute_sna_atom.cpp's own
convention, confirmed by compare_descriptors.py's exact bispectrum match). LAMMPS's
compute snad/atom deliberately uses the opposite convention (rij = r_neighbor - r_center)
so it can accumulate its per-atom row directly with no extra sign flip. Working through
the chain rule for both conventions (see the project notes) reduces to one uniform rule:

  true dB_i/dR_k = -stored[block i, row k]   (for k == i itself, or k a neighbor of i)

So LAMMPS's aggregated snad[m] = sum over every row (in any block) whose neighbor-id is m,
of -stored[...], i.e.: group ALL rows (self rows included) by their neighbor-id column,
sum the raw stored vectors, then negate once at the end.

Usage: python3 compare_derivatives.py
"""
import sys

NCOEFF = 30  # twojmax=6, nelements=1 (see param.txt / ncoeff= printed by the .msp)

def read_positions(path, is_lammps):
    """id -> (x,y,z), from either exastamp_descriptors.txt or a LAMMPS dump custom file."""
    pos = {}
    if is_lammps:
        with open(path) as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith("ITEM: ATOMS"):
                i += 1
                while i < len(lines) and not lines[i].startswith("ITEM:"):
                    vals = lines[i].split()
                    if vals:
                        pos[int(vals[0])] = tuple(float(v) for v in vals[1:4])
                    i += 1
            else:
                i += 1
    else:
        with open(path) as f:
            for line in f:
                vals = line.split()
                if vals:
                    pos[int(vals[0])] = tuple(float(v) for v in vals[1:4])
    return pos

def build_id_map(pos_a, pos_b):
    lut = { tuple(round(c, 6) for c in p): idb for idb, p in pos_b.items() }
    idmap = {}
    for ida, p in pos_a.items():
        key = tuple(round(c, 6) for c in p)
        if key in lut:
            idmap[ida] = lut[key]
    return idmap

def read_exastamp_derivatives(path):
    """yields (center_id, [(neighbor_id, [90 floats]), ...]) per block."""
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        header = lines[i].split()
        assert header[0] == "ATOM"
        center_id = int(header[1])
        nrows = int(header[2])
        rows = []
        for r in range(nrows):
            i += 1
            vals = lines[i].split()
            rows.append((int(vals[0]), [float(v) for v in vals[1:]]))
        i += 1
        yield center_id, rows

def read_lammps_snad(path):
    """id -> [90 floats] in LAMMPS's [xyz-block][coeff] layout."""
    out = {}
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: ATOMS"):
            i += 1
            while i < len(lines) and not lines[i].startswith("ITEM:"):
                vals = lines[i].split()
                if vals:
                    out[int(vals[0])] = [float(v) for v in vals[4:]]
                i += 1
        else:
            i += 1
    return out

def main():
    pos_a = read_positions("exastamp_descriptors.txt", is_lammps=False)
    pos_b = read_positions("lammps_descriptors.txt", is_lammps=True)
    idmap = build_id_map(pos_a, pos_b)
    if len(idmap) != len(pos_a):
        print(f"FAIL: only matched {len(idmap)}/{len(pos_a)} atoms by position")
        sys.exit(1)

    # group every row (self rows included) by its neighbor-id column, in LAMMPS id space
    agg = { idb: [0.0] * (3 * NCOEFF) for idb in idmap.values() }
    for center_a, rows in read_exastamp_derivatives("exastamp_derivatives.txt"):
        for nbr_a, vec in rows:
            nbr_b = idmap[nbr_a]
            acc = agg[nbr_b]
            for i, v in enumerate(vec):
                acc[i] += v

    # negate once, and reorder our [k*3+xyz] layout into LAMMPS's [xyz*ncoeff+k] layout
    agg_true = {}
    for idb, vec in agg.items():
        out = [0.0] * (3 * NCOEFF)
        for k in range(NCOEFF):
            for c in range(3):
                out[c * NCOEFF + k] = -vec[k * 3 + c]
        agg_true[idb] = out

    lammps_snad = read_lammps_snad("lammps_snad.txt")

    max_abs = 0.0
    max_rel = 0.0
    for idb, ours in agg_true.items():
        theirs = lammps_snad[idb]
        for da, db in zip(ours, theirs):
            diff = abs(da - db)
            max_abs = max(max_abs, diff)
            denom = max(abs(da), abs(db), 1e-8)
            max_rel = max(max_rel, diff / denom)

    print(f"atoms compared        = {len(agg_true)} / {len(lammps_snad)}")
    print(f"max abs component err = {max_abs:.3e}")
    print(f"max rel component err = {max_rel:.3e}")

    if max_abs > 1e-4:
        print("FAIL: aggregated derivative does not match LAMMPS compute snad/atom")
        sys.exit(1)
    print("PASS: compute_descriptor_snap's derivative matches LAMMPS compute snad/atom")

if __name__ == "__main__":
    main()
