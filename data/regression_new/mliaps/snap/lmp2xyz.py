#!/usr/bin/env python3
"""Convert LAMMPS 'write_dump custom id type x y z vx vy vz' (orthogonal box, origin 0) to extended XYZ.
usage: lmp2xyz.py init.dump init.xyz Elem1 [Elem2 ...]   (ElemN = name of LAMMPS type N)"""
import sys
lines = open(sys.argv[1]).read().splitlines()
elems = sys.argv[3:]
n = int(lines[lines.index("ITEM: NUMBER OF ATOMS") + 1])
b = lines.index(next(l for l in lines if l.startswith("ITEM: BOX BOUNDS")))
lo_hi = [lines[b + 1 + i].split() for i in range(3)]
assert all(float(lo) == 0.0 for lo, _ in lo_hi), "box origin must be 0"
L = [hi for _, hi in lo_hi]
atoms = lines[lines.index(next(l for l in lines if l.startswith("ITEM: ATOMS"))) + 1:][:n]
with open(sys.argv[2], "w") as f:
    f.write(f"{n}\n")
    f.write(f'Lattice="{L[0]} 0 0 0 {L[1]} 0 0 0 {L[2]}" Properties=species:S:1:pos:R:3:velo:R:3\n')
    for a in atoms:
        t = a.split()
        f.write(" ".join([elems[int(t[1]) - 1]] + t[2:8]) + "\n")
