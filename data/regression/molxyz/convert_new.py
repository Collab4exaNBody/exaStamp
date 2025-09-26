#!/usr/bin/python3

import sys

lines = open(sys.argv[1]).readlines()
natoms = int(lines[0])
print(natoms)
print(lines[1],end="")
id=1
for l in lines[2:]:
    d = l.split(' ')
    typeAtom = d[0]
    X = float(d[1])
    Y = float(d[2])
    Z = float(d[3])
    at_id = int(d[4])
    typeMol = d[5]
    mol_id = int(d[6])
    c0 = int(d[8])
    c1 = int(d[9])
    c2 = int(d[10])
    c3 = int(d[11])
    c4 = int(d[12])
    charge = float(d[13])
    print(typeMol,typeAtom,mol_id,at_id,X,Y,Z,c0,c1,c2,c3,c4,charge)

