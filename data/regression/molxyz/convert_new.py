#!/usr/bin/python3

import sys

lines = open(sys.argv[1]).readlines()
natoms = int(lines[0])
print(natoms)
print(lines[1],end="")
id=1
for l in lines[2:]:
    d = l.split(' ')
    typeMol = d[0].split('_')[0]
    typeAtom = d[1]
    molId = d[2]
    X = float(d[4])
    Y = float(d[5])
    Z = float(d[6])
    c0 = int(d[7])
    c1 = int(d[8])
    c2 = int(d[9])
    c3 = int(d[10])
    print(typeMol,typeAtom,molId,id,X,Y,Z,c0,c1,c2,c3)
    id = id+1

