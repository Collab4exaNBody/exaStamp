#!/usr/bin/env python3
# Writes igar_field.vtk : smooth periodic scalar field 0.5*(1+cos(2pi x/L) cos(2pi y/L) cos(2pi z/L)) sampled at the
# centers of n^3 voxels (ASCII legacy VTK structured points, doubles, F order), read by igar_read_field with
# storage: double. usage: make_igar_field.py [n=16]
import math, sys
n = int(sys.argv[1]) if len(sys.argv) > 1 else 16
with open("igar_field.vtk", "w") as f:
    f.write("# vtk DataFile Version 3.0\nIGAR regression field\nASCII\nDATASET STRUCTURED_POINTS\n")
    f.write(f"DIMENSIONS {n} {n} {n}\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA {n**3}\nSCALARS ScalarData double\nLOOKUP_TABLE default\n")
    for k in range(n):
        for j in range(n):
            for i in range(n):
                c = lambda m: math.cos(2.0 * math.pi * (m + 0.5) / n)
                f.write("%.12f\n" % (0.5 * (1.0 + c(i) * c(j) * c(k))))
