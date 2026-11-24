"""Water configuration for the DeePMD tests (generated, no third-party data): n^3 molecules on a jittered cubic lattice, random orientations.
    python water_frame.py n seed out   -> out.npz (coord, box, atype; O=0 H=1) and out.xyz (extended xyz, species O/H)"""
import sys, numpy as np
n, seed, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rng = np.random.default_rng(seed)
nmol = n ** 3
L = (nmol * 18.015 / 0.6022 / 1.0) ** (1 / 3)   # density 1 g/cm3
a = L / n
th = np.deg2rad(104.5) / 2
local = np.array([[0, 0, 0], [0.9572 * np.sin(th), 0.9572 * np.cos(th), 0], [-0.9572 * np.sin(th), 0.9572 * np.cos(th), 0]])
coord, atype = [], []
for i in range(n):
    for j in range(n):
        for k in range(n):
            q = rng.normal(size=4); q /= np.linalg.norm(q)
            w, x, y, z = q
            R = np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)], [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)], [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])
            c = (np.array([i, j, k]) + 0.5) * a + rng.uniform(-0.15, 0.15, 3) * a
            for p, t in zip(local, (0, 1, 1)):
                coord.append(c + R @ p); atype.append(t)
coord = np.mod(np.array(coord), L); atype = np.array(atype)
box = np.diag([L, L, L])
np.savez(out + ".npz", coord=coord, box=box, atype=atype)
with open(out + ".xyz", "w") as fh:
    fh.write(f"{len(coord)}\nLattice=\"{L} 0 0 0 {L} 0 0 0 {L}\" Properties=species:S:1:pos:R:3\n")
    for c, t in zip(coord, atype):
        fh.write(f"{'OH'[t]} {c[0]:.10f} {c[1]:.10f} {c[2]:.10f}\n")
print(out, len(coord), "atoms, L =", L)
