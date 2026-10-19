"""Generate an FAtomes file of trans-1,4-polybutadiene chains (CH3-[CH2-CH=CH-CH2]n-CH3, types CT/CM/HT/HM),
using an existing PB FAtomes file as template (species, pair potentials, force field, intramolecular weights).

Chains are grown carbon by carbon from the force field equilibrium bond lengths and angles, with random backbone
torsions. Hydrogens are placed as soon as both neighbors of a carbon are known, and every new atom is checked against
all atoms already placed (periodic box), except atoms within 3 backbone bonds. Dead ends are handled by stepping
back a few carbons. The configuration is not equilibrated : it is meant for molecule numbering / unwrapping tests
and short runs.

usage: python make_pb_chains.py <template FAtomes> <output FAtomes> [n_chains] [n_units] [box_ang] [seed]
"""
import sys
import numpy as np

tmpl, out = sys.argv[1], sys.argv[2]
n_chains = int(sys.argv[3]) if len(sys.argv) > 3 else 20
n_units = int(sys.argv[4]) if len(sys.argv) > 4 else 60
box = float(sys.argv[5]) if len(sys.argv) > 5 else 60.0
rng = np.random.default_rng(int(sys.argv[6]) if len(sys.argv) > 6 else 1234)

# minimum non-bonded distances (ang) between heavy (C) and light (H) atoms
DMIN = {("C", "C"): 3.2, ("C", "H"): 2.6, ("H", "C"): 2.6, ("H", "H"): 2.2}
CUT = 3.2

# force field equilibrium values (bond lengths in ang, angles in degrees)
B = {("CT", "CT"): 1.529, ("CT", "CM"): 1.510, ("CM", "CM"): 1.340, ("CT", "HT"): 1.090, ("CM", "HM"): 1.080}
A = {("CT", "CM", "CM"): 124.0, ("CT", "CT", "CM"): 111.1}
CHARGE = {"CT": -0.120, "CT3": -0.180, "HT": 0.060, "CM": -0.115, "HM": 0.115}
TETRA = np.radians(109.47)


def bond(a, b):
    return B.get((a, b), B.get((b, a)))


def angle(a, b, c):
    return A.get((a, b, c), A.get((c, b, a)))


def unit(v):
    return v / np.linalg.norm(v)


def place(a, b, c, length, angle_deg, torsion_deg):
    """NeRF : position of d with |cd| = length, angle(b,c,d) = angle, dihedral(a,b,c,d) = torsion"""
    th, ph = np.radians(angle_deg), np.radians(torsion_deg)
    bc = unit(c - b)
    n = unit(np.cross(b - a, bc))
    m = np.cross(n, bc)
    return c - length * np.cos(th) * bc + length * np.sin(th) * (np.cos(ph) * m + np.sin(ph) * n)


def hydrogens(t, b, neighbors):
    """hydrogen positions of carbon b (type t) bonded to the given carbon neighbor positions"""
    if t == "CM":                                            # sp2 : opposite to the bisector of the two neighbors
        u = unit(unit(b - neighbors[0]) + unit(b - neighbors[1]))
        return [("HM", b + bond("CM", "HM") * u)]
    if len(neighbors) == 1:                                  # CH3 chain end : 3 H around the C-C axis
        ax = unit(b - neighbors[0])
        p = np.cross(ax, [1.0, 0.0, 0.0])
        if np.linalg.norm(p) < 0.1: p = np.cross(ax, [0.0, 1.0, 0.0])
        p = unit(p); q = np.cross(ax, p)
        return [("HT", b + bond("CT", "HT") * (np.cos(np.pi - TETRA) * ax + np.sin(np.pi - TETRA) * (np.cos(f) * p + np.sin(f) * q)))
                for f in (0.0, 2 * np.pi / 3, 4 * np.pi / 3)]
    a1, a2 = unit(b - neighbors[0]), unit(b - neighbors[1])  # CH2 : 2 H in the plane bisecting the C-C-C angle
    u = unit(a1 + a2); n = unit(np.cross(a1, a2)); half = np.radians(107.8 / 2)
    return [("HT", b + bond("CT", "HT") * (np.cos(half) * u + s * np.sin(half) * n)) for s in (1, -1)]


class Atoms:
    """atoms already placed (all chains), with a cell list for distance checks"""
    def __init__(self):
        self.nc = max(1, int(box // CUT))
        self.cells = {}

    def key(self, r):
        return tuple((np.floor((r % box) / box * self.nc).astype(int)) % self.nc)

    def add(self, r, cls):
        self.cells.setdefault(self.key(r), []).append((r % box, cls))

    def clash(self, r, cls):
        k = self.key(r)
        w = r % box
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for (q, c2) in self.cells.get(((k[0] + dx) % self.nc, (k[1] + dy) % self.nc, (k[2] + dz) % self.nc), ()):
                        d = w - q; d -= box * np.round(d / box)
                        if d @ d < DMIN[(cls, c2)] ** 2: return True
        return False


def chain_clash(r, cls, bb, hs, i):
    """new atom (belonging to backbone index i) against atoms of its own chain more than 3 backbone bonds away"""
    def d2(q):                                               # minimum image : the chain may wrap around the box
        d = r - q; d -= box * np.round(d / box); return d @ d
    for j in range(0, max(0, i - 3)):
        if d2(bb[j]) < DMIN[(cls, "C")] ** 2: return True
        for (_, h) in hs[j]:
            if d2(h) < DMIN[(cls, "H")] ** 2: return True
    return False


def ok(r, cls, bb, hs, i, world):
    return not world.clash(r, cls) and not chain_clash(r, cls, bb, hs, i)


def torsion_choices(types, i):
    b, c = types[i - 2], types[i - 1]
    if b == "CM" and c == "CM": return [180.0]               # trans double bond
    if b == "CT" and c == "CT": return [180.0, 60.0, -60.0]  # CM-CT-CT-CM : trans / gauche
    return [120.0, -120.0, 180.0]                            # around CT-CM single bonds : skew


world = Atoms()
chains = []
for ic in range(n_chains):
    types = ["CT", "CM", "CM", "CT"] * n_units
    nb = len(types)
    for attempt in range(1000):
        r0 = rng.uniform(0, box, 3)
        u = unit(rng.normal(size=3)); v = unit(np.cross(u, rng.normal(size=3)))
        r1 = r0 + bond(types[0], types[1]) * u
        th = np.radians(180.0 - angle(types[0], types[1], types[2]))
        r2 = r1 + bond(types[1], types[2]) * (np.cos(th) * u + np.sin(th) * v)
        bb, hs = [r0, r1, r2], [[], [], []]
        if any(world.clash(r, "C") for r in bb): continue
        hs[0] = hydrogens(types[0], r0, [r1]); hs[1] = hydrogens(types[1], r1, [r0, r2])
        if any(world.clash(h, "H") for (_, h) in hs[0] + hs[1]): continue
        failures = 0
        while len(bb) < nb and failures <= 400:
            i = len(bb)
            placed = False
            for t in torsion_choices(types, i) * 3:
                if not (types[i - 2] == "CM" and types[i - 1] == "CM"): t += rng.normal(0, 8.0)
                r = place(bb[i - 3], bb[i - 2], bb[i - 1], bond(types[i - 1], types[i]), angle(types[i - 2], types[i - 1], types[i]), t)
                if not ok(r, "C", bb, hs, i, world): continue
                h_prev = hydrogens(types[i - 1], bb[i - 1], [bb[i - 2], r])      # H of the previous carbon
                h_last = hydrogens(types[i], r, [bb[i - 1]]) if i == nb - 1 else []
                if all(ok(h, "H", bb, hs, i - 1, world) for (_, h) in h_prev + h_last):
                    hs[i - 1] = h_prev; bb.append(r); hs.append(h_last); placed = True; break
            if not placed:                                     # dead end : step back a few carbons
                failures += 1
                back = max(3, len(bb) - int(rng.integers(2, 9)))
                del bb[back:]; del hs[back:]
                hs[back - 1] = []                              # H of the new last carbon are recomputed
        if len(bb) == nb: break
    else:
        sys.exit("could not place chain %d, lower the density" % ic)
    for k in range(nb):
        world.add(bb[k], "C")
        for (_, h) in hs[k]: world.add(h, "H")
    chains.append((types, bb, hs))
    print("chain %d placed (attempt %d)" % (ic, attempt), file=sys.stderr)

# atom list and connectivity (ids follow chain order : carbons and their hydrogens)
atoms, con = [], []


def add_atom(t, r, q):
    atoms.append((t, r, q)); con.append([]); return len(atoms) - 1


def link(i, j):
    con[i].append(j); con[j].append(i)


for types, bb, hs in chains:
    nb = len(types)
    prev = None
    for k in range(nb):
        end = (k == 0 or k == nb - 1)
        c = add_atom(types[k], bb[k], CHARGE["CT3"] if (types[k] == "CT" and end) else CHARGE[types[k]])
        if prev is not None: link(prev, c)
        for (ht, h) in hs[k]: link(c, add_atom(ht, h, CHARGE[ht]))
        prev = c

n = len(atoms)
assert max(len(c) for c in con) <= 4
print("atoms:", n, " chains:", len(chains), " total charge: %.2e" % sum(a[2] for a in atoms), file=sys.stderr)

# write from template
L = open(tmpl, errors="replace").read().split("\n")
o = []
i = 0
while i < len(L):
    t = L[i].split()
    key = t[0] if t else ""
    if key == "maille_long":
        o.append("maille_long\t%.6f %.6f %.6f ang" % (box, box, box)); i += 1; continue
    if key == "maille_angle":
        o.append("maille_angle\t90.0 90.0 90.0 degre"); i += 1; continue
    if key == "PositionDesAtomesCart":
        o.append(L[i]); o.append(str(n))
        for (ty, r, q) in atoms:
            rw = r % box
            o.append("  %-4s %16.10f %16.10f %16.10f" % (ty, rw[0], rw[1], rw[2]))
        i += 2 + int(L[i + 1]); continue
    if key == "VitesseDesAtomes":
        o.append(L[i]); o.extend(["0.0 0.0 0.0"] * n); i += 1
        while i < len(L) and len(L[i].split()) == 3 and not L[i].startswith("*"): i += 1
        continue
    if key == "ModificationChargeDesAtomes":
        o.append(L[i]); o.append(str(n))
        o.extend("%d %.6f" % (k, a[2]) for k, a in enumerate(atoms))
        i += 2 + int(L[i + 1]); continue
    if key == "Zmatrice":
        o.append(L[i]); o.append(str(n))
        o.extend(" ".join(str(x) for x in [k] + c) for k, c in enumerate(con))
        i += 2 + int(L[i + 1]); continue
    if key == "Reactivite":                                  # polymerisation data of the template, not needed
        while o and o[-1].startswith("*"): o.pop()
        break
    o.append(L[i]); i += 1
open(out, "w").write("\n".join(o) + "\n")
