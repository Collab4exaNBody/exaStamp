# makes the slab validation systems from UO2_disturbed_charge.lmp : same atoms, z shifted by DZ in a box of height LZ
# (vacuum above and below, z non periodic), orthogonal and with an xy tilt (same fractional x,y coordinates)
import sys
src, out = sys.argv[1], sys.argv[2]
L, LZ, DZ, XY = 54.5, 66.0, 5.75, 5.0
lines = open(f"{src}/UO2_disturbed_charge.lmp").read().splitlines()
i0 = next(i for i,l in enumerate(lines) if l.startswith("Atoms")) + 2
atoms = [l.split() for l in lines[i0:] if l.strip()]
masses = lines[next(i for i,l in enumerate(lines) if l.startswith("Masses")):i0-2]
names = {"1":"U", "2":"O"}
for tag, xy in (("slab",0.0), ("slab_tri",XY)):
    with open(f"{out}/UO2_{tag}.lmp","w") as f:
        f.write(f"# disturbed UO2 (UO2_disturbed_charge.lmp) as a slab : z shifted by {DZ} ang in a box of height {LZ} ang"
                + (f", xy tilt {xy} ang (same fractional x,y)" if xy else "") + "\n\n")
        f.write(f"{len(atoms)} atoms\n2 atom types\n\n0.0 {L} xlo xhi\n0.0 {L} ylo yhi\n0.0 {LZ} zlo zhi\n")
        if xy: f.write(f"{xy} 0.0 0.0 xy xz yz\n")
        f.write("\n" + "\n".join(masses) + "\n\nAtoms  # charge\n\n")
        xyz = []
        for a in atoms:
            x, y, z = float(a[3]), float(a[4]), float(a[5])
            sx, sy = x/L, y/L
            X, Y, Z = sx*L + sy*xy, sy*L, z + DZ
            f.write(f"{a[0]} {a[1]} {a[2]} {X:.10f} {Y:.10f} {Z:.10f}\n")
            xyz.append(f"{names[a[1]]} {X:.10f} {Y:.10f} {Z:.10f}")
    with open(f"{out}/UO2_{tag}_ext.xyz","w") as f:
        f.write(f"{len(atoms)}\nLattice=\"{L} 0.0 0.0 {xy} {L} 0.0 0.0 0.0 {LZ}\" Properties=species:S:1:pos:R:3 pbc=\"T T F\"\n")
        f.write("\n".join(xyz) + "\n")
print("ok")
