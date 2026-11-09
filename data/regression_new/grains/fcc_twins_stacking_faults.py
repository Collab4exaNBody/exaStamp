# Writes fcc_twins_sf.xyz for compute_grain_segmentation_twins.msp: FCC Cu stacked along [111]
# (x = [1-10], y = [11-2], z = [111], periodic), two twins (two coherent twin planes, one is the
# periodic boundary) and one intrinsic stacking fault inside the first twin, small random noise.
# Without coherent-interface handling: 6 grains (3 FCC slabs + 3 HCP planes); with it: 2 grains.
import numpy as np

a = 3.615                      # Cu lattice parameter, ang
d = a / np.sqrt(2)             # nearest-neighbor distance
h = a / np.sqrt(3)             # (111) layer spacing
nx, ny = 12, 7
stacking = "ABC" * 8 + "BAC" * 8              # twin planes at layers 23 and 47
stacking = stacking[:11] + stacking[12:]      # one C layer removed: intrinsic stacking fault
shift = { 'A': 0.0, 'B': d / np.sqrt(3), 'C': 2 * d / np.sqrt(3) }

box = np.array( [ nx * d, ny * d * np.sqrt(3), len(stacking) * h ] )
pos = np.array( [ ( (i + u) * d, (j + u) * d * np.sqrt(3) + shift[c], l * h )
                  for l, c in enumerate(stacking) for i in range(nx) for j in range(ny) for u in (0.0, 0.5) ] )
pos += np.random.default_rng(7).normal( 0.0, 0.03, pos.shape )
pos %= box

with open( "fcc_twins_sf.xyz", "w" ) as f:
    f.write( f"{len(pos)}\nLattice=\"{box[0]} 0 0 0 {box[1]} 0 0 0 {box[2]}\" Properties=species:S:1:pos:R:3\n" )
    for p in pos: f.write( f"Cu {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n" )
