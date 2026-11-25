# Writes fcc_planar_defects.xyz for the sfd_gmm_*.msp inputs: FCC Cu stacked along [111]
# (x = [1-10], y = [11-2], z = [111]), periodic in all directions by construction: two coherent
# twin planes (one at the periodic boundary), one intrinsic stacking fault in the first twin,
# 24 vacancies away from the planar defects, and Gaussian noise (0.05 ang, thermal-like).
#   python3 fcc_planar_defects.py      # writes fcc_planar_defects.xyz
import numpy as np

a = 3.615                      # Cu lattice parameter, ang
d = a / np.sqrt(2)             # nearest-neighbor distance
h = a / np.sqrt(3)             # (111) layer spacing
nx, ny = 16, 10
stacking = "ABC" * 12 + "BAC" * 12
stacking = stacking[:17] + stacking[18:]      # one layer removed: intrinsic stacking fault (layers 16-17)
# twin planes (hcp-like layers): 34 and 70 (next to the periodic boundary)
shift = { 'A': 0.0, 'B': d / np.sqrt(3), 'C': 2 * d / np.sqrt(3) }

box = np.array( [ nx * d, ny * d * np.sqrt(3), len(stacking) * h ] )
pos = np.array( [ ( (i + u) * d, (j + u) * d * np.sqrt(3) + shift[c], l * h )
                  for l, c in enumerate(stacking) for i in range(nx) for j in range(ny) for u in (0.0, 0.5) ] )
layer = np.repeat( np.arange( len(stacking) ), nx * ny * 2 )

rng = np.random.default_rng( 11 )
# vacancies in layers at least 4 planes away from the hcp-like layers 16, 17, 34 and 70
nl = len(stacking)
far = [ l for l in range( nl ) if min( min( abs(l - k), nl - abs(l - k) ) for k in ( 16, 17, 34, 70 ) ) >= 4 ]
cand = np.flatnonzero( np.isin( layer, far ) )
keep = np.ones( len(pos), bool )
keep[ rng.choice( cand, 24, replace=False ) ] = False
pos = pos[keep]
pos += rng.normal( 0.0, 0.05, pos.shape )
pos %= box

with open( "fcc_planar_defects.xyz", "w" ) as f:
    f.write( f"{len(pos)}\nLattice=\"{box[0]} 0 0 0 {box[1]} 0 0 0 {box[2]}\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n" )
    for p in pos: f.write( f"Cu {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n" )
