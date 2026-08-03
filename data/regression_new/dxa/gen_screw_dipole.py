import numpy as np

a = 3.3024  # BCC Ta lattice parameter, matches quadrupole_dislo.xyz (132.096/40)

# orthogonal frame for a <111> screw dislocation: line along z=[111]
e1 = np.array([1., -1., 0.]); e1 /= np.linalg.norm(e1)   # x
e2 = np.array([1., 1., -2.]); e2 /= np.linalg.norm(e2)   # y
e3 = np.array([1., 1., 1.]);  e3 /= np.linalg.norm(e3)   # z (line direction)
R = np.array([e1, e2, e3])   # rotated = R @ cubic

ax = a * np.sqrt(2.0)
ay = a * np.sqrt(6.0)
az = a * np.sqrt(3.0)   # = 2*|b|, minimal periodic repeat along the line

nx, ny, nz = 30, 18, 20
Lx, Ly, Lz = nx*ax, ny*ay, nz*az
print("box:", Lx, Ly, Lz)

b_mag = a * np.sqrt(3.0) / 2.0   # |a/2<111>|
print("burgers magnitude:", b_mag)

# generate BCC lattice (2-atom basis, cubic units of a) over a generous cubic-index range,
# rotate into the orthogonal frame, keep only points landing inside the target box
pad = 60
idx = np.arange(-pad, pad+1)
I, J, K = np.meshgrid(idx, idx, idx, indexing='ij')
I = I.ravel(); J = J.ravel(); K = K.ravel()

basis = np.array([[0.0,0.0,0.0],[0.5,0.5,0.5]])
pts = []
for bx,by,bz in basis:
    cubic = np.stack([I+bx, J+by, K+bz], axis=1) * a
    rot = cubic @ R.T
    pts.append(rot)
rot_all = np.concatenate(pts, axis=0)

tol = 1e-6 * a
mask = (rot_all[:,0] >= -tol) & (rot_all[:,0] < Lx-tol) & \
       (rot_all[:,1] >= -tol) & (rot_all[:,1] < Ly-tol) & \
       (rot_all[:,2] >= -tol) & (rot_all[:,2] < Lz-tol)
atoms = rot_all[mask]
print("n atoms (perfect lattice):", len(atoms))

# screw dislocation dipole, isotropic elasticity exact solution: u_z = b/(2*pi) * theta
# two lines along z, opposite sign, same y0, separated along x -> periodic-compatible dipole
x0off = 0.13 * a
y0 = Ly/2.0 + x0off
x1 = Lx*0.25 + x0off
x2 = Lx*0.75 + x0off

x = atoms[:,0]; y = atoms[:,1]; z = atoms[:,2]
theta1 = np.arctan2(y - y0, x - x1)
theta2 = np.arctan2(y - y0, x - x2)
uz = b_mag/(2.0*np.pi) * (theta1 - theta2)
z_new = (z + uz) % Lz

out = np.stack([x, y, z_new], axis=1)

with open('/home/lafourcadep/CODES/ATOMISTIC/EXANBODY/exaStamp/data/regression_new/delaunay/screw_dislo_dipole.xyz', 'w') as f:
    f.write(f"{len(out)}\n")
    f.write(f'Lattice="{Lx} 0.0 0.0 0.0 {Ly} 0.0 0.0 0.0 {Lz}" Properties=species:S:1:pos:R:3\n')
    for p in out:
        f.write(f"Ta {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}\n")

print("wrote", len(out), "atoms")
print("dislocation cores at (x,y):", (x1,y0), "and", (x2,y0), "separation:", x2-x1)
