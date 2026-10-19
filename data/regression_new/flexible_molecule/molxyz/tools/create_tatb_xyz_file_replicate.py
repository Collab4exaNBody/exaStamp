import numpy as np
from math import sqrt,sin,cos,tan,acos,atan,pi,floor,ceil,log

# TATB crystal supercell in read_xyz_molecules format.  usage (from this folder) : python create_tatb_xyz_file_replicate.py Na Nb Nc
import os
here=os.path.dirname(os.path.abspath(__file__))
datamol=np.loadtxt(os.path.join(here,'tatb_isolated_molecule.dat'),dtype='str')
conn=np.loadtxt(os.path.join(here,'connectivity_table.dat'))

a=8.900
b=8.937
c=6.947
alpha=111.3
beta=86.7
gamma=120.0

alpha *= (pi/180.)
beta  *= (pi/180.)
gamma *= (pi/180.)

H = np.array([[ a,     b*cos(gamma),     c*cos(beta)],
              [ 0.,    b*sin(gamma),     c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)],
              [ 0.,              0.,     c*sqrt(1.-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)+2.*cos(alpha)*cos(beta)*cos(gamma))/sin(gamma)]])

Ht=np.transpose(H)
vec_a=Ht[0,:]
vec_b=Ht[1,:]
vec_c=Ht[2,:]

px=datamol[:,1].astype(float)
py=datamol[:,2].astype(float)
pz=datamol[:,3].astype(float)

Nat_per_mol=24
Nmols_ucell=2
pos_ucell=np.zeros((Nmols_ucell,Nat_per_mol,3))
pos_ucell[0,:,0]=px[:Nat_per_mol]
pos_ucell[0,:,1]=py[:Nat_per_mol]
pos_ucell[0,:,2]=pz[:Nat_per_mol]
pos_ucell[1,:,0]=px[Nat_per_mol:]
pos_ucell[1,:,1]=py[Nat_per_mol:]
pos_ucell[1,:,2]=pz[Nat_per_mol:]

import sys

Na=int(sys.argv[1])
Nb=int(sys.argv[2])
Nc=int(sys.argv[3])

# atom charges of the TATB force field (../tatb_ff.cfg)
charge = { 'Cno2': 0.34068, 'Cnh2': -0.243, 'No2': 0.6798, 'Nh2': -0.693, 'O': -0.51024, 'H': 0.468 }

# read_xyz_molecules format : cell matrix (row major, columns are the cell vectors), then one line per atom :
# molecule type, atom type, molecule id, atom id, x y z, 5 bonded atom ids (-1 = none), charge
ofile=open('tatb_triclinic_%sx%sx%s.xyz' %(str(Na), str(Nb), str(Nc)),'w')
ofile.write('%d\n' %(Nat_per_mol*Nmols_ucell*Na*Nb*Nc))
S = H * np.array([Na,Nb,Nc])   # supercell : column k is Nk times cell vector k
ofile.write(' '.join('%.10f' %(S[r,s]) for r in range(3) for s in range(3)) + '\n')
molid=1
Ncell=0
for m in range(Na):
    for n in range(Nb):
        for p in range(Nc):
            shift = m*vec_a + n*vec_b + p*vec_c
            for i in range(Nmols_ucell):
                offset_conn=(Ncell*Nmols_ucell+i)*Nat_per_mol
                for j in range(Nat_per_mol):
                    pos=pos_ucell[i,j,:]+shift
                    atom=datamol[i*Nat_per_mol+j,0].replace('TATB_','')
                    # connectivity table : 1-based atom ids within the molecule -> global 1-based atom ids
                    c=[ int(conn[j,k])+offset_conn if conn[j,k] != -1 else -1 for k in range(4) ]
                    ofile.write('TATB %s %d %d %.8f %.8f %.8f %d %d %d %d -1 %.6f\n' %(atom, molid, offset_conn+j+1, pos[0], pos[1], pos[2], c[0], c[1], c[2], c[3], charge[atom]))
                molid+=1
            Ncell+=1
ofile.close()
