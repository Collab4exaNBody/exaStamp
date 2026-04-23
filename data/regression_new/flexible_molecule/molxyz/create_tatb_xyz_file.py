import numpy as np
from math import sqrt,sin,cos,tan,acos,atan,pi,floor,ceil,log

datamol=np.loadtxt('tatb_isolated_molecule.dat',dtype='str')
conn=np.loadtxt('connectivity_table.dat')

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

ofile=open('triclinic_unit_cell.xyz','w')
ofile.write('48\n')
#ofile.write('Lattice=\"%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\"\n' %(Ht[0,0], Ht[0,1], Ht[0,2], Ht[1,0], Ht[1,1], Ht[1,2], Ht[2,0], Ht[2,1], Ht[2,2]))
ofile.write('%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n' %(Ht[0,0], Ht[0,1], Ht[0,2], Ht[1,0], Ht[1,1], Ht[1,2], Ht[2,0], Ht[2,1], Ht[2,2]))

for i in range(Nmols_ucell):
    for j in range(Nat_per_mol):
        pos=pos_ucell[i,j,:]
        invpos=np.dot(np.linalg.inv(H),pos)
        
        if(invpos[0]<0.):
            invpos[0]+=1.
        elif(invpos[0]>1.):
            invpos[0]-=1.
            
        if(invpos[1]<0.):
            invpos[1]+=1.
        elif(invpos[1]>1.):
            invpos[1]-=1.
            
        if(invpos[2]<0.):
            invpos[2]+=1.
        elif(invpos[2]>1.):
            invpos[2]-=1.
            
        pos=np.dot(H,invpos)
        px=pos[0]
        py=pos[1]
        pz=pos[2]
        specie=datamol[i*Nat_per_mol+j,0]
        if (conn[j,0] != -1):
            c1=conn[j,0]+i*Nat_per_mol
        else:
            c1=conn[j,0]
        if (conn[j,1] != -1):
            c2=conn[j,1]+i*Nat_per_mol
        else:
            c2=conn[j,1]
        if (conn[j,2] != -1):
            c3=conn[j,2]+i*Nat_per_mol
        else:
            c3=conn[j,2]
        if (conn[j,3] != -1):
            c4=conn[j,3]+i*Nat_per_mol
        else:
            c4=conn[j,3]            
        ofile.write('%s %s %d %d %5.8f %5.8f %5.8f %d %d %d %d\n' %(specie, specie.replace('TATB_', ''), i+1, i, px, py, pz, c1, c2, c3, c4))
ofile.close()
