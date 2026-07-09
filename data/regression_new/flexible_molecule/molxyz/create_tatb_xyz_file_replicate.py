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

ofile=open('triclinic_unit_cell_%sx%sx%s_normal.xyz' %(str(Na), str(Nb), str(Nc)),'w')
ofile.write('%d\n' %(Nat_per_mol*Nmols_ucell*Na*Nb*Nc))
#ofile.write('Lattice=\"%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\"\n' %(Na*Ht[0,0], Na*Ht[0,1], Na*Ht[0,2], Nb*Ht[1,0], Nb*Ht[1,1], Nb*Ht[1,2], Nc*Ht[2,0], Nc*Ht[2,1], Nc*Ht[2,2]))
#ofile.write('100. 0. 0. 0. 100. 0. 0. 0. 100.\n')
ofile.write('%5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f %5.4f\n' %(Na*Ht[0,0], Na*Ht[0,1], Na*Ht[0,2], Nb*Ht[1,0], Nb*Ht[1,1], Nb*Ht[1,2], Nc*Ht[2,0], Nc*Ht[2,1], Nc*Ht[2,2]))
molid=1
Ncell=0
for m in range(Na):
    for n in range(Nb):
        for p in range(Nc):
            shift_a=m*vec_a
            shift_b=n*vec_b
            shift_c=p*vec_c
            offset_idmol=((m+1)*(n+1)*(p+1)-1)*Nmols_ucell
            Nmol=0
            for i in range(Nmols_ucell):
                moltype=i
                offset_conn=(Ncell*Nmols_ucell+i)*Nat_per_mol                
                for j in range(Nat_per_mol):
                    
                    pos=pos_ucell[i,j,:]
                    # invpos=np.dot(np.linalg.inv(H),pos)
        
                    # if(invpos[0]<0.):
                    #     invpos[0]+=1.
                    # elif(invpos[0]>1.):
                    #     invpos[0]-=1.
            
                    # if(invpos[1]<0.):
                    #     invpos[1]+=1.
                    # elif(invpos[1]>1.):
                    #     invpos[1]-=1.
            
                    # if(invpos[2]<0.):
                    #     invpos[2]+=1.
                    # elif(invpos[2]>1.):
                    #     invpos[2]-=1.
                    # pos=np.dot(H,invpos)
                    px=pos[0]+shift_a[0]+shift_b[0]++shift_c[0]
                    py=pos[1]+shift_a[1]+shift_b[1]++shift_c[1]
                    pz=pos[2]+shift_a[2]+shift_b[2]++shift_c[2]
                    
                    specie=datamol[i*Nat_per_mol+j,0]
                    if (conn[j,0] != -1):
                        c1=conn[j,0]+offset_conn
                    else:
                        c1=conn[j,0]
                    if (conn[j,1] != -1):
                        c2=conn[j,1]+offset_conn
                    else:
                        c2=conn[j,1]
                    if (conn[j,2] != -1):
                        c3=conn[j,2]+offset_conn
                    else:
                        c3=conn[j,2]
                    if (conn[j,3] != -1):
                        c4=conn[j,3]+offset_conn
                    else:
                        c4=conn[j,3]
                    
                    ofile.write('%s %s %d %d %5.8f %5.8f %5.8f %d %d %d %d\n' %(specie, specie.replace('TATB_', ''), molid, moltype, px, py, pz, c1, c2, c3, c4))
#                    ofile.write('%s %s %d %d %5.8f %5.8f %5.8f %d %d %d %d\n' %(specie, specie.replace('TATB_', ''), molid, moltype, px+40., py+40., pz+50., c1, c2, c3, c4))
#                    ofile.write('%s %5.8f %5.8f %5.8f\n' %(specie, px+40., py+40., pz+50.))                                                            
                molid+=1
                Nmol+=1
            Ncell+=1
#                    ofile.write('%s %5.8f %5.8f %5.8f %d\n' %(specie, px, py, pz, i+offset_idmol))
ofile.close()
