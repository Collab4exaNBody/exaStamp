import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.labelsize']=25
mpl.rcParams['xtick.labelsize']=20
mpl.rcParams['ytick.labelsize']=20
mpl.rcParams['legend.fontsize']=20
mpl.rcParams['axes.titlesize']=24
mpl.rcParams['axes.linewidth']=2.0

def joules_to_ev(joules):
    """Convert energy from joules to electronvolts (eV)."""
    eV_per_joule = 1 / 1.6021892e-19  # exact CODATA value
    return joules * eV_per_joule

class EAMSuttonChenParameters:
    def __init__(self,rcut,c,epsilon,a0,n,m,Nrho,drho,Nr,dr,Zmat,mass,a0bis,lattice,nmats,specy):
        # EAM potential parameters
        self.rcut    = rcut
        self.c       = c
        self.epsilon = epsilon
        self.a0      = a0
        self.n       = n
        self.m       = m

        # Spline + eam/alloy format parameters
        self.Nrho=Nrho
        self.drho=drho
        self.Nr=Nr
        self.dr=dr
        self.Zmat=Zmat
        self.mass=mass
        self.a0bis=a0bis
        self.lattice=lattice
        self.nmats=nmats
        self.specy=specy
        
    def generate_file(self,filename):
        ofile=open(filename,'w')
        ofile.write('DATE: 2024-04-25 CONTRIBUTOR: Paul Lafourcade\n')
        ofile.write('Generated from read_eam_file.py\n')
        ofile.write('XX\n')
        ofile.write('%d %s\n' %(self.nmats,self.specy))
        ofile.write('%d %5.12e %d %5.12e %5.12e\n' %(self.Nrho,self.drho,self.Nr,self.dr,self.rcut))
        ofile.write('%d %5.4f %5.4f %s\n' %(self.Zmat,self.mass,self.a0,self.lattice))
        
        rho_x=np.linspace(0.,self.Nrho*self.drho,self.Nrho)
        r_x=np.linspace(0.,self.rcut,self.Nr)
        
        Frho=np.nan_to_num(self.Frho(rho_x),0.)
        for i in range(len(Frho)):
            ofile.write('%5.16e\n' %(Frho[i]))
            
        rhor=np.nan_to_num(self.rho(r_x),0.)
        for j in range(len(rhor)):
            ofile.write('%5.16e\n' %(rhor[j]))
            
        rphir=np.nan_to_num(self.rphi(r_x),0.)
        for k in range(len(rphir)):
            ofile.write('%5.16e\n' %(rphir[k]))
            
        ofile.close()
        
    def Frho(self, rho):
        sqrtRho = np.sqrt(rho)
        f = -1. * self.c * self.epsilon * sqrtRho
        return f
    
    def rho(self, r):
        ratio = self.a0/r
        rhoValue = np.power(ratio,self.m)
        return rhoValue
    
    def rphi(self, r):
        ratio = self.a0/r
        phiValue = self.epsilon * np.power(ratio,self.n)
        return r*phiValue

nmats=1
specy='Cu'
lattice='fcc'
Zmat=29
mass=63.546
a0bis=5.410

rcut = 7.29
c    = 33.17
epsilon = joules_to_ev(3.605e-21)
a0 = 3.27
n = 9.050
m = 5.005

Nrho=5000
rhomax=200.
drho=rhomax/Nrho

Nr=5000
rcutmax=rcut
dr=rcutmax/Nr

rho_x=np.linspace(0,Nrho*drho,Nrho)
r_x=np.linspace(0,Nr*dr,Nr)

EAM_SuttonChen = EAMSuttonChenParameters(rcut,
                                         c,
                                         epsilon,
                                         a0,
                                         n,
                                         m,
                                         Nrho,
                                         drho,
                                         Nr,
                                         dr,
                                         Zmat,
                                         mass,
                                         a0,
                                         lattice,
                                         nmats,
                                         specy)

EAM_SuttonChen.generate_file("Cu.eam.alloy")
