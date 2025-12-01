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

class EAMVniitfParameters:
    def __init__(self,rcut,rmax,rmin,rt0,Ecoh,E0,beta,A,Z,n,alpha,D,eta,mu,Nrho,drho,Nr,dr,Zmat,mass,a0,lattice,nmats,specy):
        # EAM potential parameters
        self.rcut  = rcut
        self.rmax  = rmax
        self.rmin  = rmin
        self.rt0   = rt0
        self.Ecoh  = Ecoh
        self.E0    = E0
        self.beta  = beta
        self.A     = A
        self.Z     = Z
        self.n     = n
        self.alpha = alpha
        self.D     = D
        self.eta   = eta
        self.mu    = mu

        # Spline + eam/alloy format parameters
        self.Nrho=Nrho
        self.drho=drho
        self.Nr=Nr
        self.dr=dr
        self.Zmat=Zmat
        self.mass=mass
        self.a0=a0
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
        a = np.power(rho, self.n);
        b = self.A*self.Ecoh * a;
        c = np.log(a);
        f  = b * (c-1.);
        f[np.where(rho<0.)]=0.
        return f

    def Spline_S3(self, x):
        x2 = x*x;
        return x2*x2 * ( -20*x2*x + 70*x2 - 84*x + 35 )
    
    def rho(self, r):
        irt0 = 1/self.rt0;

        F  = np.exp(-self.beta * (r*irt0 - 1.0) ) / self.Z;

        drS = (self.rmax-r)/(self.rmax-self.rmin);
        S  = self.Spline_S3(drS);
        
        S[np.where(drS<0)]=0.
        S[np.where(drS>1)]=1.

        return F*S
    
    def rphi(self, r):
        ir   = 1/r;
        irt0 = 1/self.rt0;
        dr   = r*irt0 - 1.0;
        dr2  = dr*dr;
        a    = -2*self.Ecoh/self.Z;
        b    = self.alpha*self.alpha*self.alpha*self.D*self.rt0;
        
        f1 = a*( 1 + self.alpha*dr + self.eta*dr2 + (self.mu+b*ir)*dr2*dr )
        f2 = np.exp(-self.alpha*dr)

        drS = (self.rmax-r)/(self.rmax-self.rmin)
        
        S = self.Spline_S3(drS)
        S[np.where(drS<0)]=0.
        S[np.where(drS>1)]=1.
        phi = (self.E0 + f1*f2) * S
        return r*phi

nmats=1
specy='Sn'
lattice='bcc'
Zmat=50
mass=118.69
a0=5.410

rcut = 8.0
rmax = 5.599
rmin = 1.000
rt0 = 3.437

Ecoh = 1.8450094310887324
E0 = 0.3214401234364775

beta = 6.000
A = 1.401
Z = 7.618
n = 0.724
alpha = 3.072
D = 0.1450
eta = 2.720
mu = -1.870

Nrho=5000
rhomax=10.
drho=rhomax/Nrho

Nr=5000
rcutmax=rcut
dr=rcutmax/Nr

rho_x=np.linspace(0,Nrho*drho,Nrho)
r_x=np.linspace(0,Nr*dr,Nr)

EAM_Sapozhnikov = EAMVniitfParameters(rcut,
                                      rmax,
                                      rmin,
                                      rt0,
                                      Ecoh,
                                      E0,
                                      beta,
                                      A,
                                      Z,
                                      n,
                                      alpha,
                                      D,
                                      eta,
                                      mu,
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

fig,ax=plt.subplots(1,4,figsize=(25,9))
ax[0].plot(rho_x,EAM_Sapozhnikov.Frho(rho_x),'r--',linewidth=3.0)
ax[1].plot(r_x,EAM_Sapozhnikov.rho(r_x),'r--',linewidth=3.0)
ax[2].plot(r_x,EAM_Sapozhnikov.rphi(r_x),'r--',linewidth=3.0)
ax[3].plot(r_x,(rmax-r_x)/(rmax-rmin),'r--',linewidth=3.0)

ax[0].set_xlabel(r'$\rho$')
ax[0].set_ylabel(r'$F(\rho)$')

ax[1].set_xlabel(r'r (ang)')
ax[1].set_ylabel(r'$\rho(r)$')

ax[2].set_xlabel(r'r (ang)')
ax[2].set_ylabel(r'r$\phi(r)$')

plt.savefig('EAM_Sapozhnikov.png')

EAM_Sapozhnikov.generate_file('Sn.eam.alloy')
