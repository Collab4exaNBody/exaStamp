#pragma once

#include <cmath>
#include <yaml-cpp/yaml.h>
#include <functional>
#include <vector>
#include <algorithm>
#include <string>
#include <mpi.h>

#include <onika/physics/units.h>
#include <onika/physics/units.h>
#include <exaStamp/particle_species/particle_specie.h>
//#include "exanb/container_utils.h"
#include <onika/log.h>
#include <onika/file_utils.h>

#include "lchbop_potpar.h"
#include "lchbop_FATP.h"

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

//#include "exanb/units/quantity_yaml.h"
//#include "exanb/units/unityConverterHelper.h"
//#include "exanb/particle_specie.h"
//#include "exanb/container_utils.h"
//#include "exanb/log.h"
//#include "exanb/file_utils.h"

//#include "lchbop_inits.h"
//#include "lchbop_struct.h"
//#include "lchbop_FATP.h"

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

using std::cin;


//=======================================================================//

struct VLCHBOP
{
  size_t newparfile;
//  int nspc,ispcc,ispch;
  size_t nspc,iispcc,iispch;
  int ispcc,ispch;
  double r2srm,rcmrm;
  double nelkimin=11.0/12.0,nelkibnd=5.0/6.0;

  RCUT pRC[2][2];
  VSR pVSR[2][2];
  VLR pVLR[2][2];
  VMR pVMR[2];
  GHc pGHc; 
  GHh pGHh; 
  FATP pFATP; 

  inline void initVCC(MPI_Comm);
  inline void initVHH(MPI_Comm);
  inline void initVCH(MPI_Comm);
  inline void printVLR();
  inline void printGc();
  inline void printGh();
  inline void printH();
};

//-----------------------------------------------------------------------//

inline void VLCHBOP::initVCC(MPI_Comm comm)
{
  const std::string vCC_filename = onika::data_file_path( "vCC.par" );

  newparfile=1;

  FILE *VCC;
  char texte[256];

  int err;
  int i1,i2,nfnnb;
  double r1,dr,drsq,drcu,drqu,drsx;
  double gc1[4],gc2[6],gc3[3];
  double y,Py,dPy,d2Py;
  double fdr,gdr,dfdr;
  double term1,dterm1,d2term1;
  double term2,dterm2,d2term2;
  double vsrr,dvsrr,d2vsrr;
  double vsra,dvsra,d2vsra;
  double xnij,xnji,xmijel,xmjiel,xmeleff,xmijelmin,xmijelmax,xmeleffmin;
  double vlr,dvlr,ddvlr;

  RCUT prc;
  VSR pvsr;
  VLR pvlr;
  VMR pvmr;

  // MPI Initialization
  int rank=0;
  MPI_Comm_rank(comm, &rank);

  if (rank==0) lout << "using vCC file "<<vCC_filename<<std::endl;
  VCC=fopen(vCC_filename.c_str(),"r");
  assert(VCC!=nullptr);

//C Reading SR parameters
  err=fscanf(VCC,"%s\n",texte);
  err=fscanf(VCC,"%s\n",texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&prc.r1sr,&prc.r2sr,&prc.p2z,texte);
  err=fscanf(VCC,"%lf%lf%s\n",&pvsr.Asr,&pvsr.alphasr,texte);			
//  printf("Asr :%20.12lf\n",pvsr.Asr);
  err=fscanf(VCC,"%lf%lf%s\n",&pvsr.Bsr1,&pvsr.betasr1,texte);	
  err=fscanf(VCC,"%lf%lf%s\n",&pvsr.Bsr2,&pvsr.betasr2,texte);
  err=fscanf(VCC,"%hhu%hhu%hhu%s\n",&pGHc.ngc[0],&pGHc.ngc[1],&pGHc.ngc[2],texte);
  err=fscanf(VCC,"%lf%lf%s\n",&gc1[0],&gc1[1],texte);
  err=fscanf(VCC,"%lf%lf%s\n",&gc1[2],&gc1[3],texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&gc2[0],&gc2[1],&gc2[2],texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&gc2[3],&gc2[4],&gc2[5],texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&gc3[0],&gc3[1],&gc3[2],texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&pGHc.zm,&pGHc.Ay0,&pGHc.By0,texte);	
  err=fscanf(VCC,"%lf%lf%lf%lf%s\n",&pGHc.Agzmax,&pGHc.Bgzmax,&pGHc.Cgzmax,&pGHc.Dgzmax,texte);
  err=fscanf(VCC,"%lf%lf%s\n",&pGHc.Egz2,&pGHc.Fgz2,texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&pGHc.du,&pGHc.C1H,&pGHc.C4H,texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&pGHc.Ak,&pGHc.Bk,&pGHc.Ck,texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&pFATP.cf1,&pFATP.p2sm1,&pFATP.p2sm2,texte);
  for (size_t nc=0;nc<2;nc++) 
  {
    err=fscanf(VCC,"%s\n",texte);
    for (size_t i2=0;i2<4;i2++)
    {
      err=fscanf(VCC,"%lf%lf%lf%lf\n",
      &pFATP.fcnj[i2][0][nc],&pFATP.fcnj[i2][1][nc],&pFATP.fcnj[i2][2][nc],&pFATP.fcnj[i2][3][nc]);
    }
  }
  err=fscanf(VCC,"%lf%lf%lf%lf%s\n",&pFATP.aab1,&pFATP.aab2[1][1],&pFATP.aab2[2][2],&pFATP.aab2[1][2],texte);
  err=fscanf(VCC,"%lf%lf%lf%lf%s\n",&pFATP.At10,&pFATP.At11,&pFATP.At20,&pFATP.At21,texte);
  err=fscanf(VCC,"%lf%lf%lf%lf%s\n",&pFATP.Bt1,&pFATP.Bt2,&pFATP.Bt3,&pFATP.Bt4,texte);

//C Reading LR parameters
  err=fscanf(VCC,"%s\n",texte);
  err=fscanf(VCC,"%s\n",texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&pvlr.r1vlr,&pvlr.r2vlr,&pvlr.rsvlr,texte);
  err=fscanf(VCC,"%lf%lf%s\n",&pvlr.epsLJ,&pvlr.sigLJ,texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&pvlr.C1lr,&pvlr.D1lr,&pvlr.E1lr,texte);

//C Reading MR parameters
  err=fscanf(VCC,"%s\n",texte);
  err=fscanf(VCC,"%s\n",texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&prc.rcmr0,&prc.rcmr,&pvsr.b230,texte);
  err=fscanf(VCC,"%lf%lf%lf%s\n",&pvmr.Amr1,&pvmr.Amr2,&pvmr.Amr3,texte);
  err=fscanf(VCC,"%lf%lf%s\n",&pvmr.Bmr1,&pvmr.Bmr2,texte);
  fclose(VCC);
//  printf("err :%6d\n",err);
//  printf("Parameters correctly read\n");
//  cin.get();

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C Initilisation of old potential parameters
  pGHc.Cy0=0.0;
  pvmr.Bmr3=0.0;
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

// Code for unit conversion
//  printf("Asr :%20.12lf\n",pvsr.Asr);
//  double cfe=1.0/1.03641882007443324881e-04;
//  printf("Asr :%20.12lf\n",cfe*pvsr.Asr);
#ifndef EnIneV
  pvsr.Asr=UnityConverterHelper::convert(pvsr.Asr,"eV");
  pvsr.Bsr1=UnityConverterHelper::convert(pvsr.Bsr1,"eV");
  pvsr.Bsr2=UnityConverterHelper::convert(pvsr.Bsr2,"eV");
  pvlr.epsLJ=UnityConverterHelper::convert(pvlr.epsLJ,"eV");
  pvlr.C1lr=UnityConverterHelper::convert(pvlr.C1lr,"eV");
  pvlr.D1lr=UnityConverterHelper::convert(pvlr.D1lr,"eV");
  pvlr.E1lr=UnityConverterHelper::convert(pvlr.E1lr,"eV");
#endif
//  printf("Asr :%20.12lf\n",pvsr.Asr);
//  prc.r1sr=UnityConverterHelper::convert(prc.r1sr,"ang");
//  prc.r2sr=UnityConverterHelper::convert(prc.r2sr,"ang");
//  pvsr.alphasr=UnityConverterHelper::convert(pvsr.alphasr,"1/ang");
//  pvsr.betasr1=UnityConverterHelper::convert(pvsr.betasr1,"1/ang");
//  pvsr.betasr2=UnityConverterHelper::convert(pvsr.betasr2,"1/ang");
//  pGHc.du=UnityConverterHelper::convert(pGHc.du,"ang");
//  pGHc.C1H=UnityConverterHelper::convert(pGHc.C1H,"1/ang");
//  pGHc.C4H=UnityConverterHelper::convert(pGHc.C4H,"1/ang^4");
//  pvlr.r1vlr=UnityConverterHelper::convert(pvlr.r1vlr,"ang");
//  pvlr.r2vlr=UnityConverterHelper::convert(pvlr.r2vlr,"ang");
//  pvlr.rsvlr=UnityConverterHelper::convert(pvlr.rsvlr,"ang");
//  pvlr.sigLJ=UnityConverterHelper::convert(pvlr.sigLJ,"ang");
//  pvlr.C1lr=UnityConverterHelper::convert(pvlr.C1lr,"1/ang^3");
//  pvlr.D1lr=UnityConverterHelper::convert(pvlr.D1lr,"1/ang^4");
//  pvlr.E1lr=UnityConverterHelper::convert(pvlr.E1lr,"1/ang^5");
//  prc.rcmr0=UnityConverterHelper::convert(prc.rcmr0,"ang");
//  prc.rcmr=UnityConverterHelper::convert(prc.rcmr,"ang");
//  cin.get();

//C parameters for cut-off radii and functions
  prc.r21inv=1./(prc.r2sr-prc.r1sr);
  prc.r1srsq=prc.r1sr*prc.r1sr;
  prc.r2srsq=prc.r2sr*prc.r2sr;
  if(prc.r2sr>r2srm)r2srm=prc.r2sr;
  if(prc.r2sr>r2srm)r2srm=prc.r2sr;
  
  if(prc.rcmr>prc.rcmr0){prc.rcmrmax=prc.rcmr;}
  else{prc.rcmrmax=prc.rcmr0;}
  if(prc.rcmrmax>rcmrm)rcmrm=prc.rcmrmax;
  prc.rcmrmaxsq=prc.rcmrmax*prc.rcmrmax;
  prc.drcmr=prc.rcmr-prc.rcmr0;

//C Parameters for the polynomial cut-off off the SMR potentials
  pvsr.r1sr=prc.r1sr;
  pvsr.r2sr=prc.r2sr;
  pvsr.r21inv=prc.r21inv;
  r1=prc.r1sr;
  vsrr=pvsr.Asr*exp(-pvsr.alphasr*r1);
  dvsrr=-pvsr.alphasr*vsrr;
  d2vsrr=-pvsr.alphasr*dvsrr;
  
  term1=pvsr.Bsr1*exp(-pvsr.betasr1*r1);
  term2=pvsr.Bsr2*exp(-pvsr.betasr2*r1);
  dterm1=-pvsr.betasr1*term1;
  dterm2=-pvsr.betasr2*term2;
  d2term1=-pvsr.betasr1*dterm1;
  d2term2=-pvsr.betasr2*dterm2;
  vsra=term1+term2;
  dvsra=dterm1+dterm2;
  d2vsra=d2term1+d2term2;
  
  dr=pvsr.r2sr-pvsr.r1sr;
  drsq=dr*dr;
  pvsr.a10=vsrr;
  pvsr.b10=vsra;
  pvsr.a11=dr*dvsrr+3.*vsrr;
  pvsr.b11=dr*dvsra+3.*vsra;
  pvsr.a12=0.5*drsq*d2vsrr-3.*(vsrr-pvsr.a11);
  pvsr.b12=0.5*drsq*d2vsra-3.*(vsra-pvsr.b11);
  pvsr.b200=vsra;
  pvsr.db200=dvsra;
  pvsr.d2b200=d2vsra;

//C Parameters for G(y,z)
  for(uint8_t i=0;i<pGHc.ngc[0];i++){pGHc.gc[0].push_back(gc1[i]);}
  for(uint8_t i=0;i<pGHc.ngc[1];i++){pGHc.gc[1].push_back(gc2[i]);}
  for(uint8_t i=0;i<pGHc.ngc[2];i++){pGHc.gc[2].push_back(gc3[i]);}

//C Determination of gmin by imposing continuity of G1 at y=-1/2
  y=-0.5;
  pGHc.PdPy(0,y,Py,dPy,d2Py);
  pGHc.gmin=-Py*(1.+y);
  pGHc.PdPy(1,y,Py,dPy,d2Py);
  pGHc.gmin=pGHc.gmin+Py;

//C determination of gmin by imposing continuity of G1 at y=-1/3
  y=-1./3.;
  pGHc.PdPy(1,y,Py,dPy,d2Py);
  pGHc.gmax=Py;
  pGHc.PdPy(2,y,Py,dPy,d2Py);
  pGHc.gmax=pGHc.gmax-Py*(1.-y)*(1.-y);

//C Parameters for H(u,z)
  dr=-pGHc.du;
  drsq=dr*dr;
  drqu=drsq*drsq;
  drsx=drqu*drsq;
  pGHc.C2H=0.5*pGHc.C1H*pGHc.C1H;
  pGHc.C6H=-(2.*pGHc.C2H+12.*pGHc.C4H*drsq)/(30.*drqu);
  
  fdr=1.+pGHc.C1H*dr+pGHc.C2H*drsq+pGHc.C4H*drqu+pGHc.C6H*drsx;
  gdr=2.*dr*(pGHc.C2H+2.*pGHc.C4H*drsq+3.*pGHc.C6H*drqu);
  dfdr=pGHc.C1H+gdr;
  pGHc.xLH=fdr;
  pGHc.xkappa=dfdr/pGHc.xLH;
  
  dr=pGHc.du;
  fdr=1.+pGHc.C1H*dr+pGHc.C2H*drsq+pGHc.C4H*drqu+pGHc.C6H*drsx;
  dfdr=pGHc.C1H-gdr;
  pGHc.R0H=fdr;
  pGHc.R1H=dfdr;

//C Symmetrization of anti-bonding
  pFATP.aab2[2][1]=pFATP.aab2[1][2];

//C matrix used in the conjugation routine
  pFATP.imat[0][1]=1;
  pFATP.imat[1][1]=1;
  pFATP.imat[0][2]=1;
  pFATP.imat[1][2]=2;
  pFATP.imat[2][2]=1;
  for(nfnnb=3;nfnnb<40;nfnnb++)
  {
    pFATP.imat[0][nfnnb]=1;
    pFATP.imat[1][nfnnb]=nfnnb;
    pFATP.imat[2][nfnnb]=nfnnb*(nfnnb-1)/2;
    pFATP.imat[3][nfnnb]=nfnnb*(nfnnb-1)*(nfnnb-2)/6;
  }

// tables used in the conjugation routine
  for(i1=0;i1<3;i1++)
  {
    xnij=i1;
    xmijelmin=4./(xnij+1.);
    xmijelmax=4.-xnij;
    pFATP.xmelmin1[i1]=xmijelmin;
    pFATP.dxmel1[i1]=xmijelmax-xmijelmin;
    for (i2=0;i2<3;i2++)
    {
      xnji=i2;
      xmijel=xmijelmin;
      xmjiel=4./(xnji+1.);
      xmeleff=xmijel+xmjiel;
      xmeleffmin=xmeleff;
      pFATP.xmelmin2[i1][i2]=xmeleff;
      xmijel=xmijelmax;
      xmjiel=4.-xnji;
      xmeleff=xmijel+xmjiel;
      pFATP.dxmel2[i1][i2]=xmeleff-xmeleffmin;
    }
  }
  pFATP.dxmel1[0]=1.;

//C Initialisations for the long range potential parameters 
//C for cut-off radii and functions
  pvlr.r2vlrsq=pvlr.r2vlr*pvlr.r2vlr;
  pvlr.dr12vlri=1./(pvlr.r2vlr-pvlr.r1vlr);

//C LJ/poly1 junction (at LJ minimum)
  pvlr.fepsLJ=4.*pvlr.epsLJ;
  pvlr.r0vlr=pow(2.,1./6.)*pvlr.sigLJ;
  pvlr.A1lr=-pvlr.epsLJ;
  pvlr.B1lr=144.*pow(2.,-4./3.)*pvlr.epsLJ/(pvlr.sigLJ*pvlr.sigLJ);

//C poly0/poly1 junction (at r2clr=2A)
  r1=pvlr.rsvlr;
  pvlr.rc0vlr=prc.r1sr;
  dr=r1-pvlr.rc0vlr;
  drsq=dr*dr;
  drcu=dr*drsq;
  pvlr.poly1(r1,vlr,dvlr,ddvlr);
  pvlr.A0lr=(vlr-r1*dvlr+3.*r1*vlr/dr)/drcu;
  pvlr.B0lr=(dvlr-3.*vlr/dr)/drcu;

//C LJ/poly2 junction (at r1vlr)
  r1=pvlr.r1vlr;
  dr=pvlr.r2vlr-pvlr.r1vlr;
  pvlr.dr12vlri=1./dr;
  pvlr.vlj(r1,vlr,dvlr,ddvlr);
  pvlr.A2lr=vlr;
  pvlr.B2lr=dr*dvlr+2.*pvlr.A2lr;
  pvlr.C2lr=0.5*dr*dr*ddvlr-pvlr.A2lr+2.*pvlr.B2lr;

//C Usefull parameter for compute MR interaction
  pvmr.imrnbm=2;
  pRC[ispcc][ispcc]=prc;
  pVSR[ispcc][ispcc]=pvsr;
  pVLR[ispcc][ispcc]=pvlr;
  pVMR[ispcc]=pvmr;
}

//=========================================================================//

inline void VLCHBOP::initVHH(MPI_Comm comm)
{
  const std::string vHH_filename = onika::data_file_path( "vHH.par" );

  FILE *VHH;
  char texte[256];

  int err;
  double r1,dr,drsq,drcu;
  double xnum,xden;
  double fdr,dfdr;
  double vsrr,dvsrr,d2vsrr;
  double vsra,dvsra,d2vsra;
  double vlr,dvlr,ddvlr;

  RCUT prc;
  VSR pvsr;
  VLR pvlr;
  VMR pvmr;

  // MPI Initialization
  int rank=0;
  MPI_Comm_rank(comm, &rank);

  if (rank==0) lout << "using vHH file "<<vHH_filename<<std::endl;
  VHH=fopen(vHH_filename.c_str(),"r");
// Reading SR parameters
  err=fscanf(VHH,"%s\n",texte);
  err=fscanf(VHH,"%s\n",texte);
  err=fscanf(VHH,"%lf%lf%lf%s\n",&prc.r1sr,&prc.r2sr,&prc.p2z,texte);
  err=fscanf(VHH,"%lf%lf%s\n",&pvsr.Asr,&pvsr.alphasr,texte);			
  err=fscanf(VHH,"%lf%lf%s\n",&pvsr.Bsr1,&pvsr.betasr1,texte);	

  err=fscanf(VHH,"%hhd%lf%s\n",&pGHh.ngh,&pGHh.ghmin,texte);	
  err=fscanf(VHH,"%lf%lf%lf%s\n",&pGHh.gh[0],&pGHh.gh[1],&pGHh.gh[2],texte);
  err=fscanf(VHH,"%lf%lf%lf%s\n",&pGHh.du,&pGHh.C1H,&pGHh.C2H,texte);
  err=fscanf(VHH,"%lf%lf%s\n",&pFATP.FcnjH2,&pFATP.pHH0,texte);

// Reading LR parameters
  err=fscanf(VHH,"%s\n",texte);
  err=fscanf(VHH,"%s\n",texte);
  err=fscanf(VHH,"%lf%lf%lf%s\n",&pvlr.r1vlr,&pvlr.r2vlr,&pvlr.rsvlr,texte);
  err=fscanf(VHH,"%lf%lf%s\n",&pvlr.epsLJ,&pvlr.sigLJ,texte);
  err=fscanf(VHH,"%lf%lf%lf%s\n",&pvlr.C1lr,&pvlr.D1lr,&pvlr.E1lr,texte);

// Reading MR parameters
  err=fscanf(VHH,"%s\n",texte);
  err=fscanf(VHH,"%s\n",texte);
  err=fscanf(VHH,"%lf%lf%lf%s\n",&prc.rcmr0,&prc.rcmr,&pvsr.b230,texte);
  err=fscanf(VHH,"%lf%lf%lf%s\n",&pvmr.Amr1,&pvmr.Amr2,&pvmr.Amr3,texte);
  fclose(VHH);
//  printf("err :%6d\n",err);
//  assert(err==0);

// Code for unit conversion
#ifndef EnIneV
  pvsr.Asr=UnityConverterHelper::convert(pvsr.Asr,"eV");
  pvsr.Bsr1=UnityConverterHelper::convert(pvsr.Bsr1,"eV");
  pvlr.epsLJ=UnityConverterHelper::convert(pvlr.epsLJ,"eV");
  pvlr.C1lr=UnityConverterHelper::convert(pvlr.C1lr,"eV");
  pvlr.D1lr=UnityConverterHelper::convert(pvlr.D1lr,"eV");
  pvlr.E1lr=UnityConverterHelper::convert(pvlr.E1lr,"eV");
#endif
//  prc.r1sr=UnityConverterHelper::convert(prc.r1sr,"ang");
//  prc.r2sr=UnityConverterHelper::convert(prc.r2sr,"ang");
//  pvsr.alphasr=UnityConverterHelper::convert(pvsr.alphasr,"1/ang");
//  pvsr.betasr1=UnityConverterHelper::convert(pvsr.betasr1,"1/ang");
//  pGHh.du=UnityConverterHelper::convert(pGHh.du,"ang");
//  pGHh.C1H=UnityConverterHelper::convert(pGHh.C1H,"1/ang");
//  pGHh.C2H=UnityConverterHelper::convert(pGHh.C2H,"1/ang^2");
//  pvlr.r1vlr=UnityConverterHelper::convert(pvlr.r1vlr,"ang");
//  pvlr.r2vlr=UnityConverterHelper::convert(pvlr.r2vlr,"ang");
//  pvlr.rsvlr=UnityConverterHelper::convert(pvlr.rsvlr,"ang");
//  pvlr.sigLJ=UnityConverterHelper::convert(pvlr.sigLJ,"ang");
//  pvlr.C1lr=UnityConverterHelper::convert(pvlr.C1lr,"1/ang^3");
//  pvlr.D1lr=UnityConverterHelper::convert(pvlr.D1lr,"1/ang^4");
//  pvlr.E1lr=UnityConverterHelper::convert(pvlr.E1lr,"1/ang^5");
//  prc.rcmr0=UnityConverterHelper::convert(prc.rcmr0,"ang");
//  prc.rcmr=UnityConverterHelper::convert(prc.rcmr,"ang");

//C Initilisation of old potential parameters
   pvmr.Bmr1=0.0;
   pvmr.Bmr2=0.0;
   pvmr.Bmr3=0.0;
  
//C parameters for cut-off radii and functions
  prc.r21inv=1./(prc.r2sr-prc.r1sr);
  prc.r1srsq=prc.r1sr*prc.r1sr;
  prc.r2srsq=prc.r2sr*prc.r2sr;
  if(prc.r2sr>r2srm)r2srm=prc.r2sr;
  
  if(prc.rcmr>prc.rcmr0){prc.rcmrmax=prc.rcmr;}
  else{prc.rcmrmax=prc.rcmr0;}
  prc.rcmrmaxsq=prc.rcmrmax*prc.rcmrmax;
  prc.drcmr=prc.rcmr-prc.rcmr0;

//C Parameters for the polynomial cut-off off the SMR potentials
  pvsr.r1sr=prc.r1sr;
  pvsr.r2sr=prc.r2sr;
  pvsr.r21inv=prc.r21inv;

  r1=prc.r1sr;
  vsrr=pvsr.Asr*exp(-pvsr.alphasr*r1);
  dvsrr=-pvsr.alphasr*vsrr;
  d2vsrr=-pvsr.alphasr*dvsrr;
  
  pvsr.Bsr2=0.;
  pvsr.betasr2=0.;
  vsra=pvsr.Bsr1*exp(-pvsr.betasr1*r1);
  dvsra=-pvsr.betasr1*vsra;
  d2vsra=-pvsr.betasr1*dvsra;
  
  dr=pvsr.r2sr-pvsr.r1sr;
  drsq=dr*dr;
  pvsr.a10=vsrr;
  pvsr.b10=vsra;
  pvsr.a11=dr*dvsrr+3.*vsrr;
  pvsr.b11=dr*dvsra+3.*vsra;
  pvsr.a12=0.5*drsq*d2vsrr-3.*(vsrr-pvsr.a11);
  pvsr.b12=0.5*drsq*d2vsra-3.*(vsra-pvsr.b11);
  pvsr.b200=vsra;
  pvsr.db200=dvsra;
  pvsr.d2b200=d2vsra;

//C Parameters for H(u,z)
  dr=-pGHh.du;
  drsq=dr*dr;
  xnum=-2.*pGHh.C2H*dr;
  xden=12.*drsq*dr;
  pGHh.C4H=xnum/xden;
  xnum=-pGHh.C2H-6.*pGHh.C4H*drsq;
  xden=3.*dr;
  pGHh.C3H=xnum/xden;
  
  dr=-pGHh.du;
  fdr=1.+dr*(pGHh.C1H+pGHh.C2H*dr+pGHh.C3H*drsq+pGHh.C4H*drsq*dr);
  dfdr=pGHh.C1H+dr*(2.*pGHh.C2H+3.*pGHh.C3H*dr+4.*pGHh.C4H*drsq);
  pGHh.xLH=fdr;
  pGHh.xkappa=dfdr/pGHh.xLH;
  
  dr=pGHh.du;
  fdr=1.+dr*(pGHh.C1H+pGHh.C2H*dr+pGHh.C3H*drsq+pGHh.C4H*drsq*dr);
  dfdr=pGHh.C1H+dr*(2.*pGHh.C2H+3.*pGHh.C3H*dr+4.*pGHh.C4H*drsq);
  pGHh.R0H=fdr;
  pGHh.R1H=dfdr;

//C Initialisations for the long range potential parameters for 
//C cut-off radii and functions
  pvlr.r2vlrsq=pvlr.r2vlr*pvlr.r2vlr;
  pvlr.dr12vlri=1./(pvlr.r2vlr-pvlr.r1vlr);

//C LJ/poly1 junction (at LJ minimum)
  pvlr.fepsLJ=4.*pvlr.epsLJ;
  pvlr.r0vlr=pow(2.,1./6.)*pvlr.sigLJ;
  pvlr.A1lr=-pvlr.epsLJ;
  pvlr.B1lr=144.*pow(2.,-4./3.)*pvlr.epsLJ/(pvlr.sigLJ*pvlr.sigLJ);

//C poly0/poly1 junction (at r2clr=2A)
  r1=pvlr.rsvlr;
  pvlr.rc0vlr=prc.r1sr;
  dr=r1-pvlr.rc0vlr;
  drsq=dr*dr;
  drcu=dr*drsq;
  pvlr.poly1(r1,vlr,dvlr,ddvlr);
  pvlr.A0lr=(vlr-r1*dvlr+3.*r1*vlr/dr)/drcu;
  pvlr.B0lr=(dvlr-3.*vlr/dr)/drcu;

//C LJ/poly2 junction (at r1vlr)
  r1=pvlr.r1vlr;
  dr=pvlr.r2vlr-pvlr.r1vlr;
  pvlr.dr12vlri=1./dr;
  pvlr.vlj(r1,vlr,dvlr,ddvlr);
  pvlr.A2lr=vlr;
  pvlr.B2lr=dr*dvlr+2.*pvlr.A2lr;
  pvlr.C2lr=0.5*dr*dr*ddvlr-pvlr.A2lr+2.*pvlr.B2lr;

// Usefull parameter for compute MR interaction
  pvmr.imrnbm=1;
  pRC[ispch][ispch]=prc;
  pVSR[ispch][ispch]=pvsr;
  pVLR[ispch][ispch]=pvlr;
  pVMR[ispch]=pvmr;
}

//=========================================================================//

inline void VLCHBOP::initVCH(MPI_Comm comm)
{
  const std::string vCH_filename = onika::data_file_path( "vCH.par" );

  FILE *VCH;
  char texte[256];

  int err;
  int i1,i2;
  double r1,dr,drsq,drcu;
  double vsrr,dvsrr,d2vsrr;
  double vsra,dvsra,d2vsra;
  double vlr,dvlr,ddvlr;

  RCUT prc;
  VSR pvsr;
  VLR pvlr;

  // MPI Initialization
  int rank=0;
  MPI_Comm_rank(comm, &rank);

  if (rank==0) lout << "using vCH file "<<vCH_filename<<std::endl;
  VCH=fopen(vCH_filename.c_str(),"r");
//C Reading SR parameters
  err=fscanf(VCH,"%s\n",texte);
  err=fscanf(VCH,"%s\n",texte);
  err=fscanf(VCH,"%lf%lf%lf%s\n",&prc.r1sr,&prc.r2sr,&prc.p2z,texte);
  err=fscanf(VCH,"%lf%lf%s\n",&pvsr.Asr,&pvsr.alphasr,texte);			
  err=fscanf(VCH,"%lf%lf%s\n",&pvsr.Bsr1,&pvsr.betasr1,texte);	
  err=fscanf(VCH,"%lf%lf%lf%lf%s\n",&pGHc.gC,&pGHc.gCbc0[ispcc][ispch],&pGHc.gCbc0[ispch][ispcc],&pGHc.gCbc0[ispch][ispch],texte);
  err=fscanf(VCH,"%lf%lf%lf%s\n",&pGHc.gCbc1[ispcc][ispch],&pGHc.gCbc1[ispch][ispcc],&pGHc.gCbc1[ispch][ispch],texte);
  err=fscanf(VCH,"%lf%s\n",&pFATP.phc0,texte);

  err=fscanf(VCH,"%s\n",texte);
  for(i1=0;i1<4;i1++)
  {
    err=fscanf(VCH,"%lf%lf%lf%lf\n",&pFATP.pcc[i1][0],&pFATP.pcc[i1][1],&pFATP.pcc[i1][2],&pFATP.pcc[i1][3]);
  }
  err=fscanf(VCH,"%s\n",texte);
  for(i1=0;i1<4;i1++)
  {
    err=fscanf(VCH,"%lf%lf%lf%lf\n",&pFATP.pch[i1][0],&pFATP.pch[i1][1],&pFATP.pch[i1][2],&pFATP.pch[i1][3]);
  }

//C Reading LR parameters
  err=fscanf(VCH,"%s\n",texte);
  err=fscanf(VCH,"%s\n",texte);
  err=fscanf(VCH,"%lf%lf%lf%s\n",&pvlr.r1vlr,&pvlr.r2vlr,&pvlr.rsvlr,texte);
  err=fscanf(VCH,"%lf%lf%s\n",&pvlr.epsLJ,&pvlr.sigLJ,texte);
  err=fscanf(VCH,"%lf%lf%lf%s\n",&pvlr.C1lr,&pvlr.D1lr,&pvlr.E1lr,texte);

//C Reading MR parameters
  err=fscanf(VCH,"%s\n",texte);
  err=fscanf(VCH,"%s\n",texte);
  err=fscanf(VCH,"%lf%lf%lf%s\n",&prc.rcmr0,&prc.rcmr,&pvsr.b230,texte);
  fclose(VCH);
//  printf("err :%6d\n",err);
//  assert(err==0);

// Code for unit conversion
#ifndef EnIneV
  pvsr.Asr=UnityConverterHelper::convert(pvsr.Asr,"eV");
  pvsr.Bsr1=UnityConverterHelper::convert(pvsr.Bsr1,"eV");
  pvlr.epsLJ=UnityConverterHelper::convert(pvlr.epsLJ,"eV");
  pvlr.C1lr=UnityConverterHelper::convert(pvlr.C1lr,"eV");
  pvlr.D1lr=UnityConverterHelper::convert(pvlr.D1lr,"eV");
  pvlr.E1lr=UnityConverterHelper::convert(pvlr.E1lr,"eV");
#endif
//  prc.r1sr=UnityConverterHelper::convert(prc.r1sr,"ang");
//  prc.r2sr=UnityConverterHelper::convert(prc.r2sr,"ang");
//  pvsr.alphasr=UnityConverterHelper::convert(pvsr.alphasr,"1/ang");
//  pvsr.betasr1=UnityConverterHelper::convert(pvsr.betasr1,"1/ang");
//  pvlr.r1vlr=UnityConverterHelper::convert(pvlr.r1vlr,"ang");
//  pvlr.r2vlr=UnityConverterHelper::convert(pvlr.r2vlr,"ang");
//  pvlr.rsvlr=UnityConverterHelper::convert(pvlr.rsvlr,"ang");
//  pvlr.sigLJ=UnityConverterHelper::convert(pvlr.sigLJ,"ang");
//  pvlr.C1lr=UnityConverterHelper::convert(pvlr.C1lr,"1/ang^3");
//  pvlr.D1lr=UnityConverterHelper::convert(pvlr.D1lr,"1/ang^4");
//  pvlr.E1lr=UnityConverterHelper::convert(pvlr.E1lr,"1/ang^5");
//  prc.rcmr0=UnityConverterHelper::convert(prc.rcmr0,"ang");
//  prc.rcmr=UnityConverterHelper::convert(prc.rcmr,"ang");

//C match Jan's MC code for pch
  pFATP.phc0*=0.5;
  for(i1=0;i1<4;i1++){
  for(i2=0;i2<4;i2++){pFATP.pch[i2][i1]*=0.5;}}

//iC Parameters for cut-off radii and functions
  prc.r21inv=1./(prc.r2sr-prc.r1sr);
  prc.r1srsq=prc.r1sr*prc.r1sr;
  prc.r2srsq=prc.r2sr*prc.r2sr;
  if(prc.r2sr>r2srm)r2srm=prc.r2sr;
  
  if(prc.rcmr>prc.rcmr0){prc.rcmrmax=prc.rcmr;}
  else{prc.rcmrmax=prc.rcmr0;}
  prc.rcmrmaxsq=prc.rcmrmax*prc.rcmrmax;
  prc.drcmr=prc.rcmr-prc.rcmr0;

//C Parameters for the polynomial cut-off off the SMR potentials
  pvsr.r1sr=prc.r1sr;
  pvsr.r2sr=prc.r2sr;
  pvsr.r21inv=prc.r21inv;

  r1=prc.r1sr;
  vsrr=pvsr.Asr*exp(-pvsr.alphasr*r1);
  dvsrr=-pvsr.alphasr*vsrr;
  d2vsrr=-pvsr.alphasr*dvsrr;
  
  pvsr.Bsr2=0.;
  pvsr.betasr2=0.;
  vsra=pvsr.Bsr1*exp(-pvsr.betasr1*r1);
  dvsra=-pvsr.betasr1*vsra;
  d2vsra=-pvsr.betasr1*dvsra;
  
  dr=pvsr.r2sr-pvsr.r1sr;
  drsq=dr*dr;
  pvsr.a10=vsrr;
  pvsr.b10=vsra;
  pvsr.a11=dr*dvsrr+3.*vsrr;
  pvsr.b11=dr*dvsra+3.*vsra;
  pvsr.a12=0.5*drsq*d2vsrr-3.*(vsrr-pvsr.a11);
  pvsr.b12=0.5*drsq*d2vsra-3.*(vsra-pvsr.b11);
  pvsr.b200=vsra;
  pvsr.db200=dvsra;
  pvsr.d2b200=d2vsra;

//C Initialisations for the Gyz function      
  pGHc.gCbc0[ispcc][ispcc]=1.;
  pGHc.gCbc1[ispcc][ispcc]=0.;

//C Parameters for H(u,z)
  pGHc.drCbc[ispcc][ispcc]=0.;
  pGHc.drCbc[ispch][ispch]=0.;
  pGHc.drCbc[ispch][ispcc]=-pGHc.drCbc[ispcc][ispch];
  pGHh.drHbc[ispcc][ispcc]=0.;
  pGHh.drHbc[ispch][ispch]=0.;
  pGHh.drHbc[ispch][ispcc]=-pGHh.drHbc[ispcc][ispch];

//C Initialisations for the long range potential parameters 
//C for cut-off radii and functions
  pvlr.r2vlrsq=pvlr.r2vlr*pvlr.r2vlr;
  pvlr.dr12vlri=1./(pvlr.r2vlr-pvlr.r1vlr);

//C LJ/poly1 junction (at LJ minimum)
  pvlr.fepsLJ=4.*pvlr.epsLJ;
  pvlr.r0vlr=pow(2.,1./6.)*pvlr.sigLJ;
  pvlr.A1lr=-pvlr.epsLJ;
  pvlr.B1lr=144.*pow(2.,-4./3.)*pvlr.epsLJ/(pvlr.sigLJ*pvlr.sigLJ);

//C poly0/poly1 junction (at r2clr=2A)
  r1=pvlr.rsvlr;
  pvlr.rc0vlr=prc.r1sr;
  dr=r1-pvlr.rc0vlr;
  drsq=dr*dr;
  drcu=dr*drsq;
  pvlr.poly1(r1,vlr,dvlr,ddvlr);
  pvlr.A0lr=(vlr-r1*dvlr+3.*r1*vlr/dr)/drcu;
  pvlr.B0lr=(dvlr-3.*vlr/dr)/drcu;

//C LJ/poly2 junction (at r1vlr)
  r1=pvlr.r1vlr;
  dr=pvlr.r2vlr-pvlr.r1vlr;
  pvlr.dr12vlri=1./dr;
  pvlr.vlj(r1,vlr,dvlr,ddvlr);
  pvlr.A2lr=vlr;
  pvlr.B2lr=dr*dvlr+2.*pvlr.A2lr;
  pvlr.C2lr=0.5*dr*dr*ddvlr-pvlr.A2lr+2.*pvlr.B2lr;

  pRC[ispcc][ispch]=prc;
  pRC[ispch][ispcc]=prc;
  pVSR[ispcc][ispch]=pvsr;
  pVSR[ispch][ispcc]=pvsr;
  pVLR[ispcc][ispch]=pvlr;
  pVLR[ispch][ispcc]=pvlr;
}

//=========================================================================//
//=========================================================================//

//-------------------------------------------------------------------------//

inline void VLCHBOP::printVLR()
{
  FILE *VLRDAT;
  size_t ir,nr=100,ispci=0,ispcj=0;
  double r1,r2,r,dr;
  double vlr,dvlr,ddvlr;

  VLRDAT=fopen("TMPDAT/vlr.dat","w");
  r1=pRC[ispci][ispcj].r1sr;
  r2=pVLR[ispci][ispcj].r2vlr;
  dr=(r2-r1)/nr;
  r=r1-dr;

  for(ir=0; ir <= nr; ir++)
  {
    r=r+dr;
    pVLR[ispci][ispcj].vlr(r,vlr,dvlr,ddvlr);
    fprintf(VLRDAT,"%16.8lf%16.8lf%16.8lf\n",r,vlr,dvlr);
//    fprintf(VLRDAT,"%16.8lf%16.8lf%16.8lf%16.8lf\n",r,vlr,dvlr,ddvlr);
  }
  fclose(VLRDAT);
}

//-----------------------------------------------------------------------//

inline void VLCHBOP::printGc()
{
  FILE *GcDAT;
  size_t iy,iz,ny=100,nz=6;
  double ymin=-1.0,ymax=1.0,y,dy,z;
  double Gyz,dGyzdy,dGyzdz;

  GcDAT=fopen("TMPDAT/Gc.dat","w");

  dy=(ymax-ymin)/ny;
  z=6.;

  for(iz=0;iz<=nz;iz++)
  {
    z=static_cast<double>(iz);
    y=ymin-dy;
    for(iy=0;iy<=ny;iy++)
    {
      y=y+dy;
      pGHc.GdGyz(y,z,Gyz,dGyzdy,dGyzdz);
      fprintf(GcDAT,"%16.8lf%16.8lf%16.8lf%16.8lf\n",y,Gyz,dGyzdy,dGyzdz);
    }
    fprintf(GcDAT,"\n");
  }  
  fclose(GcDAT);
}

//-----------------------------------------------------------------------//

inline void VLCHBOP::printGh()
{
  FILE *GhDAT;
  int iy,ny=100;
  double ymin=-1.0,ymax=1.0,y,dy,z=0.;
  double Gyz,dGyzdy,dGyzdz;

  GhDAT=fopen("TMPDAT/Gh.dat","w");

  dy=(ymax-ymin)/ny;
  y=ymin-dy;
  for(iy=0;iy<=ny;iy++)
  {
    y=y+dy;
    pGHh.GdGyz(y,z,Gyz,dGyzdy,dGyzdz);
    fprintf(GhDAT,"%16.8lf%16.8lf%16.8lf%16.8lf\n",y,Gyz,dGyzdy,dGyzdz);
  }
  fclose(GhDAT);
}

//-----------------------------------------------------------------------//

inline void VLCHBOP::printH()
{
  FILE *HDAT;
  size_t iu,nu=100;
  double umin=-2.0,umax=2.0,u,du;
  double Hdu,dHdu;

  HDAT=fopen("TMPDAT/H.dat","w");

  du=(umax-umin)/nu;
  u=umin-du;

  for(iu=0;iu<=nu;iu++)
  {
    u=u+du;
//    pGHc.HdHu(u,Hdu,dHdu);
    pGHh.HdHu(u,Hdu,dHdu);
    fprintf(HDAT,"%16.8lf%16.8lf%16.8lf\n",u,Hdu,dHdu);
  }
  fclose(HDAT);
}

//=========================================================================//

struct NbCells
{
  int iprcb,nnbc[4];
  IJK nbclst[172];
  inline void getnbclst();
};

//-----------------------------------------------------------------------//

inline void NbCells::getnbclst()
{
  IJK dims{7,7,7};
  ssize_t nijk=grid_cell_count(dims); 
  ssize_t iclc=(nijk-1)/2;
//  cout << "nijk,icl: " << nijk << "  " << icl << endl;
  IJK ijkc=grid_index_to_ijk(dims,iclc);
  IJK ijkmin=ijkc;
  IJK ijkmax=ijkc;
  IJK ijk0{0,0,0};

//  printf("ijkc   :%8zu%8zu%8zu\n",ijkc.i,ijkc.j,ijkc.k);
//  printf("ijkmin :%8zu%8zu%8zu\n",ijkmin.i,ijkmin.j,ijkmin.k);
//  printf("ijkmax :%8zu%8zu%8zu\n",ijkmax.i,ijkmax.j,ijkmax.k);
//  std::cin.get();
//  cout << "ijkmin  : " << ijkmin << endl;
//  cout << "ijkmax  : " << ijkmax << endl;

  int inbl=0,inbc=0;
  IJK dijk{0,0,0};
  nbclst[inbc++]={0,0,0};
  nnbc[inbl]=inbc;
  for(inbl=1;inbl<4;inbl++)
  {
    ijkmin=ijkmin-1; // operator - defined for IJK struct //
    ijkmax=ijkmax+1; // operator + defined for IJK struct //
//    printf("ijkmin :%8zd%8zd%8zd\n",ijkmin.i,ijkmin.j,ijkmin.k);
//    printf("ijkmax :%8zd%8zd%8zd\n",ijkmax.i,ijkmax.j,ijkmax.k);

    for(ssize_t jclk=ijkmin.k;jclk<ijkmax.k+1;jclk++)
    for(ssize_t jclj=ijkmin.j;jclj<ijkmax.j+1;jclj++)
    for(ssize_t jcli=ijkmin.i;jcli<ijkmax.i+1;jcli++)
    {
//      printf("jclijk:%6zd%6zd%6zd\n",jcli,jclj,jclk);
//      printf("ijkc :%6zd%6zd%6zd\n",ijkc.i,ijkc.j,ijkc.k);
//      printf("jclijk-ijkc :%6zd%6zd%6zd\n",jcli-ijkc.i,jclj-ijkc.j,jclk-ijkc.k);
      dijk={jcli-ijkc.i,jclj-ijkc.j,jclk-ijkc.k};
//      printf("dijk :%6zd%6zd%6zd\n",dijk.i,dijk.j,dijk.k);
      if(static_cast<int>(abs(dijk.i))<inbl&& 
         static_cast<int>(abs(dijk.j))<inbl&&
         static_cast<int>(abs(dijk.k))<inbl)continue;
      if(ijk0<dijk)continue;
      nbclst[inbc++]=dijk;
//      printf("inbc,nbclst :%6d,%6zd%6zd%6zd\n",inbc-1,nbclst[i].i,nbclst[i].j,nbclst[i].k);
//      cin.get();
    }  
    nnbc[inbl]=inbc;
  }  
}

//=========================================================================//
