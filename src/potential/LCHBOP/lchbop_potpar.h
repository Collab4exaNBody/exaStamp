/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

//#include "exanb/container_utils.h"
#include <onika/log.h>

#include "lchbop_getEF.h"
//#include "lchbop_utils.h"


using std::pow;
//using std::pow(x,4) (doule) -ffast-math;
//using std::pow(x,6) (doule) -ffast-math;
using std::cin;
using std::cout;
using std::endl;

void sqdsqdown(double,double,double,double,double&,double&);
void sqdsqup(double,double,double,double,double&,double&);

//=======================================================================//

struct RCUT
{
  double r1sr,r2sr,p2z;
  double rcmr0,rcmr;
  double rcmrmax,rcmrmaxsq;
  double r21inv,r1srsq,r2srsq; 
  double drcmr,r2vlrsq; 
  inline void rcmrij(int,int,int,double,double&,double&);
};

//---------------------------------------------------------------------//

void RCUT::rcmrij(int ispci,int ispcj,int ispcc,double ndb,
double& rcmr,double& drcmrdndb)
{
  double sndb,dsndb;
  if(ispci==ispcc&&ispcj==ispcc)
  {
    sqdsqdown(ndb,1.,2.,1.,sndb,dsndb);
    rcmr=rcmr0+sndb*drcmr;
    drcmrdndb=dsndb*drcmr;
  }  
  else{rcmr=rcmr0;drcmrdndb=0.0;}
}

//=======================================================================//

struct VSR
{
  double r1sr,r2sr,r21inv;
  double Asr,alphasr;
  double Bsr1,betasr1;
  double Bsr2,betasr2;
  double b230;
  double a10,a11,a12;
  double b10,b11,b12;
  double b200,db200,d2b200;
  inline void vdvsrra(double,double,double,
  double&,double&,double&,double&,double&,double&);
  inline void vdvmra(double,double,double&,double&,double&);
};

//-----------------------------------------------------------------------//

void VSR::vdvsrra(double rij,double rceff,double smrloc,
double& vsrr,double& dvsrrdr,double& vsmra,double& dvsmradr,
double& dvsmradrc,double& dvsmradsmr)
{      
  if(rij<=r1sr)
  {
//C printf("vdvsrra; pure \n");
    vsrr=Asr*exp(-alphasr*rij);
    dvsrrdr=-alphasr*vsrr;
    double term1=Bsr1*exp(-betasr1*rij);
    double term2=Bsr2*exp(-betasr2*rij);
    vsmra=term1+term2;
    dvsmradr=-betasr1*term1-betasr2*term2;
    dvsmradrc=0.;
    dvsmradsmr=0.;
  }
  else
  {
    if(smrloc==0.)
    {
      double q=r21inv*(rij-r1sr);
      double qsq=q*q;
      double onemq=1.-q;
      double onemqsq=onemq*onemq;
      double onemqcu=onemq*onemqsq;
      double donemqcu=-3.*onemqsq;
      double frq=a10+a11*q+a12*qsq;
      double dfrq=a11+2.*a12*q;
      double faq=b10+b11*q+b12*qsq;
      double dfaq=(b11+2.*b12*q);

      vsrr=frq*onemqcu;
      dvsrrdr=(donemqcu*frq+onemqcu*dfrq)*r21inv;
      vsmra=faq*onemqcu;
      dvsmradr=(donemqcu*faq+onemqcu*dfaq)*r21inv;
      dvsmradrc=0.;
      dvsmradsmr=0.;
    }  
    else if(smrloc==1.)
    {
      if(rij<=r2sr)
      {
//C printf("vdvsrra; frac full mr\n");
        double q=r21inv*(rij-r1sr);
        double qsq=q*q;
        double onemq=1.-q;
        double onemqsq=onemq*onemq;
        double onemqcu=onemq*onemqsq;
        double donemqcu=-3.*onemqsq;
        double frq=a10+a11*q+a12*qsq;
        double dfrq=a11+2.*a12*q;
        double faq=b10+b11*q+b12*qsq;
        double vsmra0=faq*onemqcu;

        vsrr=frq*onemqcu;
        dvsrrdr=(donemqcu*frq+onemqcu*dfrq)*r21inv;
        vdvmra(rij,rceff,vsmra,dvsmradr,dvsmradrc);
        dvsmradsmr=vsmra-vsmra0;
      }
      else
      {
//C printf("vdvsrra; pure full mr\n");
        vsrr=0.;
        dvsrrdr=0.;
        vdvmra(rij,rceff,vsmra,dvsmradr,dvsmradrc);
        dvsmradsmr=vsmra;
      }
    }  
    else
    { 
      if(rij<=r2sr)
      {
//C printf("vdvsrra; frac mr\n");
        double q=r21inv*(rij-r1sr);
        double qsq=q*q;
        double onemq=1.-q;
        double onemqsq=onemq*onemq;
        double onemqcu=onemq*onemqsq;
        double donemqcu=-3.*onemqsq;
        double frq=a10+a11*q+a12*qsq;
        double dfrq=a11+2.*a12*q;
        double faq=b10+b11*q+b12*qsq;
        double dfaq=(b11+2.*b12*q);
        double vsmra0=faq*onemqcu;
        double dvsmra0dr=(donemqcu*faq+onemqcu*dfaq)*r21inv;

        vsrr=frq*onemqcu;
        dvsrrdr=(donemqcu*frq+onemqcu*dfrq)*r21inv;
        vdvmra(rij,rceff,vsmra,dvsmradr,dvsmradrc);
        dvsmradrc=smrloc*dvsmradrc;
        dvsmradsmr=vsmra-vsmra0;
        vsmra=vsmra0+smrloc*(vsmra-vsmra0);
        dvsmradr=dvsmra0dr+smrloc*(dvsmradr-dvsmra0dr);
      } 
      else
      {
//C printf("vdvsrra; pure mr\n");
        vsrr=0.;
        dvsrrdr=0.;
        vdvmra(rij,rceff,vsmra,dvsmradr,dvsmradrc);
  
        dvsmradrc=smrloc*dvsmradrc;
        dvsmradsmr=vsmra;
        vsmra=smrloc*vsmra;
        dvsmradr=smrloc*dvsmradr;
      }
    }
  }
}
//----------------------------------------------------------------------//

void VSR::vdvmra(double rij,double rceff,
double& vsmra, double& dvsmradr,double& dvsmradrc)
{      
//C local data
  double b20,b21,b22,b23,db20,d2b20,db21drc,db22drc;
  double rcmr1,rcmr1sq;
  double q,qsq,dqdr,dqdrc;
  double onemq,onemqsq,onemqcu,donemqcudr,donemqcudrc;
  double faq,dfaqdq,dfaqdr,dfaqdrc;
  double faqnum,dfaqnumdq,dfaqnumdrc;
  double faqden,dfaqdendq;
  
  rcmr1=rceff-r1sr;
  rcmr1sq=rcmr1*rcmr1;
  
  b20=b200;
  db20=db200;
  d2b20=d2b200;
  b23=b230;
  
  b21=rcmr1*db20+(3.+b23)*b20;
  db21drc=db20;
  
  b22=0.5*rcmr1sq*d2b20-3.*(b20-b21)-(3.*b20-b21)*b23-b20*b23*b23;
  db22drc=rcmr1*d2b20+(3.+b23)*db21drc;
  
  dqdr=1./rcmr1;
  q=(rij-r1sr)*dqdr;
  dqdrc=-(rij-r1sr)/rcmr1sq;
  
  qsq=q*q;
  faqnum=b20+b21*q+b22*qsq;
  dfaqnumdq=b21+2.*b22*q;
  dfaqnumdrc=db21drc*q+db22drc*qsq;
  faqden=1.+b23*q;
  dfaqdendq=b23;
  faq=faqnum/faqden;
  dfaqdq=(dfaqnumdq*faqden-faqnum*dfaqdendq)/(faqden*faqden);
  dfaqdr=dfaqdq*dqdr;
  dfaqdrc=dfaqdq*dqdrc+dfaqnumdrc/faqden;
  
  onemq=1.-q;
  onemqsq=onemq*onemq;
  onemqcu=onemq*onemqsq;
  donemqcudr=-3.*onemqsq*dqdr;
  donemqcudrc=-3.*onemqsq*dqdrc;
  
  vsmra=faq*onemqcu;
  dvsmradr=donemqcudr*faq+onemqcu*dfaqdr;
  dvsmradrc=donemqcudrc*faq+onemqcu*dfaqdrc;
}

//=======================================================================//

struct GHc
{
  size_t newparfile=1;
  uint8_t ngc[3];
  std::vector<double> gc[3];
  double oneby3=1.0/3.0;
  double gC,gCbc0[2][2],gCbc1[2][2];
  double zm,Ay0,By0,Cy0;
  double Agzmax,Bgzmax,Cgzmax,Dgzmax,Egz2,Fgz2;
  double gmin,gmax;
  double du,xkappa,xLH,drCbc[2][2];
  double C1H,C2H,C4H,C6H;
  double R0H,R1H;
  double Ak,Bk,Ck;
  inline void GdGyz(double,double,double&,double&,double&);
  inline void G1dG1y(double,double&,double&,double&);
  inline void PdPy(uint8_t,double,double&,double&,double&);
  inline void HdHu(double,double&,double&);
};

//-----------------------------------------------------------------------//

void GHc::GdGyz(double y,double z,double& Gyz,double& dGyzdy,double& dGyzdz)
{
  double G1,dG1dy,d2G1dy;
  
  if(z>zm)z=zm;
  double zmzm=z-zm;
#ifdef PowRepl
  double zmzmsq=zmzm*zmzm;
  double zmzmqu=zmzmsq*zmzmsq;
  double y0=Ay0*zmzm+By0*zmzmqu;
#else
  double y0=Ay0*zmzm+By0*pow(zmzm,4);
#endif
  
  if (y<y0) 
  {
    G1dG1y(y,G1,dG1dy,d2G1dy);
    Gyz=G1;
    dGyzdy=dG1dy;
    dGyzdz=0.;
  } 
  else 
  {
    double gzmax,dgzmax,gz0,gz1,gz2,dgz0,dgz1,dgz2;
    if(z<zm)
    {
      double zsq=z*z;
#ifdef PowRepl
      double zmzmcu=zmzmsq*zmzm;
#else
      double zmzmcu=zmzm*zmzm*zmzm;
#endif
      double y0mone=y0-1.;
      double y0monesq=y0mone*y0mone;
      double y0monecu=y0monesq*y0mone;
      double y0monequ=y0monecu*y0mone;
      double y0sq=y0*y0;
      double dy0=Ay0+4.0*By0*zmzmcu;
      
      G1dG1y(y0,G1,dG1dy,d2G1dy);
      
      gzmax=gmax-(Agzmax+Bgzmax*z+Cgzmax*zsq+Dgzmax*zsq*zsq)*y0sq;
      dgzmax=-(Bgzmax+2.*Cgzmax*z+4.*Dgzmax*zsq*z)*y0sq-
             2.*(Agzmax+Bgzmax*z+Cgzmax*zsq+Dgzmax*zsq*zsq)*y0*dy0;
      
      double hfz=Egz2*(36.0+12.0*z+Fgz2*zsq); 
      double dhfz=Egz2*(12.0+2.0*Fgz2*z);
#ifndef PowRepl
      double zmzmsq=zmzm*zmzm;
#endif

      gz2=gc[2][2]+hfz*zmzmsq;
      dgz2=dhfz*zmzmsq+2.0*hfz*zmzm;
      gz1=dG1dy/y0monesq-2.*(G1-gzmax)/y0monecu-2.*gz2*y0;
      dgz1=(d2G1dy/y0monesq-4.*dG1dy/y0monecu+6.*(G1-gzmax)/y0monequ-2.*gz2)*dy0+2.*dgzmax/y0monecu-2.*y0*dgz2;
      gz0=(G1-gzmax)/y0monesq-gz1*y0-gz2*y0sq;
      dgz0=( dG1dy/y0monesq-2.*(G1-gzmax)/y0monecu-gz1-2.*y0*gz2 )*dy0-dgzmax/y0monesq-y0*dgz1-y0sq*dgz2;

      double onemy=1.-y;
      double onemysq=onemy*onemy;
      double Py=gz0+gz1*y+gz2*y*y;
      double dPydy=gz1+2.*gz2*y;
      double dPydz=dgz0+dgz1*y+dgz2*y*y;
      
      Gyz=gzmax+onemysq*Py;
      dGyzdy=-2.*onemy*Py+onemysq*dPydy;
      dGyzdz=dgzmax+onemysq*dPydz;
    } 
    else 
    {
      G1dG1y(y0,G1,dG1dy,d2G1dy);
      
      gzmax=gmax;
      gz2=gc[2][2];
      gz1=dG1dy+2.*(G1-gzmax);
      gz0=(G1-gzmax);

      double onemy=1.-y;
      double onemysq=onemy*onemy;
      double Py=gz0+gz1*y+gz2*y*y;
      double dPydy=gz1+2.*gz2*y;
      
      Gyz=gzmax+onemysq*Py;
      dGyzdy=-2.*onemy*Py+onemysq*dPydy;
      dGyzdz=0.;
    }
  }

//printf("Gyz :%16.8lf%16.8lf%16.8lf\n",y,z,Gyz);
//cin.get();

}

//-----------------------------------------------------------------------//

void GHc::G1dG1y(double y,double& G1,double& d1G1,double& d2G1)
{
  double Py,dPy,d2Py;

  if(y<=-0.5)
  {
    PdPy(0,y,Py,dPy,d2Py);
    double ypone=y+1.;
    G1=gmin+ypone*Py;
    d1G1=Py+ypone*dPy;
    d2G1=2.0*dPy+ypone*d2Py;
  } 
  else if (y<=-oneby3) 
  {
    PdPy(1,y,Py,dPy,d2Py);
    G1=Py;
    d1G1=dPy;
    d2G1=d2Py;
  } 
  else 
  {
    PdPy(2,y,Py,dPy,d2Py);
    double ymone=y-1.;
    double ymonesq=ymone*ymone;
    G1=gmax+ymonesq*Py;
    d1G1=2.*ymone*Py+ymonesq*dPy;
    d2G1=2.*Py+4.*ymone*dPy+ymonesq*d2Py;
  }
}

//-----------------------------------------------------------------------//

//CC Polynome P(y) and dP(y) needed for G_C(y,z) and dG_C(y,z)
void GHc::PdPy(uint8_t igc, double y, double& Py, double& dPy, double& d2Py)
{
  int i1;
  double yn,ynp1;
  
  yn=1.;
  ynp1=y;
  Py=gc[igc][0]+gc[igc][1]*ynp1;
  dPy=gc[igc][1];
  d2Py=0.;
  for (i1=2;i1<ngc[igc];i1++)
  {
    d2Py+=(i1-1)*i1*gc[igc][i1]*yn;
    yn=ynp1;
    dPy+=i1*gc[igc][i1]*yn;
    ynp1=yn*y;
    Py+=gc[igc][i1]*ynp1;
  }
}

//-----------------------------------------------------------------------//

void GHc::HdHu(double u,double& Hdu,double& dHdu)
{
//  double usq,ucu,uqu;
//  double y,yqu,f1dr,f2dr;
  
  if(u<-du) 
  {
    double y=xkappa*(u+du);
#ifdef PowRepl
    double yqu=y*y;
    yqu*=yqu;
#else
    double yqu=pow(y,4);
#endif
    double f1dr=1./(1.+yqu);
    double f2dr=pow(f1dr,0.25);
    Hdu=xLH*(1.+y*f2dr);
    dHdu=xLH*xkappa*f1dr*f2dr;
  } 
  else if(u<du) 
  {
    double usq=u*u;
    double ucu=usq*u;
    double uqu=usq*usq;
    Hdu=1.+C1H*u+C2H*usq+C4H*uqu+C6H*uqu*usq;
    dHdu=C1H+2.*C2H*u+4.*C4H*ucu+6.*C6H*uqu*u;
  } 
  else 
  {
    Hdu=R0H+R1H*(u-du);
    dHdu=R1H;
//    printf("term1,term2 :%18.10lf%18.10lf\n",R0H,R1H*(u-du));
//    printf("u,du        :%18.10lf%18.10lf\n",u,du);
  }
}

//=======================================================================//

struct GHh
{
  uint8_t ngh;
  double ghmin,gh[3];
  double du,xkappa,xLH,drHbc[2][2];
  double C1H,C2H,C3H,C4H,C6H;
  double R0H,R1H;
  inline void GdGyz(double,double,double&,double&,double&);
  inline void HdHu(double,double&,double&);
};

//-----------------------------------------------------------------------//

void GHh::GdGyz(double y,double z,double& Gyz,double& dGyzdy,double& dGyzdz)
{
  double yn=1.;
  double ynp1=y;
  double Py=gh[0]+gh[1]*ynp1;
  double dPy=gh[1];
  double d2Py=0.;
  for(int i1=2;i1<ngh;i1++)
  {
    d2Py+=(i1-1)*i1*gh[i1]*yn;
    yn=ynp1;
    dPy+=i1*gh[i1]*yn;
    ynp1=yn*y;
    Py+=gh[i1]*ynp1;
  }
  double ypone=y+1.;
  double yponesq=ypone*ypone;

  Gyz=ghmin+yponesq*Py;
  dGyzdy=2.*ypone*Py+yponesq*dPy;
  dGyzdz=0.;
}
//-----------------------------------------------------------------------//

void GHh::HdHu(double u,double& Hdu,double& dHdu)
{
  if(u<-du) 
  {
    double y=xkappa*(u+du);
#ifdef PowRepl
    double yqu=y*y;
    yqu*=yqu;
#else
    double yqu=pow(y,4);
#endif
    double f1dr=1./(1.+yqu);
    double f2dr=pow(f1dr,0.25);
    Hdu=xLH*(1.+y*f2dr);
    dHdu=xLH*xkappa*f1dr*f2dr;
  } 
  else if(u<du) 
  {
    double usq=u*u;
    double ucu=usq*u;
    double uqu=usq*usq;
    Hdu=1.+C1H*u+C2H*usq+C3H*ucu+C4H*uqu;
    dHdu=C1H+2.*C2H*u+3.*C3H*usq+4.*C4H*ucu;
  } 
  else 
  {
    Hdu=R0H+R1H*(u-du);
    dHdu=R1H;
  }
}

//=======================================================================//

struct VLR
{
  double r1vlr,r2vlr,rsvlr;
  double epsLJ,sigLJ;
  double C1lr,D1lr,E1lr;
  double r2vlrsq,rc0vlr,r0vlr,dr12vlri,fepsLJ;
  double A0lr,B0lr,A1lr,B1lr,A2lr,B2lr,C2lr;
  inline void vlr  (double,double&,double&,double&);
  inline void poly0(double,double&,double&,double&);
  inline void poly1(double,double&,double&,double&);
  inline void vlj  (double,double&,double&,double&);
  inline void poly2(double,double&,double&,double&);
};

//-------------------------------------------------------------------------//

//C Compute LR energy and derivatives
void VLR::vlr(double rij,double& vlr,double& dvlr,double& ddvlr)
{
  if     (rij<rsvlr){poly0(rij,vlr,dvlr,ddvlr);}
  else if(rij<r0vlr){poly1(rij,vlr,dvlr,ddvlr);}
  else if(rij<r1vlr){vlj  (rij,vlr,dvlr,ddvlr);}
  else              {poly2(rij,vlr,dvlr,ddvlr);}
}

//-------------------------------------------------------------------------//

void VLR::poly0(double rij, double& vlr, double& dvlr, double& ddvlr)
{
  double rmrc0,rmrc0sq,rmrc0cu;
  double pref;
  
  rmrc0=rij-rc0vlr;
  rmrc0sq=rmrc0*rmrc0;
  rmrc0cu=rmrc0*rmrc0sq;
  
  pref=A0lr+B0lr*rij;
  vlr=pref*rmrc0cu;
  dvlr=B0lr*rmrc0cu+3.*pref*rmrc0sq;
  ddvlr=6.*B0lr*rmrc0sq+6.*pref*rmrc0;
}

//-------------------------------------------------------------------------//

void VLR::poly1(double rij, double& vlr, double& dvlr, double& ddvlr)
{
  double r0mr,r0mrsq,r0mrcu,r0mrqu,r0mrfi;
  
  r0mr=r0vlr-rij;
  r0mrsq=r0mr*r0mr;
  r0mrcu=r0mr*r0mrsq;
  r0mrqu=r0mr*r0mrcu;
  r0mrfi=r0mr*r0mrqu;
  
  vlr=A1lr+0.5*B1lr*r0mrsq+C1lr*r0mrcu+D1lr*r0mrqu+E1lr*r0mrfi;
  dvlr=-B1lr*r0mr-3.*C1lr*r0mrsq-4.*D1lr*r0mrcu- 5.*E1lr*r0mrqu;
  ddvlr=B1lr+6.*C1lr*r0mr+12.*D1lr*r0mrsq+20.*E1lr*r0mrcu;
}

//-------------------------------------------------------------------------//
#ifdef PowRepl
  template <typename T>
  constexpr auto pow6(T const x) noexcept -> T 
  {
    auto const y=x*x*x;
    return y*y;
  }
#endif

void VLR::vlj(double rij, double& vlr, double& dvlr, double& ddvlr)
{
  double rijsq;
  double sigbyr,sigbyrsx,sigbyrsxsq;
  
  double rijinv=1.0/rij;
//  rijsq=rij*rij;
  
  sigbyr=rijinv*sigLJ;
//  sigbyr=sigLJ/rij;
#ifdef PowRepl
  sigbyrsx=pow6(sigbyr);
#else
  sigbyrsx=pow(sigbyr,6);
#endif
  sigbyrsxsq=sigbyrsx*sigbyrsx;
  
  vlr=fepsLJ*(sigbyrsxsq-sigbyrsx);
  dvlr=-6.*rijinv*fepsLJ*(2.*sigbyrsxsq-sigbyrsx);
  ddvlr=6.*rijinv*rijinv*fepsLJ*(26.*sigbyrsxsq-7.*sigbyrsx);
//  dvlr=-6.*fepsLJ/rij*(2.*sigbyrsxsq-sigbyrsx);
//  ddvlr=6.*fepsLJ/rijsq*(26.*sigbyrsxsq-7.*sigbyrsx);
}

//-------------------------------------------------------------------------//

void VLR::poly2(double rij, double& vlr, double& dvlr, double& ddvlr)
{
  double q,qsq,omq,omqsq,prefq;
  double prefvlr,dprefvlr,ddprefvlr;
  
  prefq=dr12vlri;
  q=prefq*(rij-r1vlr);
  qsq=q*q;
  omq=1.-q;
  omqsq=omq*omq;
  
  prefvlr=A2lr+B2lr*q+C2lr*qsq;
  dprefvlr=B2lr+2.*C2lr*q;
  ddprefvlr=2.*C2lr;
  
  vlr=prefvlr*omqsq;
  dvlr=prefq*(dprefvlr*omqsq-2.*prefvlr*omq);
  ddvlr=prefq*prefq*(ddprefvlr*omqsq-4.*dprefvlr*omq+2.*prefvlr);
}

//=========================================================================//

struct VMR
{
  double Amr1,Amr2,Amr3;
  double Bmr1,Bmr2,Bmr3;
  int imrnbm;
  inline void gamma(size_t,double,double,double,int&,VmrLocDat&);
};

//=======================================================================//

void VMR::gamma(size_t nnbij,double nij,double sumk,double ndb,
int& imrnb,VmrLocDat& vmrloc)
{
  sqdsqup(ndb,0.0,1.0,1.0,vmrloc.sndb,vmrloc.dsndbdndb);
  
  if(nnbij==0) // ip is a free atom 
  {
    imrnb=imrnbm;
    vmrloc.gij=1.;
    vmrloc.dgijdndb=0.;
    vmrloc.dgijdsumk=0.;
    vmrloc.dgijdnij=0.;
  }  
  else      // general case
  {
//C calculate gammaij and derivatives
    double gnum=1.+(Amr1+Amr2*vmrloc.sndb)*nij+Amr3*sumk;
    double gden=1.+(Bmr1+Bmr2*vmrloc.sndb)*(nij+Bmr3*sumk);
    vmrloc.gij=gnum/gden;
    if(vmrloc.gij>0.0)
    {
      if(vmrloc.gij<1.0)
      {
        imrnb=1;
        double pref=1./gden;
        vmrloc.dgijdndb=pref*(Amr2*nij-vmrloc.gij*Bmr2*(nij+Bmr3*sumk))*vmrloc.dsndbdndb;
        vmrloc.dgijdsumk=pref*(Amr3-vmrloc.gij*Bmr3*(Bmr1+Bmr2*vmrloc.sndb));
        vmrloc.dgijdnij=pref*(Amr1+Amr2*vmrloc.sndb-vmrloc.gij*(Bmr1+Bmr2*vmrloc.sndb));
      }  
      else
      {
        imrnb=imrnbm;
        vmrloc.dgijdnij=0.;
      }
    }  
    else{imrnb=0;}
  }
}

//=======================================================================//
  
inline void sqdsqdown(double q,double q1,double q2,double dqinv,double& Sq,double& dSq)
{
  double xq,onemxq;
  
  if(q<=q1)
  {
    Sq=1.;
    dSq=0.;
  }
  else if(q>=q2)
  {
    Sq=0.;
    dSq=0.;
  }
  else 
  {
    xq=(q-q1)*dqinv;
    onemxq=(1.-xq);
    Sq=(1.+2.*xq)*onemxq*onemxq;
    dSq=-6.*dqinv*xq*onemxq;
  }
}

//=======================================================================//

inline void sqdsqup(double q,double q1,double q2,double dqinv,double& Sq,double& dSq)
{
  if(q<=q1){Sq=0.;dSq=0.;}
  else if(q>=q2){Sq=1.;dSq=0.;}
  else
  {
    double xq,onemxq;
    xq=(q2-q)*dqinv;
    onemxq=(1.-xq);
    Sq=(1.+2.*xq)*onemxq*onemxq;
    dSq=6.*dqinv*xq*onemxq;
  }
}

//=======================================================================//
