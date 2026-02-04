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

//#include "main.h"
//#include "structures.h"
#include "lchbop_utils.h"

//inline void sclprod(double v1[3], double v2[3], double& v1v2)
//{
//v1v2=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
//}

//----------------------------------------------------------------------//

//inline void vecprod(double v1[3], double v2[3], double (&v3)[3])
//{
//  v3[0]=v1[1]*v2[2]-v1[2]*v2[1];
//  v3[1]=v1[2]*v2[0]-v1[0]*v2[2];
//  v3[2]=v1[0]*v2[1]-v1[1]*v2[0];
//}

//----------------------------------------------------------------------//

void fcpdown(double r,double r1,double r2,double drinv,double p2,
double& fc,double& dfc)
{
  double xr=(r-r1)*drinv;
  double xrsq=xr*xr;
  double onemxr=(1.-xr);
  fc=(1.+2.*xr+p2*xrsq)*onemxr*onemxr;
  dfc=((2.+2.*p2*xr)*onemxr*onemxr-2.*(1.+2.*xr+p2*xrsq)*onemxr)*drinv;
  if(fc<1e-15)fc=1e-15;
  else if(fc>0.999999999999999e+0)fc=0.999999999999999e+0;
}

//----------------------------------------------------------------------//
  
void sqdsqdown0(double q,double q1,double q2,double dqinv,
double& Sq,double& dSq)
{
  double xq=(q-q1)*dqinv;
  double onemxq=(1.-xq);
  Sq=(1.+2.*xq)*onemxq*onemxq;
  dSq=-6.*dqinv*xq*onemxq;
}

////----------------------------------------------------------------------//
//  
//void SqdSqpdown(double q,double q1,double q2,double dqinv,double p2,double& Sq,double& dSq)
//{
//  double xq,onemxq,twoxq;
//  
//  if(q<=q1)
//  {
//    Sq=1.;
//    dSq=0.;
//  }
//  else if(q>=q2)
//  {
//    Sq=0.;
//    dSq=0.;
//  }
//  else
//  {
//    xq=(q-q1)*dqinv;
//    onemxq=(1.-xq);
//    twoxq=2.*xq;
//    Sq=(1.+2.*xq+p2*xq*xq)*onemxq*onemxq;
//    dSq=-dqinv*twoxq*(p2*twoxq-p2+3.)*onemxq;
//  }
//}
//
////----------------------------------------------------------------------//
//  
//void SqdSqpdown_alt(double q,double q1,double q2,double dqinv,double p2,double*Sq,double *dSq)
//{
//  double xq,onemxq,twoxq;
//  
//  if(q<=q1)
//  {
//    *Sq=1.;
//    *dSq=0.;
//  }
//  else if(q>=q2)
//  {
//    *Sq=0.;
//    *dSq=0.;
//  }
//  else
//  {
//    xq=(q-q1)*dqinv;
//    onemxq=(1.-xq);
//    *Sq=(1.+3.*xq+6.*xq*xq)*onemxq*onemxq*onemxq;
//    *dSq=-30.*xq*xq*onemxq*onemxq*dqinv;
//  }
//}
//
//----------------------------------------------------------------------//

void SqdSqup(double q,double q1,double q2,double dqinv,double& Sq,double& dSq)
{
  if(q<=q1)
  {
    Sq=0.;
    dSq=0.;
  }
  else if(q>=q2)
  {
    Sq=1.;
    dSq=0.;
  }
  else
  {
    double xq=(q2-q)*dqinv;
    double onemxq=(1.-xq);
    Sq=(1.+2.*xq)*onemxq*onemxq;
    dSq=6.0*dqinv*xq*onemxq;
  }

}

//----------------------------------------------------------------------//

void SqdSqpup(double q,double q1,double q2,double dqinv,double p2,double& Sq,double& dSq)
{
  if(q<=q1)
  {
    Sq=0.;
    dSq=0.;
  }
  else if(q>=q2)
  {
    Sq=1.;
    dSq=0.;
  }
  else
  {
    double xq=(q2-q)*dqinv;
    double onemxq=(1.-xq);
    double twoxq=2.*xq;
    Sq=(1.+2.*xq+p2*xq*xq)*onemxq*onemxq;
    dSq=dqinv*twoxq*(p2*twoxq-p2+3.)*onemxq;
  }

}

////----------------------------------------------------------------------//
//
//void SqdSqpup_alt(double q,double q1,double q2,double dqinv,double p2,double *Sq,double *dSq)
//{
//  double xq,onemxq;
//  
//  if(q<=q1)
//  {
//    *Sq=0.;
//    *dSq=0.;
//  }
//  else if(q>=q2)
//  {
//    *Sq=1.;
//    *dSq=0.;
//  }
//  else
//  {
//    xq=(q2-q)*dqinv;
//    onemxq=(1.-xq);
//    *Sq=(1.+3.*xq+6.*xq*xq)*onemxq*onemxq*onemxq;
//    *dSq=30.*xq*xq*onemxq*onemxq*dqinv;
//  }
//}
//
////----------------------------------------------------------------------//
//
//void fcup(double r,double r1,double r2,double drinv,double *fc)
//{
//  double xr,onemxr;
//  
//  if(r<=r1)
//  {
//    *fc=0.;
//  }
//  else if(r>=r2)
//  {
//    *fc=1.;
//  }
//  else
//  {
//    xr=(r2-r)*drinv;
//    onemxr=(1.-xr);
//    *fc=(1.+2.*xr)*onemxr*onemxr;
//  }
//}
//
////----------------------------------------------------------------------//
//
//void fcdown(double r,double r1,double r2,double drinv,double *fc)
//{
//  double xr,onemxr;
//  
//  if(r<=r1)
//  {
//    *fc=1.;
//  }
//  else if(r>=r2)
//  {
//    *fc=0.;
//  }
//  else
//  {
//    xr=(r-r1)*drinv;
//    onemxr=(1.-xr);
//    *fc=(1.+2.*xr)*onemxr*onemxr;
//  }
//}
//
////----------------------------------------------------------------------//
//
//void fcpup(double r,double r1,double r2,double drinv,double p2,double *fc)
//{
//  double xr,xrsq,onemxr;
//  
//  if(r<=r1)
//  {
//    *fc=0.;
//  }
//  else if(r>=r2)
//  {
//    *fc=1.;
//  }
//  else
//  {
//    xr=(r2-r)*drinv;
//    xrsq=xr*xr;
//    onemxr=(1.-xr);
//    *fc=(1.+2.*xr+p2*xrsq)*onemxr*onemxr;
//  }
//}
//
////----------------------------------------------------------------------//
//
//void fcpdown(double r,double r1,double r2,double drinv,double p2,double *fc,double *dfc)
//{
//  double xr,xrsq,onemxr;
//  
//  if(r<=r1)
//  {
//    *fc=1.;
//    *dfc=0;
//  }
//  else if(r>=r2)
//  {
//    *fc=0.;
//    *dfc=0;
//  }
//  else
//  {
//    xr=(r-r1)*drinv;
//    xrsq=xr*xr;
//    onemxr=(1.-xr);
//    *fc=(1.+2.*xr+p2*xrsq)*onemxr*onemxr;
//    *dfc=((2.+2.*p2*xr)*onemxr*onemxr-2.*(1.+2.*xr+p2*xrsq)*onemxr)*drinv;
//  }
//
//}
//
////----------------------------------------------------------------------//

void fcpdownrc(double r,double r1,double r2,double drinv,double p2,
double& fc,double& dfcdr,double& dfcdrc)
{
  if(r>=r2)
  {
    fc=0.;
    dfcdr=0;
    dfcdrc=0;
  }
  else
  {
    double xr=(r-r1)*drinv;
    double xrsq=xr*xr;
    double onemxr=(1.-xr);
    double onemxrsq=onemxr*onemxr;
    double dfc=(2.+2.*p2*xr)*onemxrsq-2.*(1.+2.*xr+p2*xrsq)*onemxr;
    fc=(1.+2.*xr+p2*xrsq)*onemxrsq;
//    if (p2<4.)
//    {
//      double dfc=(2.+2.*p2*xr)*onemxrsq-2.*(1.+2.*xr+p2*xrsq)*onemxr;
//      fc=(1.+2.*xr+p2*xrsq)*onemxrsq;
//    }
//    else
//    {
//      double dfc=-30.*xrsq*onemxrsq;
//      fc=(1.+3.*xr+6.*xrsq)*onemxrsq*onemxr;
//    }
    dfcdr=dfc*drinv;
    dfcdrc=-dfc*xr*drinv;
  }
}

////----------------------------------------------------------------------//
//
//void fcvdown(double r,double r1,double r2,double drinv,double *fc)
//{
//  double p2=3.,p3=3.;
//  double xr,xrsq,xrcu,onemxr,onemxrsq;
//  
//  if(r<=r1)
//  {
//    *fc=1.;
//  }
//  else if(r>=r2)
//  {
//    *fc=0.;
//  }
//  else
//  {
//    xr=(r-r1)*drinv;
//    xrsq=xr*xr;
//    xrcu=xr*xrsq;
//    onemxr=(1.-xr);
//    onemxrsq=onemxr*onemxr;
//    *fc=(1.+2.*xr+p2*xrsq+p3*xrcu)*onemxrsq;
//  }
//}
//
//======================================================================//
//
//  template<typename T>
//  T typeconv(T x){return x;}  
//  inline int typeconv(int8_t x){return static_cast<int>(x);}  
//  inline int typeconv(uint8_t x){return static_cast<int>(x);}  
//
////======================================================================//
//     
//  template<typename T,size_t N>
//  void wrtptbl(size_t n,int iw1,int iw2,T (&table)[N])
//  {
//    wrtptbl(n,iw1,iw2,vector<std::remove_const_t<T>>(table,table+N));
//  }
//
////======================================================================//
//     
//  template<typename T>
//  void wrtptbl(size_t n,int iw1,int iw2,const vector<T>& table)
//  {
//    std::ofstream TBLDAT("TMPDAT/tbl.dat",std::ios_base::out);
////    std::ofstream TBLDAT("TMPDAT/tbl.dat",std::ios_base::app);
////    TBLDAT << std::setprecision(6);
////    TBLDAT << std::setprecision(15);
//    assert(table.size()>=n);
//    for(size_t i=0;i<n;i++){TBLDAT << std::setw(iw1) << std::setprecision(iw1) << i << "    " 
//                                   << std::setw(iw2) << std::setprecision(iw2) << typeconv(table[i]) << '\n';}  
////    for(size_t i=0;i<n;i++){TBLDAT << i << '\t' << typeconv(table[i]) << '\n';}  
//  } // end of routine wrtptbl
//
//======================================================================//

