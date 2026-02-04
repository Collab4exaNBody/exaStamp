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

#include <cstdint>
#include <vector>
#include <cstdlib>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_def.h>

using std::vector;
using std::pair;
using std::cin;
using std::cout;
using std::endl;
using std::sqrt;

#include "lchbop_common.h"



namespace exaStamp
{


//=======================================================================//

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C Suggested new implementation
//struct VCTbl
//{
//  vector<<pair<size_t,int>> pmap;
//  vector<<pair<size_t,size_t>> pmapinv;
//  vector<Vec3d> atpos;
//};
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

//=======================================================================//

struct VCStr
{
//  static constexpr uint8_t npcm=10;
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C    FOR TESTING ONLY
  int rank,iprint=0,icnt=0; 
//  int ipmap[41984],ipmapinv[41984];
//  size_t iclt=22,ipct=51,jclt=25,jpct=7,kclt=22,kpct=51,lclt=22,lpct=51;
//  size_t ivclt,jvclt,kvclt,lvclt,ipt=14243,jpt=2824,kpt,lpt;
//  int ipvct,jpvct,kpvct,lpvct;
//  vector<int> clown;    // FOR TESTING ONLY
//  clock_t time1,time2;  
//  double time1,time2;  
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  bool NPT,mkvcs,newgrid;
  int npvcm;
  size_t ncl,nvcl,nvclcb,mvclxy,mvclcbxy,np,npcb;
  size_t nnbm,nnbmrcnm,npcbp2,nglst1m,nglst2m;
  ssize_t ngls0;
  IJK dims,dimscb,ijkcb,mvcl,mvclcb,mvclg;
  Mat3d hmat;
  double decl,rcutinc0;
//  Mat3d hmat={1.,0.,0.,0.,1.,0.,0.,0.,1.};
  int nnbc[4],ijcl[172];
  Vec3d xyz0,dvcl,dvcli;
  vector<int8_t> incb;
  vector<size_t> ip0vcl;
  vector<int8_t> istatcl;
  vector<size_t> ivclpcl;
  vector<size_t> ivclp;
  vector<uint8_t> ispc;
  vector<int8_t> istat;
//  vector<size_t> ip0vcl,ivclp;
  vector<pair<size_t,int>>icllst;

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C Present implementation
  vector<vector<pair<size_t,int>>> pmap;
  vector<vector<pair<size_t,size_t>>> pmapinv;
  vector<vector<Vec3d>> atpos;
//  vector<vector<Vec3d>> atposold;
//C Suggested new implementation
//  vector<VCTbl> vct;
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  vector<Vec3d> atposcl;

  inline void initvcs1(size_t);
  inline void initvcs2(int*,IJK*);
//  int iprcb,nnbc[4];
//  IJK nbclst[172];
//  inline void initvcs2();
//  inline void initvcs3(size_t);
};

//---------------------------------------------------------------------//

inline void VCStr::initvcs1(size_t npcm)
{
  npcb=0;
  npvcm=0;
  istatcl.resize(ncl);
  ispc.resize(np);
  istat.resize(np);
  ivclp.resize(np);
  ivclpcl.resize(npcm);
  atposcl.resize(npcm);

//  clown.resize(ncl);    // FOR TESTING ONLY

  pmap.resize(ncl);
  for(size_t icl=0;icl<ncl;icl++)pmap[icl].clear();
}

//---------------------------------------------------------------------//

inline void VCStr::initvcs2(int* nnbc0,IJK* nbclst)
//inline void VCStr::initvcs2()
{
//  printf("nvcl,nvclcb :%6zu%6zu\n",nvcl,nvclcb);
  incb.resize(nvcl);
  ip0vcl.resize(nvcl);
  icllst.resize(nvclcb);

  std::fill(incb.begin(),incb.end(),0);

  pmapinv.resize(nvcl);
  atpos.resize(nvcl);
  for(size_t ivcl=0;ivcl<nvcl;ivcl++)
  {
    pmapinv[ivcl].clear();
    atpos[ivcl].clear();
  }  

  mvclxy=mvcl.i*mvcl.j;
  mvclcbxy=mvclcb.i*mvclcb.j;

  for(ssize_t iclk=mvclg.k;iclk<mvcl.k-mvclg.k;iclk++)
  for(ssize_t iclj=mvclg.j;iclj<mvcl.j-mvclg.j;iclj++)
  for(ssize_t icli=mvclg.i;icli<mvcl.i-mvclg.i;icli++)
  {
    IJK iclijk{icli,iclj,iclk};
    size_t icl=grid_ijk_to_index(mvcl,iclijk);
    if(icli==mvclg.i||icli==mvcl.i-mvclg.i-1||
       iclj==mvclg.j||iclj==mvcl.j-mvclg.j-1||
       iclk==mvclg.k||iclk==mvcl.k-mvclg.k-1){incb[icl]=1;}
//       iclk==mvclg.k||iclk==mvcl.k-mvclg.k-1){incb.at(icl)=1;}
    else if(icli==mvclg.i+1||icli==mvcl.i-mvclg.i-2||
            iclj==mvclg.j+1||iclj==mvcl.j-mvclg.j-2||
            iclk==mvclg.k+1||iclk==mvcl.k-mvclg.k-2){incb[icl]=2;}
    else if(icli==mvclg.i+2||icli==mvcl.i-mvclg.i-3||
            iclj==mvclg.j+2||iclj==mvcl.j-mvclg.j-3||
            iclk==mvclg.k+2||iclk==mvcl.k-mvclg.k-3){incb[icl]=3;}
    else{incb[icl]=4;}
  }
    
//C Building table vcs.ijcl
  ssize_t icli=mvclg.i;
  ssize_t iclj=mvclg.j;
  ssize_t iclk=mvclg.k;
  int inbc=0;
  IJK iclijk{icli,iclj,iclk};
  size_t icl=grid_ijk_to_index(mvcl,iclijk);
  for(int inbl=0;inbl<4;inbl++) 
  {
    nnbc[inbl]=nnbc0[inbl];
    for(;inbc<nnbc[inbl];inbc++) 
    {
//      printf("inbc,nbclst :%6d,%6zd%6zd%6zd\n",i,nbclst[i].i,nbclst[i].j,nbclst[i].k);
//      cin.get();
      IJK jclijk=iclijk+nbclst[inbc];
      size_t jcl=grid_ijk_to_index(mvcl,jclijk);
      ijcl[inbc]=jcl-icl;
    }  
  }  
}

//---------------------------------------------------------------------//
//
//inline void VCStr::initvcs3(size_t npcm)
//{
//  ivclpcl.resize(npcm);
//  atposcl.resize(npcm);
//}  
//
//  if(icli==2&&iclj==3&&iclk==5)
//  {
//    printf("iclijk,icl :%8zd%8zd%8zd%8zu\n",icli,iclj,iclk,icl);
//  }
//  IJK iclijk{2,3,5};
//  size_t ivcl=grid_ijk_to_index(mvcl,iclijk);
//  printf("iclijk    :%8zd%8zd%8zd\n",iclijk.i,iclijk.j,iclijk.k);
//  printf("ivcl,incb :%6zd%6d\n",ivcl,incb[ivcl]);
//  cin.get();
//
//=======================================================================//

  struct NBStr
  {
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//C For testing only
    size_t npt,iclt,ipct,jclt,jpct,kclt,kpct,lclt,lpct;
    size_t ivclt,jvclt,kvclt,lvclt,ipt,jpt,kpt,lpt;
    int ipvct,jpvct,kpvct,lpvct;
//    int ipmap[41984],ipmapinv[41984];
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    size_t ncl,nnbsrm,inbm,inbmrm,nnnbm,nbndmrm,nbndm,nbndmrcnm;
// size = number of particles
    std::vector<size_t> nnb0,nnbfrac,nnbfull,nnbmrcn0,nnbmrcn;
    std::vector<size_t> iclp;
    std::vector<uint8_t> ispc;
// size = half time total number of neighbors of particles in central box
    std::vector<size_t> srfl,srfr,mrcn;
    std::vector<Vec3d> vijsrfl,vijsrfr,vijmrcn;
    std::vector<double> rijsrfl,rijsrfr,rijmrcn;
// size = total number of neighbors of all relevant particles
    std::vector<size_t> nbtbl,indji,isr2mr; 
    std::vector<double> rstor,crd,ndb;
    std::vector<Vec3d> sigstor;
    std::vector<vector<double>> snstor,crdspc;
    std::vector<double> esrpp,elrpp;
    std::vector<Vec3d> force;

    size_t nmr;
    std::vector<size_t>nnbmr,indmr,indmrji,nextind,nbtblmr,nderndb,ndersumk;
    std::vector<size_t>iderndb0,nderatndb,idersumk0,nderatsumk;
    std::vector<double>rstormr;
    std::vector<double>gij,smr,rcmr,dsmrdgij,drcmrdndb,dgijdnij,dgijdndb,dgijdsumk;
    std::vector<double>snsmr,dsnsmrdr,dsnsmrdndb,dsnsmrdnik,dsnsmrdsumk;
    std::vector<Vec3d>frcndb,frcsumk,sigstormr;
//#ifdef VIRIEL
    std::vector<Mat3d>viriel,vrlndb,vrlsumk;
//#endif    

    void initnbs1(size_t,size_t,size_t,size_t);
    void initnbs2(size_t);
  };

//  using NBStr=std::vector<NbtblCell>;

//----------------------------------------------------------------------//

  inline void NBStr::initnbs1(size_t np,size_t npcbp2,size_t nnbm,size_t nnbmrcnm)
  {
    nnb0.resize(np);
    nnbfull.resize(np);
    nnbfrac.resize(np);
    nnbmrcn0.resize(np);
    nnbmrcn.resize(np);
    crd.resize(np);
    crdspc.resize(np,vector<double>(2));
    ndb.resize(np);
    esrpp.resize(np);
    elrpp.resize(np);
    force.resize(np);
#ifdef VIRIEL
#ifdef VIRIELLOCAL
    viriel.resize(np);
#else
    viriel.resize(1);
#endif
#endif

    nbndm=3*npcbp2*nnbm/5;
    nbndmrcnm=npcbp2*nnbmrcnm;
    srfl.resize(nbndm);
    srfr.resize(nbndm);
    rijsrfl.resize(nbndm);
    rijsrfr.resize(nbndm);
    vijsrfl.resize(nbndm);
    vijsrfr.resize(nbndm);

    nnnbm=npcbp2*nnbm;
    nbtbl.resize(nnnbm);
    indji.resize(nnnbm);
    isr2mr.resize(nnnbm);
    rstor.resize(nnnbm);
    sigstor.resize(nnnbm);
    snstor.resize(nnnbm,vector<double>(3));

    nnbmr.resize(np);
    indmr.resize(np);

    nbndmrcnm=npcbp2*nnbmrcnm;
    mrcn.resize(nbndmrcnm);
    rijmrcn.resize(nbndmrcnm);
    vijmrcn.resize(nbndmrcnm);

    nbndmrm=npcbp2*nnbm;
    indmrji.resize(nbndmrm);
    nextind.resize(nbndmrm);
    nbtblmr.resize(nbndmrm);
    rstormr.resize(nbndmrm);
    sigstormr.resize(nbndmrm);
    snsmr.resize(nbndmrm);
    dsnsmrdr.resize(nbndmrm);
    dsnsmrdnik.resize(nbndmrm);
    dsnsmrdsumk.resize(nbndmrm);
    dsnsmrdndb.resize(nbndmrm);
    gij.resize(nbndmrm);
    smr.resize(nbndmrm);
    rcmr.resize(nbndmrm);
    dsmrdgij.resize(nbndmrm);
    drcmrdndb.resize(nbndmrm);
    dgijdnij.resize(nbndmrm);
    dgijdndb.resize(nbndmrm);
    dgijdsumk.resize(nbndmrm);
    nderndb.resize(nbndmrm);
    ndersumk.resize(nbndmrm);
    nderatndb.resize(10*nbndmrm);
    nderatsumk.resize(10*nbndmrm);
    frcndb.resize(10*nbndmrm);
    frcsumk.resize(10*nbndmrm);
#ifdef VIRIEL
#ifdef VIRIELLOCAL
    vrlndb.resize(10*nbndmrm);
    vrlsumk.resize(10*nbndmrm);
#else
    vrlndb.resize(nbndmrm);
    vrlsumk.resize(nbndmrm);
#endif
#endif
//    frcndb.resize(10*nbndmrm,{0.,0.,0.});  // ??????????????? //
//    frcsumk.resize(10*nbndmrm,{0.,0.,0.});
  }

//----------------------------------------------------------------------//

  inline void NBStr::initnbs2(size_t np)
  {
    nmr=0;
    nderndb[0]=1;
    ndersumk[0]=1;
    std::vector<double> crdspc0(2,0.0);
    Vec3d force0={0.0,0.0,0.0};
    Mat3d viriel0={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
// Setting for maximal (re-)initialisation
//    std::fill(nnb0.begin(),nnb0.end(),0);
//    std::fill(nnbfull.begin(),nnbfull.end(),0);
//    std::fill(nnbfrac.begin(),nnbfrac.end(),0);
//    std::fill(nnbmrcn.begin(),nnbmrcn.end(),0);
//    std::fill(nnbmrcn0.begin(),nnbmrcn0.end(),0);
//    std::fill(nnbmr.begin(),nnbmr.end(),0);
//    std::fill(indmr.begin(),indmr.end(),0);
//    std::fill(nextind.begin(),nextind.end(),0);
//    std::fill(ndb.begin(),ndb.end(),0.0);
//    std::fill(crd.begin(),crd.end(),0.0);
//    std::fill(crdspc.begin(),crdspc.end(),crdspc0);
//    std::fill(esrpp.begin(),esrpp.end(),0.0);
//    std::fill(elrpp.begin(),elrpp.end(),0.0);
//    std::fill(force.begin(),force.end(),force0);
//    std::fill(viriel.begin(),viriel.end(),viriel0);
//
//    std::fill(nbtbl.begin(),nbtbl.end(),0);
//    std::fill(indji.begin(),indji.end(),0);
//    std::fill(isr2mr.begin(),isr2mr.end(),0);
//
//    std::fill(mrcn.begin(),mrcn.end(),0);
//    std::fill(indmrji.begin(),indmrji.end(),0);
//
//    std::fill(snsmr.begin(),snsmr.end(),0.0);
//    std::fill(dsnsmrdr.begin(),dsnsmrdr.end(),0.0);
//    std::fill(dsnsmrdnik.begin(),dsnsmrdnik.end(),0.0);      // <---- this one
//    std::fill(dsnsmrdsumk.begin(),dsnsmrdsumk.end(),0.0);
//    std::fill(dsnsmrdndb.begin(),dsnsmrdndb.end(),0.0);
//    std::fill(gij.begin(),gij.end(),0.0);
//    std::fill(smr.begin(),smr.end(),0.0);
//    std::fill(rcmr.begin(),rcmr.end(),0.0);
//    std::fill(dsmrdgij.begin(),dsmrdgij.end(),0.0);          // <---- this one
//    std::fill(drcmrdndb.begin(),drcmrdndb.end(),0.0);
//    std::fill(dgijdnij.begin(),dgijdnij.end(),0.0);
//    std::fill(dgijdndb.begin(),dgijdndb.end(),0.0);
//    std::fill(dgijdsumk.begin(),dgijdsumk.end(),0.0);
//
//    std::fill(frcndb.begin(),frcndb.end(),force0);
//    std::fill(frcsumk.begin(),frcsumk.end(),force0);

// Additional initializations
//#ifdef VIRIEL
//    std::fill(vrlndb.begin(),vrlndb.end(),viriel0);
//    std::fill(vrlsumk.begin(),vrlsumk.end(),viriel0);
//#endif

// Setting for minimal (re-)initialisation
//-----------------------------------------------------------------------
    std::fill(nnb0.begin(),nnb0.end(),0);
    std::fill(nnbfull.begin(),nnbfull.end(),0);
    std::fill(nnbfrac.begin(),nnbfrac.end(),0);
//    std::fill(nnbmrcn.begin(),nnbmrcn.end(),0);
//    std::fill(nnbmrcn0.begin(),nnbmrcn0.end(),0);
    std::fill(nnbmr.begin(),nnbmr.end(),0);
//    std::fill(indmr.begin(),indmr.end(),0);
//    std::fill(nextind.begin(),nextind.end(),0);
//    std::fill(ndb.begin(),ndb.end(),0.0);
    std::fill(crd.begin(),crd.end(),0.0);
    std::fill(crdspc.begin(),crdspc.end(),crdspc0);
    std::fill(esrpp.begin(),esrpp.end(),0.0);
    std::fill(elrpp.begin(),elrpp.end(),0.0);
    std::fill(force.begin(),force.end(),force0);
#ifdef VIRIEL
    std::fill(viriel.begin(),viriel.end(),viriel0);
#endif

//    std::fill(nbtbl.begin(),nbtbl.end(),0);
//    std::fill(indji.begin(),indji.end(),0);
//    std::fill(isr2mr.begin(),isr2mr.end(),0);
  
//    std::fill(mrcn.begin(),mrcn.end(),0);
    std::fill(indmrji.begin(),indmrji.end(),0);

//    std::fill(snsmr.begin(),snsmr.end(),0.0);
//    std::fill(dsnsmrdr.begin(),dsnsmrdr.end(),0.0);
//    std::fill(dsnsmrdnik.begin(),dsnsmrdnik.end(),0.0);      // <---- this one
//    std::fill(dsnsmrdsumk.begin(),dsnsmrdsumk.end(),0.0);
//    std::fill(dsnsmrdndb.begin(),dsnsmrdndb.end(),0.0);
//    std::fill(gij.begin(),gij.end(),0.0);
//    std::fill(smr.begin(),smr.end(),0.0);
//    std::fill(rcmr.begin(),rcmr.end(),0.0);
//    std::fill(dsmrdgij.begin(),dsmrdgij.end(),0.0);          // <---- this one
//    std::fill(drcmrdndb.begin(),drcmrdndb.end(),0.0);
//    std::fill(dgijdnij.begin(),dgijdnij.end(),0.0);
//    std::fill(dgijdndb.begin(),dgijdndb.end(),0.0);
//    std::fill(dgijdsumk.begin(),dgijdsumk.end(),0.0);

//    std::fill(frcndb.begin(),frcndb.end(),force0);
//    std::fill(frcsumk.begin(),frcsumk.end(),force0);

// Additional initializations
#ifdef VIRIEL
    std::fill(vrlndb.begin(),vrlndb.end(),viriel0);
    std::fill(vrlsumk.begin(),vrlsumk.end(),viriel0);
#endif
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

//    for(size_t ip=0;ip<np;ip++)nnb0[ip]=0;
//    for(size_t ip=0;ip<np;ip++)nnbfull[ip]=0;
//    for(size_t ip=0;ip<np;ip++)nnbfrac[ip]=0;
//    for(size_t ip=0;ip<np;ip++)crd[ip]=0.0;
//    for(size_t ip=0;ip<np;ip++){crdspc[0][ip]=0.0;crdspc[1][ip]=0.0;}
//    for(size_t ip=0;ip<np;ip++)esrpp[0][ip]=0.0;
//    for(size_t ip=0;ip<np;ip++)elrpp[0][ip]=0.0;

//    for(size_t ip=0;ip<np;ip++)
//    {
//      nnb0[ip]=0;
//      nnbfull[ip]=0;
//      nnbfrac[ip]=0;
//      crd[ip]=0.0;
//      crdspc[ip][0]=0.0;
//      crdspc[ip][1]=0.0;
//      esrpp[ip]=0.0;
//      elrpp[ip]=0.0;
//      force[ip]={0.,0.,0.};
//      viriel[ip]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
//
//      nnbmr[ip]=0;
//    }  

//    for(size_t i=0;i<nbndmrm;i++)indmrji[i]=0;

//      printf("nnb0           :%8zu%8zu%8zu\n",nnb0[0],nnb0[1],nnb0[np-1]);
//      printf("nnbfull        :%8zu%8zu%8zu\n",nnbfull[0],nnbfull[1],nnbfull[np-1]);
//      printf("nnbfrac        :%8zu%8zu%8zu\n",nnbfrac[0],nnbfrac[1],nnbfrac[np-1]);
//      printf("nnbmr          :%8zu%8zu%8zu\n",nnbmr[0],nnbmr[1],nnbmr[np-1]);
//      printf("indmr          :%8zu%8zu%8zu\n",indmr[0],indmr[1],indmr[np-1]);
//      printf("crd            :%18.10lf%18.10lf%18.10lf\n",crd[0],crd[1],crd[np-1]);
//      printf("crdspc[0]      :%18.10lf%18.10lf%18.10lf\n",crdspc[0][0],crdspc[1][0],crdspc[np-1][0]);
//      printf("crdspc[1]      :%18.10lf%18.10lf%18.10lf\n",crdspc[0][1],crdspc[1][1],crdspc[np-1][1]);
//      printf("esrpp          :%18.10lf%18.10lf%18.10lf\n",esrpp[0],esrpp[1],esrpp[np-1]);
//      printf("elrpp          :%18.10lf%18.10lf%18.10lf\n",elrpp[0],elrpp[1],elrpp[np-1]);
//      printf("force[0]       :%18.10lf%18.10lf%18.10lf\n",force[0].x,force[0].y,force[0].z);
//      printf("force[np-1]    :%18.10lf%18.10lf%18.10lf\n",force[np-1].x,force[np-1].y,force[np-1].z);
//      printf("virielx[0]     :%18.10lf%18.10lf%18.10lf\n",viriel[0].m11,viriel[0].m11,viriel[0].m11);
//      printf("viriely[0]     :%18.10lf%18.10lf%18.10lf\n",viriel[0].m21,viriel[0].m21,viriel[0].m21);
//      printf("virielz[0]     :%18.10lf%18.10lf%18.10lf\n",viriel[0].m31,viriel[0].m31,viriel[0].m31);
//      printf("virielx[np-1]  :%18.10lf%18.10lf%18.10lf\n",viriel[np-1].m11,viriel[np-1].m11,viriel[np-1].m11);
//      printf("viriely[np-1]  :%18.10lf%18.10lf%18.10lf\n",viriel[np-1].m21,viriel[np-1].m21,viriel[np-1].m21);
//      printf("virielz[np-1]  :%18.10lf%18.10lf%18.10lf\n",viriel[np-1].m31,viriel[np-1].m31,viriel[np-1].m31);
//      cin.get();
  }

//======================================================================//
  
  struct LNBStr
  {
    LNBStr(size_t);
    ~LNBStr(){}
    int nnbsrm,nbrm,nnbmloc;     // for testing only 
    int nnbfl,nnbfr,nnb,inbj; 
//    size_t nnbsrm,nbrm,nnbmloc;     // for testing only 
//    size_t nnbfl,nnbfr,nnb,inbj; 
    double crdmr,smr,rcmr,pfeff,rcmreff,drcmreffdndb,drcmreffdgij,dsmreffdgij;
    double dvsmradndb,dvsmradsumk,dvsmradnij;
    vector<uint8_t> ispc,isdbij;
    vector<size_t> nbtbl,isr2mr;
    std::vector<double> rstor,dbijdndb,dbijdsumk;
    std::vector<Vec3d> sigstor,dbij;
    std::vector<vector<double>> snstor;
#ifdef VIRIEL
    std::vector<Mat3d> virielbij;
#endif
  };

//----------------------------------------------------------------------//
  
  inline LNBStr::LNBStr(size_t nnbmloc)
  {
    nbtbl.resize(nnbmloc);
    ispc.resize(nnbmloc);
    isr2mr.resize(nnbmloc,0);
    isdbij.resize(nnbmloc);
    rstor.resize(nnbmloc);
    dbijdndb.resize(nnbmloc);
    dbijdsumk.resize(nnbmloc);
    sigstor.resize(nnbmloc);
    dbij.resize(nnbmloc);
    snstor.resize(nnbmloc,vector<double>(3));
#ifdef VIRIEL
#ifdef VIRIELLOCAL
    virielbij.resize(nnbmloc);
#else
    virielbij.resize(1);
#endif
#endif
//    printf("getEF; viriel :%20.12lf%20.12lf%20.12lf\n",virielbij[0].m11,virielbij[0].m12,virielbij[0].m13);
//    printf("getEF; viriel :%20.12lf%20.12lf%20.12lf\n",virielbij[0].m21,virielbij[0].m22,virielbij[0].m23);
//    printf("getEF; viriel :%20.12lf%20.12lf%20.12lf\n",virielbij[0].m31,virielbij[0].m32,virielbij[0].m33);
//    cin.get();
  }

//======================================================================//
  
  struct FATPLoc
  {
    FATPLoc(size_t);
    ~FATPLoc(){}
    vector<int> iscnb;
    vector<vector<size_t>> lstfrc3;
    vector<double> prwij;
    vector<vector<Vec3d>> frc1,frc2,frc3;
#ifdef VIRIEL
    vector<vector<Mat3d>> virielfrc1,virielfrc2,virielfrc3;
#endif
  };

//----------------------------------------------------------------------//
  
  inline FATPLoc::FATPLoc(size_t nnbsrm)
  {
    size_t nnbsrmsq=nnbsrm*nnbsrm;
    iscnb.resize(nnbsrm);
    lstfrc3.resize(2,vector<size_t>(nnbsrmsq));
    prwij.resize(nnbsrm+1);
    frc1.resize(2,vector<Vec3d>(nnbsrm));
    frc2.resize(2,vector<Vec3d>(nnbsrm));
    frc3.resize(2,vector<Vec3d>(nnbsrmsq));
#ifdef VIRIEL
#ifdef VIRIELLOCAL
    virielfrc1.resize(2,vector<Mat3d>(nnbsrm));
    virielfrc2.resize(2,vector<Mat3d>(nnbsrm));
    virielfrc3.resize(2,vector<Mat3d>(nnbsrmsq));
#else
    virielfrc1.resize(2,vector<Mat3d>(1));
    virielfrc2.resize(2,vector<Mat3d>(1));
    virielfrc3.resize(2,vector<Mat3d>(1));
#endif
#endif
//    for(int ij=0;ij<2;ij++)
//    {
//      for(int inb=0;inb<nnbsrm;inb++){frc1[ij][inb]={0.,0.,0.};}
//      for(int inb=0;inb<nnbsrm;inb++){frc2[ij][inb]={0.,0.,0.};}
//      for(int inb=0;inb<nnbsrmsq;inb++){frc3[ij][inb]={0.,0.,0.};}
//    }  
  }

//======================================================================//
  
  struct VmrLocDat
  {
    size_t nmr0p1,nmr0p2;
    double rcmr[2],drcmrdndb[2];
    double sndb,dsndbdndb,gij,dgijdndb,dgijdsumk,dgijdnij;
  };

//======================================================================//
  
  struct GhostTbls
  {
    GhostTbls(size_t,size_t,size_t);
    ~GhostTbls(){}
    size_t nglst1m,nglst2m;
    size_t nglst1,nglst1sr,nglst2,ngbsrfl1,ngbsrfr1,ngbsrfl2,ngbsrfr2,ngbmrcn; 
    vector<size_t> glst1,glst2;
    vector<pair<size_t,size_t>>gbsrfl1,gbsrfr1,gbsrfl2,gbsrfr2,gbmrcn;
    vector<size_t> gp2cp;   // for testing
  };

//----------------------------------------------------------------------//
  
  inline GhostTbls::GhostTbls(size_t np,size_t ng1m,size_t ng2m)
  {
    nglst1m=ng1m;
    nglst2m=ng2m;
    nglst1=0;
    nglst2=0;
    ngbsrfl1=0;
    ngbsrfr1=0;
    ngbsrfl2=0;
    ngbsrfr2=0;
    ngbmrcn=0;
    glst1.resize(nglst1m);
    glst2.resize(2*nglst2m);
    gp2cp.resize(np);   // for testing
    gbsrfl1.resize(nglst1m);
    gbsrfr1.resize(nglst1m);
    gbsrfl2.resize(nglst2m);
    gbsrfr2.resize(nglst2m);
    gbmrcn.resize(nglst2m);
  }

//======================================================================//

}

inline IJK getmvclcb(VCStr&);

inline IJK getmvclg(VCStr&);

inline IJK getmvcl(VCStr&);

//inline void getdcut(double,Mat3d&,Vec3d&,Vec3d&,Vec3d&);

template<typename CellT>
inline void mapcl2vcl(CellT&,size_t,size_t&,VCStr&);

template<typename CellT>
inline void getvclpcl(CellT&,size_t,size_t,VCStr&);

//======================================================================//
//C Function getvcsgrid builds a new verlet cell structure when this 
//C is required, either after an update of the local grid by exaStamp
//C or, in the case of NPT simulation, when the size the old verlet cell
//C structure does not match the interaction cut-off criteria anymore
//C due to the changing box variables.
     
  static inline void getvcsgrid(VCStr& vcs)
  {
    vcs.mvclcb=getmvclcb(vcs);
    vcs.mvclg=getmvclg(vcs);
    vcs.mvcl=getmvcl(vcs);

    if(vcs.iprint>0)
    {
      printf("rank      :%6d\n",vcs.rank);
      printf("np        :%8zu\n",vcs.np);
      printf("decl      :%14.8lf\n",vcs.decl);
      printf("vcs.dims  :%8zu%8zu%8zu\n",vcs.dims.i,vcs.dims.j,vcs.dims.k);
      printf("mvclcb    :%8zu%8zu%8zu\n",vcs.mvclcb.i,vcs.mvclcb.j,vcs.mvclcb.k);
      printf("mvclg     :%8zu%8zu%8zu\n",vcs.mvclg.i,vcs.mvclg.j,vcs.mvclg.k);
      printf("mvcl      :%8zu%8zu%8zu\n",vcs.mvcl.i,vcs.mvcl.j,vcs.mvcl.k);
      printf("nvclcb    :%8zu\n",vcs.nvclcb);
      printf("nvcl      :%8zu\n",vcs.nvcl);
//      printf("hcb       :%14.8lf%14.8lf%14.8lf\n",hcb.x,hcb.y,hcb.z);
//      printf("rcutsr    :%14.8lf%14.8lf%14.8lf\n",rcutsr.x,rcutsr.y,rcutsr.z);
//      printf("rcutinc   :%14.8lf%14.8lf%14.8lf\n",rcutinc.x,rcutinc.y,rcutinc.z);
//      printf("dghost    :%14.8lf%14.8lf%14.8lf\n",dghost.x,dghost.y,dghost.z);
      printf("dvcl      :%14.8lf%14.8lf%14.8lf\n",vcs.dvcl.x,vcs.dvcl.y,vcs.dvcl.z);
      printf("dvcli     :%14.8lf%14.8lf%14.8lf\n",vcs.dvcli.x,vcs.dvcli.y,vcs.dvcli.z);
      printf("xyz0      :%14.8lf%14.8lf%14.8lf\n",vcs.xyz0.x,vcs.xyz0.y,vcs.xyz0.z);
      cin.get();
    }
  } // end of routine getvcsgrid

//======================================================================//
//C Function chkvcsgrid checks whether there is a change in the
//C verlet cell grid dimensions during NPT simulation mvclcb or mvclg and, 
//C if so, updates the grid dimensions returns true.
     
  static inline bool chkvcsgrid(VCStr& vcs)
  {
    bool vcschg=false;
    IJK mvclcb=getmvclcb(vcs);
    IJK mvclg=getmvclg(vcs);
    if(mvclcb!=vcs.mvclcb){vcschg=true;vcs.mvclcb=mvclcb;}
    if(mvclg!=vcs.mvclg){vcschg=true;vcs.mvclg=mvclg;}
    if(vcschg){vcs.mvcl=getmvcl(vcs);}

    if(vcs.iprint>0)
    {
      printf("chkvcsgrid \n");
//      printf("rank      :%6d\n",vcs.rank);
      printf("np        :%8zu\n",vcs.np);
      printf("decl      :%14.8lf\n",vcs.decl);
      printf("vcs.dims  :%8zu%8zu%8zu\n",vcs.dims.i,vcs.dims.j,vcs.dims.k);
      printf("mvclcb    :%8zu%8zu%8zu\n",vcs.mvclcb.i,vcs.mvclcb.j,vcs.mvclcb.k);
      printf("mvclg     :%8zu%8zu%8zu\n",vcs.mvclg.i,vcs.mvclg.j,vcs.mvclg.k);
      printf("mvcl      :%8zu%8zu%8zu\n",vcs.mvcl.i,vcs.mvcl.j,vcs.mvcl.k);
      printf("nvclcb    :%8zu\n",vcs.nvclcb);
      printf("nvcl      :%8zu\n",vcs.nvcl);
//      printf("hcb       :%14.8lf%14.8lf%14.8lf\n",hcb.x,hcb.y,hcb.z);
//      printf("rcutsr    :%14.8lf%14.8lf%14.8lf\n",rcutsr.x,rcutsr.y,rcutsr.z);
//      printf("rcutinc   :%14.8lf%14.8lf%14.8lf\n",rcutinc.x,rcutinc.y,rcutinc.z);
//      printf("dghost    :%14.8lf%14.8lf%14.8lf\n",dghost.x,dghost.y,dghost.z);
      printf("dvcl      :%14.8lf%14.8lf%14.8lf\n",vcs.dvcl.x,vcs.dvcl.y,vcs.dvcl.z);
      printf("dvcli     :%14.8lf%14.8lf%14.8lf\n",vcs.dvcli.x,vcs.dvcli.y,vcs.dvcli.z);
      printf("xyz0      :%14.8lf%14.8lf%14.8lf\n",vcs.xyz0.x,vcs.xyz0.y,vcs.xyz0.z);
      cin.get();
    }

    return vcschg;
  } // end of fonction chkvcsgrid

//======================================================================//
//C Function getmvclcb gets the dimensions in the verlet cell grid
//C wthin the Verlet Structure Central Box (VSCB).
     
  inline IJK getmvclcb(VCStr& vcs)
  {
    double rcutsr0=2.2;

    Vec3d hcb={vcs.dimscb.i*vcs.decl,vcs.dimscb.j*vcs.decl,vcs.dimscb.k*vcs.decl};
    hcb=hcb+vcs.rcutinc0;
    Vec3d vscb[3]={hcb.x*column1(vcs.hmat),hcb.y*column2(vcs.hmat),hcb.z*column3(vcs.hmat)};

    Vec3d vperp[3];
    vperp[0]=cross(vscb[1],vscb[2]);
    vperp[1]=cross(vscb[2],vscb[0]);
    vperp[2]=cross(vscb[0],vscb[1]);

    double avscb[3];
    avscb[0]=sqrt(norm2(vperp[0]));
    avscb[1]=sqrt(norm2(vperp[1]));
    avscb[2]=sqrt(norm2(vperp[2]));
//    printf("avscb :%14.8lf%14.8lf%14.8lf\n",avscb[0],avscb[1],avscb[2]);

    double volvscb=dot(vperp[0],vscb[0]);
//    printf("volvscb =%14.8lf\n",volvscb);

    Vec3d rcutsr;
    rcutsr.x=(rcutsr0/volvscb)*avscb[0]*hcb.x;
    rcutsr.y=(rcutsr0/volvscb)*avscb[1]*hcb.y;
    rcutsr.z=(rcutsr0/volvscb)*avscb[2]*hcb.z;

    IJK mvclcb{floor(hcb.x/rcutsr.x),floor(hcb.y/rcutsr.y),floor(hcb.z/rcutsr.z)};

    vcs.dvcl={hcb.x/mvclcb.i,hcb.y/mvclcb.j,hcb.z/mvclcb.k};
    vcs.dvcli={1.0/vcs.dvcl.x,1.0/vcs.dvcl.y,1.0/vcs.dvcl.z};
//    printf("dvcl      :%14.8lf%14.8lf%14.8lf\n",vcs.dvcl.x,vcs.dvcl.y,vcs.dvcl.z);
//    printf("dvcli     :%14.8lf%14.8lf%14.8lf\n",vcs.dvcli.x,vcs.dvcli.y,vcs.dvcli.z);

    return mvclcb;
  } // end of routine getmvclcb

//======================================================================//
//C Function getmvclg gets the dimensions in the verlet cell grid
//C within the ghost layer of the verlet cell structure
     
  inline IJK getmvclg(VCStr& vcs)
  {
    double dghost0=8.4;
    double dgh=dghost0+vcs.rcutinc0;
    IJK mvclg{floor(dgh*vcs.dvcli.x),floor(dgh*vcs.dvcli.y),floor(dgh*vcs.dvcli.z)};
    mvclg=mvclg+1;
    vcs.xyz0=vcs.xyz0-0.5*vcs.rcutinc0;
    vcs.xyz0.x=vcs.xyz0.x-mvclg.i*vcs.dvcl.x;
    vcs.xyz0.y=vcs.xyz0.y-mvclg.j*vcs.dvcl.y;
    vcs.xyz0.z=vcs.xyz0.z-mvclg.k*vcs.dvcl.z;
    return mvclg;
  } // end of routine getmvclg

//======================================================================//
//C Function getmvcl gets the dimensions in the verlet cell grid
//C within the total local domain 
     
  inline IJK getmvcl(VCStr& vcs)
  {
    IJK mvcl=vcs.mvclcb+vcs.mvclg*static_cast<ssize_t>(2);
    vcs.nvcl=mvcl.i*mvcl.j*mvcl.k;
    vcs.nvclcb=vcs.mvclcb.i*vcs.mvclcb.j*vcs.mvclcb.k;
    return mvcl;
  } // end of routine getmvcl

//======================================================================//
//     
//  inline void getdcut(double rcutinc0,Mat3d& hmat,Vec3d& rcutsr,Vec3d& rcutinc,Vec3d& dghost)
//  {
//    double dxinv=1.0/sqrt(hmat.m11*hmat.m11+hmat.m21*hmat.m21+hmat.m31*hmat.m31);
//    double dyinv=1.0/sqrt(hmat.m12*hmat.m12+hmat.m22*hmat.m22+hmat.m32*hmat.m32);
//    double dzinv=1.0/sqrt(hmat.m13*hmat.m13+hmat.m23*hmat.m23+hmat.m33*hmat.m33);
//    rcutsr.x=dxinv*2.2;
//    rcutsr.y=dyinv*2.2;
//    rcutsr.z=dzinv*2.2;
//    rcutinc.x=rcutinc0;
//    rcutinc.y=rcutinc0;
//    rcutinc.z=rcutinc0;
//    dghost.x=8.4+rcutinc.x;
//    dghost.y=8.4+rcutinc.y;
//    dghost.z=8.4+rcutinc.z;
////    rcutinc.x=dxinv*rcutinc0;
////    rcutinc.y=dyinv*rcutinc0;
////    rcutinc.z=dzinv*rcutinc0;
////    dghost.x=dxinv*8.4+rcutinc.x;
////    dghost.y=dyinv*8.4+rcutinc.y;
////    dghost.z=dzinv*8.4+rcutinc.z;
//  } // end of function getdcut
//
//======================================================================//
//C Function vclfill fills the verlet cells with particles belonging
//C to these cells and establishes a mapping between the exaStamp
//C cell structure and the verlet cell structure.
     
  template<typename CellT>
  static inline void vclfill(CellT& cells,VCStr& vcs)
  {
//C Mapping content of exaStamp cells inside central box to verlet cells.
    size_t npcnt=0;
    for(ssize_t iclk=vcs.ngls0;iclk<vcs.dims.k-vcs.ngls0;iclk++)
    for(ssize_t iclj=vcs.ngls0;iclj<vcs.dims.j-vcs.ngls0;iclj++)
    for(ssize_t icli=vcs.ngls0;icli<vcs.dims.i-vcs.ngls0;icli++)
    {
      IJK iclijk{icli,iclj,iclk};
      size_t icl=grid_ijk_to_index(vcs.dims,iclijk);
      assert(icl>=0&&icl<vcs.ncl);
      mapcl2vcl(cells,icl,npcnt,vcs);
    }  
//    printf("vclfill 1; npcnt :%8zu\n",npcnt);

//C Initialize icllst for central box
    size_t ilst=0;
    for(ssize_t iclk=vcs.mvclg.k;iclk<vcs.mvcl.k-vcs.mvclg.k;iclk++)
    for(ssize_t iclj=vcs.mvclg.j;iclj<vcs.mvcl.j-vcs.mvclg.j;iclj++)
    for(ssize_t icli=vcs.mvclg.i;icli<vcs.mvcl.i-vcs.mvclg.i;icli++)
    {
      IJK iclijk{icli,iclj,iclk};
      size_t ivcl=grid_ijk_to_index(vcs.mvcl,iclijk);
      int npcbcli=vcs.pmapinv[ivcl].size();
      assert(ilst<vcs.nvclcb);
      vcs.icllst[ilst++]=pair<size_t,int>(ivcl,npcbcli);
    }

//C Mapping content of exaStamp cells outside central box to verlet cells.
    ssize_t igls=vcs.ngls0;
    for(size_t il=0;il<vcs.ngls0;il++)
    {
      igls--;
      ssize_t idcli=0;
      for(ssize_t iclk=igls;iclk<vcs.dims.k-igls;iclk++)
      for(ssize_t iclj=igls;iclj<vcs.dims.j-igls;iclj++)
      {
        if(iclk==igls||iclk==vcs.dims.k-igls-1){idcli=1;}
        else if(iclj==igls||iclj==vcs.dims.j-igls-1){idcli=1;}
        else{idcli=vcs.dims.i-2*igls-1;}
        for(ssize_t icli=igls;icli<vcs.dims.i-igls;icli+=idcli)
        {
          IJK iclijk{icli,iclj,iclk};
          size_t icl=grid_ijk_to_index(vcs.dims,iclijk);
          assert(icl>=0&&icl<vcs.ncl);
          mapcl2vcl(cells,icl,npcnt,vcs);
        }  
      }
    }
//    printf("vclfill 2; npcnt,vcs.np :%8zu%8zu\n",npcnt,vcs.np);
    assert(npcnt==vcs.np);

  } // end of routine vclfill

//======================================================================//
//C Function mapcl2vcl maps the particle content of a given exaStamp
//C cell icl to the verlet cells.
     
  template<typename CellT>
  inline void mapcl2vcl(CellT& cells,size_t icl,size_t& npcnt,VCStr& vcs)
  {
    size_t npcli=cells[icl].size();
    npcnt+=npcli;
    getvclpcl(cells,icl,npcli,vcs);
    vcs.pmap[icl].resize(npcli);
    for(size_t ipc=0;ipc<npcli;ipc++)
    {
      size_t ivcl=vcs.ivclpcl[ipc];
      assert(ivcl<vcs.nvcl);
      vcs.atpos[ivcl].push_back(vcs.atposcl[ipc]);
      vcs.pmapinv[ivcl].push_back(pair<size_t,size_t>(icl,ipc));
      int ipvc=vcs.atpos[ivcl].size()-1;
      vcs.pmap[icl][ipc]=pair<size_t,int>(ivcl,ipvc);
//      if(ipc==0)
//      {
//        printf("atpos0   :%14.8lf%14.8lf%14.8lf\n",
//        vcs.atpos[ivcl][0].x,vcs.atpos[ivcl][0].y,vcs.atpos[ivcl][0].z);
//        cin.get();
//      }
    }  
  } // end of routine mapcl2vcl  

//=======================================================================//
     
template<typename GridT>
inline void initip0vcl(GridT& grid,VCStr& vcs)
{
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
//CC Testing number of particles in verlet cells
//  size_t npt=0;
//  for(size_t ivcl=0;ivcl<vcs.nvcl;ivcl++)
//  {
//    size_t npcli=vcs.pmapinv[ivcl].size();
//    npt+=npcli;
//  }
//  if(npt!=vcs.np)
//  {  
//    printf("BUG; np,npt %6zu%6zu\n",vcs.np,npt);
//    abort();
//  }  
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

  size_t ip=0;
  size_t ncntvcl=0;
  int8_t mp1[2]={-1,1};
//C Getting ip0vcl for central box verlet cells
  for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
  {
    ncntvcl++;
    size_t ivcl=vcs.icllst[ilst].first;
    int npcli=vcs.pmapinv[ivcl].size();
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
//CC Testing whether npcbcli < npcli
//    int npcbcli=vcs.icllst[ilst].second;
//    printf("ilst,ivcl,npcbcli,npcli :%8zu%8zu%8d%8d\n",ilst,ivcl,npcbcli,vcs.pmapinv[ivcl].size());
//    if(npcbcli<vcs.pmapinv[ivcl].size())
//    {
//      printf("ilst,ivcl,npcbcli,npcli :%8zu%8zu%8d%8d\n",ilst,ivcl,npcbcli,vcs.pmapinv[ivcl].size());
//      abort();
//    }
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
    assert(ivcl>=0&&ivcl<vcs.nvcl);
    vcs.ip0vcl[ivcl]=ip;
    for(int ipvc=0;ipvc<npcli;ipvc++,ip++)
    {
      vcs.ivclp[ip]=ivcl;
      size_t icl=vcs.pmapinv[ivcl][ipvc].first;
      size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
      vcs.ispc[ip]=static_cast<uint8_t>(grid.particle_type(icl,ipc));
      if(vcs.istatcl[icl]==4){vcs.istat[ip]=4;}
      else{vcs.istat[ip]=mp1[3+vcs.istatcl[icl]];}
    }
  }
  vcs.npcb=ip;
//  printf("vcs.npcb :%8zu\n",vcs.npcb);
//  printf("vcs.np   :%8zu\n",vcs.np);
//  std::cin.get();

//C Getting ip0vcl for verlet cells outside central box.
  IJK igls=vcs.mvclg;
  size_t nglmax=std::max(igls.i,std::max(igls.j,igls.k));
  for(size_t il=0;il<nglmax;il++)
  {
    igls=igls-1;
    ssize_t iglsi=std::max(static_cast<ssize_t>(0),igls.i);
    ssize_t iglsj=std::max(static_cast<ssize_t>(0),igls.j);
    ssize_t iglsk=std::max(static_cast<ssize_t>(0),igls.k);
    ssize_t idcli=1;
    ssize_t idclj=1;
    ssize_t idclk=1;
    if(igls.i<0&&igls.j<0){idclk=vcs.mvcl.k-2*iglsk-1;}
    for(ssize_t iclk=iglsk;iclk<vcs.mvcl.k-iglsk;iclk+=idclk)
    {
      if(igls.k>=0&&(iclk==iglsk||iclk==vcs.mvcl.k-iglsk-1)){idclj=1;}
      else if(igls.i<0){idclj=vcs.mvcl.j-2*iglsj-1;}
      for(ssize_t iclj=iglsj;iclj<vcs.mvcl.j-iglsj;iclj+=idclj)
      {
        if(igls.k>=0&&(iclk==iglsk||iclk==vcs.mvcl.k-iglsk-1)){idcli=1;}
        else
        {
          if(igls.j>=0&&(iclj==iglsj||iclj==vcs.mvcl.j-iglsj-1)){idcli=1;}
          else{idcli=vcs.mvcl.i-2*iglsi-1;}
        }
        for(ssize_t icli=iglsi;icli<vcs.mvcl.i-iglsi;icli+=idcli)
        {
          IJK iclijk{icli,iclj,iclk};
          ncntvcl++;
          size_t ivcl=grid_ijk_to_index(vcs.mvcl,iclijk);
          assert(ivcl>=0&&ivcl<vcs.nvcl);
          vcs.ip0vcl[ivcl]=ip;
          int npcli=vcs.pmapinv[ivcl].size();
          for(int ipvc=0;ipvc<npcli;ipvc++,ip++)
          {
            if(ip>=vcs.np)
            {
              printf("BUG; ip,vcs.np  :%8zu%8zu\n",ip,vcs.np);
              abort();
            }
            assert(ip<vcs.np);
            vcs.ivclp[ip]=ivcl;
            size_t icl=vcs.pmapinv[ivcl][ipvc].first;
            size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
            vcs.ispc[ip]=static_cast<uint8_t>(grid.particle_type(icl,ipc));
            vcs.istat[ip]=vcs.istatcl[icl];
          }
        }
      }  
    }  
    if(il==1){vcs.npcbp2=ip;}
  }
  if(vcs.iprint>0)
  {
    printf("Counted central box particles; vcs.np :%8zu\n",vcs.npcb);
    printf("Counted particles; ip,vcs.np :%8zu%8zu\n",ip,vcs.np);
    printf("Counted cells; ncntvcl,vcs.nvcl :%8zu%8zu\n",ncntvcl,vcs.nvcl);
  }
  assert(ip==vcs.np);

} // end of routine initip0vcl

//==================================================================//
//C Function getvclpcl gets the verlet cell index of all particles
//C in the exaStamp cell icl, and puts them in the list vcs.ivclpcl.
//C It also puts the atomic position wrt the origin vcs.xyz0
//C in a temporal array vsc.atposcl.
  
  template<typename CellT>
  inline void getvclpcl(CellT& cells,size_t icl,size_t npcli,VCStr& vcs)
  {
    double* pri[3]={cells[icl][field::rx],cells[icl][field::ry],cells[icl][field::rz]};
    ssize_t iclx,icly,iclz;

    for(size_t ipc=0;ipc<npcli;ipc++)
    {
      vcs.atposcl[ipc].z=pri[2][ipc]-vcs.xyz0.z;
      if(vcs.atposcl[ipc].z<0.){
        iclz=0;
      }else{
        iclz=floor(vcs.atposcl[ipc].z*vcs.dvcli.z);
        if(iclz>=vcs.mvcl.k) iclz=vcs.mvcl.k-1;
      }
//      ssize_t iclz=floor(vcs.atposcl[ipc].z*vcs.dvcli.z);   
//      assert(iclz>=0&&iclz<vcs.mvcl.k);
      vcs.ivclpcl[ipc]=iclz*vcs.mvclxy;
    }

    for(size_t ipc=0;ipc<npcli;ipc++)
    {
      vcs.atposcl[ipc].y=pri[1][ipc]-vcs.xyz0.y;
      if(vcs.atposcl[ipc].y<0.){
        icly=0;
      }else{
        icly=floor(vcs.atposcl[ipc].y*vcs.dvcli.y);
        if(icly>=vcs.mvcl.k) icly=vcs.mvcl.k-1;
      }
//      ssize_t icly=floor(vcs.atposcl[ipc].y*vcs.dvcli.y);
//      assert(icly>=0&&icly<vcs.mvcl.j);
      vcs.ivclpcl[ipc]+=icly*vcs.mvcl.i;
    }

    for(size_t ipc=0;ipc<npcli;ipc++)
    {
      vcs.atposcl[ipc].x=pri[0][ipc]-vcs.xyz0.x;
      if(vcs.atposcl[ipc].x<0.){
        iclx=0;
      }else{
        iclx=floor(vcs.atposcl[ipc].x*vcs.dvcli.x);
        if(iclx>=vcs.mvcl.k) iclx=vcs.mvcl.k-1;
      }
//      ssize_t iclx=floor(vcs.atposcl[ipc].x*vcs.dvcli.x);
//      assert(iclx>=0&&iclx<vcs.mvcl.i);
      vcs.ivclpcl[ipc]+=iclx;
    }  
  } // end of routine getvclpcl

//=======================================================================//
//C Function nbstblsizes gets the estimated minimal sizes of the tables
//C for the object nbs of type NBStr
     
  static inline void nbstblsizes(VCStr& vcs, double ntfact)
  {
    size_t ncnt=0;
    double pi=std::acos(-1.);
    double volvcl=vcs.dvcl.x*vcs.dvcl.y*vcs.dvcl.z;
    double volsph=4.*pi*pow(2.2,3)/3.0;  // HIER generalize
    double vsphbyvcl=volsph/volvcl;
    double rcmrm=4.;  // HIER generalize
    double rnnbav=0.;

//#  pragma omp parallel for
    for(size_t ivcl=0;ivcl<vcs.nvcl;ivcl++)
    {
      int npvcli=vcs.pmapinv[ivcl].size();
//      printf("npvcli   :%8zu\n",npvcli);
      if(npvcli>0)
      {
        ncnt++;
        if(npvcli>vcs.npvcm)vcs.npvcm=npvcli;
      }
    }
    rnnbav=vcs.np*vsphbyvcl/ncnt;
    double rnnbmrcnav=(4.*pi*pow(rcmrm,3)/(3.*volsph)-1.)*rnnbav;
    vcs.nnbm=std::ceil(rnnbav)-1;
    vcs.nnbmrcnm=std::ceil(rnnbmrcnav)-1;

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
//    cout << "2.2^3      = " << pow(2.2,3) << endl;
//    cout << "pi         = " << pi << endl;
//    cout << "rnnbav     = " << rnnbav << endl;
//    cout << "volsph     = " << volsph << endl;
//    cout << "rnnbmrcnav = " << rnnbmrcnav << endl;
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//

//C Determining sizes of ghost lists
    vcs.nglst1m=(vcs.mvclcb.j+2)*(vcs.mvclcb.k+2);
    vcs.nglst1m+=vcs.mvclcb.i*(vcs.mvclcb.k+2);
    vcs.nglst1m+=vcs.mvclcb.i*vcs.mvclcb.j;
//    vcs.nglst1m=2*vcs.npvcm*vcs.nglst1m;
//    vcs.nglst2m=2*vcs.nglst1m;
    vcs.nglst1m=ntfact*vcs.npvcm*vcs.nglst1m;
    vcs.nglst2m=ntfact*vcs.nglst1m;

    if(vcs.iprint>0)
    {
      cout << "npvcm     = " << vcs.npvcm << endl;
      cout << "nnbm      = " << vcs.nnbm << endl;
      cout << "nnbmrcnm  = " << vcs.nnbmrcnm << endl;
      cout << "nglst1m   = " << vcs.nglst1m << endl;
      cout << "nglst2m   = " << vcs.nglst2m << endl;
//      cin.get();
    }

  } // end of routine nbstblsizes

//=======================================================================//
     
  template<typename GridT,typename CellT>
  inline void printatpos(int icnt,size_t ncl,GridT& grid,CellT& cells,VCStr& vcs)
  {
    FILE *ATPOS;
    ATPOS=fopen("TMPDAT/atpos.xyz","a");
    fprintf(ATPOS,"icnt :%6d\n",icnt);
    fprintf(ATPOS,"%18.12lf%18.12lf%18.12lf\n",vcs.xyz0.x,vcs.xyz0.y,vcs.xyz0.z);
    for(size_t icl=0;icl<ncl;icl++)
    {
      size_t npcli=cells[icl].size();
      fprintf(ATPOS,"icl,npcli,incb :%6zu%6zu%6d\n",icl,npcli,vcs.incb[icl]);
      for(size_t ipc=0;ipc<npcli;ipc++)
      {
        int ispci=static_cast<uint8_t>(grid.particle_type(icl,ipc));
        Vec3d ri={cells[icl][field::rx][ipc],cells[icl][field::ry][ipc],cells[icl][field::rz][ipc]};
        size_t ivcl=vcs.pmap[icl][ipc].first;
        int ipvc=vcs.pmap[icl][ipc].second;
        size_t ip=vcs.ip0vcl[ivcl]+ipvc;
        if(ispci==0){fprintf(ATPOS,"C  %18.12lf%18.12lf%18.12lf%6zu%6d%6d%6d\n",ri.x,ri.y,ri.z,ivcl,ipvc,vcs.ispc[ip],vcs.istat[ip]);}
        else{fprintf(ATPOS,"H  %18.12lf%18.12lf%18.12lf%6zu%6d%6d%6d\n",ri.x,ri.y,ri.z,ivcl,ipvc,vcs.ispc[ip],vcs.istat[ip]);}
        Vec3d dri=vcs.atpos[ivcl][ipvc]+vcs.xyz0-ri;
        fprintf(ATPOS,"%18.12lf%18.12lf%18.12lf\n",dri.x,dri.y,dri.z);
      }  
    }  
    fclose(ATPOS);
    printf("printatpos done\n");
  } // end of routine printatpos

//======================================================================//
