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

//#include <memory>
//#include <iostream>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <fstream>
#include <iomanip>

#include <exanb/core/grid.h>

#include "lchbop_getEF.h"

using std::vector;


void fcpdown(double,double,double,double,double,double&,double&);
void sqdsqdown0(double,double,double,double,double&,double&);
void SqdSqup(double,double,double,double,double&,double&);
void SqdSqpup(double,double,double,double,double,double&,double&);
void fcpdownrc(double,double,double,double,double,double&,double&,double&);

template<typename T>
void wrtptbl(size_t,int,int,const vector<T>&);

template<typename T,size_t N>
void wrtptbl(size_t,int,int,T (&table)[N]);

//template<typename T>
//void wrtvectbl(size_t,vector<Vec3d>&);

template<typename CellT>
void wrtSRT(CellT&,VCStr&,NBStr&);

inline void wrtnbtip(size_t,size_t,VCStr&,NBStr&);

//======================================================================//

template<typename T>
T typeconv(T x){return x;}  
inline int typeconv(int8_t x){return static_cast<int>(x);}  
inline int typeconv(uint8_t x){return static_cast<int>(x);}  

//======================================================================//
     
template<typename T,size_t N>
void wrtptbl(size_t n,int iw1,int iw2,T (&table)[N])
{
  wrtptbl(n,iw1,iw2,vector<std::remove_const_t<T>>(table,table+N));
}

//======================================================================//
     
  template<typename T>
  void wrtptbl(size_t n,int iw1,int iw2,const vector<T>& table)
{
  std::ofstream TBLDAT("TMPDAT/tbl.dat",std::ios_base::out);
//    std::ofstream TBLDAT("TMPDAT/tbl.dat",std::ios_base::app);
//    TBLDAT << std::setprecision(6);
//    TBLDAT << std::setprecision(15);
  assert(table.size()>=n);
  for(size_t i=0;i<n;i++){TBLDAT << std::setw(iw1) << std::setprecision(iw1) << i << "    " 
                                 << std::setw(iw2) << std::setprecision(iw2) << typeconv(table[i]) << '\n';}  
//    for(size_t i=0;i<n;i++){TBLDAT << i << '\t' << typeconv(table[i]) << '\n';}  
} // end of routine wrtptbl

//======================================================================//

template<typename CellT>
void wrtSRT(CellT& cells,VCStr& vcs,NBStr& nbs)
{
  FILE *SRDAT;
  SRDAT=fopen("TMPDAT/srtcl.dat","w");
  fclose(SRDAT);

//C Mapping content of exaStamp cells inside central box to verlet cells.
  for(ssize_t iclk=vcs.ngls0;iclk<vcs.dims.k-vcs.ngls0;iclk++)
  for(ssize_t iclj=vcs.ngls0;iclj<vcs.dims.j-vcs.ngls0;iclj++)
  for(ssize_t icli=vcs.ngls0;icli<vcs.dims.i-vcs.ngls0;icli++)
  {
    IJK iclijk{icli,iclj,iclk};
    size_t icl=grid_ijk_to_index(vcs.dims,iclijk);
    assert(icl>=0&&icl<vcs.ncl);
    size_t npcli=cells[icl].size();
//    printf("icl,npcli :%8zu%8zu\n",icl,npcli);
    for(size_t ipc=0;ipc<npcli;ipc++)wrtnbtip(icl,ipc,vcs,nbs);
  }  
//    cin.get();

  GridBlock cbox;
  cbox.start={vcs.ngls0,vcs.ngls0,vcs.ngls0};
  cbox.end={vcs.dims.i-vcs.ngls0,vcs.dims.k-vcs.ngls0,vcs.dims.k-vcs.ngls0};

  for(size_t icl=0;icl<vcs.ncl;icl++)
  {
    IJK ijkcl=grid_index_to_ijk(vcs.dims,static_cast<ssize_t>(icl));
    if(inside_block(cbox,ijkcl)){vcs.istatcl[icl]=4;continue;}
    size_t npcli=cells[icl].size();
//    printf("icl,npcli,npcnt :%8zu%8zu%8zu\n",icl,npcli,npcnt);
    for(size_t ipc=0;ipc<npcli;ipc++)wrtnbtip(icl,ipc,vcs,nbs);
  }  
//    cin.get();
}

//======================================================================//

inline void wrtnbtip(size_t icl,size_t ipc,VCStr& vcs,NBStr& nbs)
{
  size_t ivcl=vcs.pmap[icl][ipc].first;
  size_t ipvc=vcs.pmap[icl][ipc].second;
  size_t ip0=vcs.ip0vcl[ivcl];
  size_t ip=ip0+ipvc;
  if(nbs.nnbfrac[ip]>nbs.nnb0[ip])
  {
    FILE *SRDAT;
    SRDAT=fopen("TMPDAT/srtcl.dat","a");
    fprintf(SRDAT,"%8zu%8zu%8zu%8zu%8zu\n",icl,ipc,nbs.nnbfull[ip]-nbs.nnb0[ip],
    nbs.nnbfrac[ip]-nbs.nnbfull[ip],nbs.nnbfrac[ip]-nbs.nnb0[ip]);
    for(size_t inb=nbs.nnb0[ip];inb<nbs.nnbfrac[ip];inb++)
    {
      size_t jp=nbs.nbtbl[inb];
      size_t jvcl=vcs.ivclp[jp];
      size_t jp0=vcs.ip0vcl[jvcl];
      int jpvc=jp-jp0;
      size_t jcl=vcs.pmapinv[jvcl][jpvc].first;
      size_t jpc=vcs.pmapinv[jvcl][jpvc].second;
//      fprintf(SRDAT,"[%6zu",jcl,"][%6zu",jpc,"]   ");
      fprintf(SRDAT,"[%6zu%2s%6zu%4s",jcl,"][",jpc,"]   ");
    }  
    fprintf(SRDAT,"\n\n");
    fclose(SRDAT);
  }
}

//======================================================================//
//
//  void printMRT0(size_t ip1,size_t ip2,NBStr& nbs)
//  {
//    printf("Pure MR bonds\n");
//    for(size_t ip=ip1;ip<ip2;ip++)
//    {
//      printf("ip,nnbmr  :%8zu%8zud\n",ip,nbs.nnbmr[ip]);
//      size_t imrij=nbs.indmr[ip];
//      for(size_t inb0=0;inb0<nbs.nnbmr[ip];inb0++)
//      {
//        size_t jp=nbs.nbtblmr[imrij];
//        size_t imrji=nbs.indmrji[imrij];
//        printf("inb0,ip,jp,imrij,imrji :%6zud%6zud%6zud%6zud%6zud\n",inb0,ip,jp,imrij,imrji);
//        imrij=nbs.nextind[imrij];
//      }  
//    }
//  }
//
//======================================================================//
