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

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/algorithm.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/physics/units.h>
#include <exanb/core/particle_id_codec.h>

#include <memory>
#include <iostream>
//#include <mpi.h>
//#include <sys/time.h>

#include "lchbop_inits.h"
#include "lchbop_getEF.h"
#include "lchbop_utils.h"
#include "lchbop_common.h"
//#include "parallel_build_dual_nbtbl.h"

using std::vector;
using std::abs;
using std::sqrt;
using std::pow;
using std::cout;
using std::cin;
using std::endl;
//using std::clock;

inline void testvalue(size_t, vector<uint8_t>&);

template<typename GridT>
void getEF(VLCHBOP&,vector<int>&,GridT&,VCStr&,NBStr&);

template<typename CellT>
inline void chkpmap(CellT&,VCStr&);

inline void chkistat(VCStr&);

template<typename GridT,typename CellT>
inline void vcsupdate(GridT&,CellT&,VCStr&);

template<typename CellT>
inline void chkvcsupdate(CellT&,VCStr&);

template<typename CellT>
inline void passEFV(CellT&,VCStr&,NBStr&);

inline void vclmigr(size_t,size_t,size_t,size_t,int,VCStr&);
inline void rmip1(size_t,int,VCStr&);
inline void rmip2(size_t,size_t,int,VCStr&);
inline void addip1(size_t,size_t,size_t,VCStr&);
inline void addip2(size_t,size_t,size_t,size_t,VCStr&);
inline void getpmap(size_t,int,pair<size_t,size_t>,vector<vector<pair<size_t,int>>>&);

void getnbcbp1(VLCHBOP&,VCStr&,NBStr&,GhostTbls&);
void getnbcbp2(VLCHBOP&,VCStr&,NBStr&,GhostTbls&);
void getnbcbp3(VLCHBOP&,VCStr&,NBStr&);

inline bool activebond(int,int);
inline pair<size_t,size_t> iclipc(size_t,VCStr&);
inline size_t ipvcl(size_t,int,VCStr&);

inline void getnbgp1(VLCHBOP&,VCStr&,NBStr&,GhostTbls&);
inline void getnbgp2(VLCHBOP&,VCStr&,NBStr&,GhostTbls&);
inline void getnbsgp(VLCHBOP&,size_t,int,int,size_t,size_t,VCStr&,NBStr&);

inline void getcbnb(VLCHBOP&,size_t,int,int,size_t,size_t,size_t,
const Vec3d&,VCStr&,NBStr&);  //  pourquoi "const" devant Vec3d ????? //
inline void getgnb1(VLCHBOP&,size_t,int,size_t,size_t,int,int,
const Vec3d,VCStr&,NBStr&); 
inline void getgnb2(VLCHBOP&,size_t,int,size_t,size_t,int,int,
const Vec3d&,VCStr&,NBStr&);
inline void getgnb3(VLCHBOP&,size_t,int,size_t,size_t,int,
const Vec3d,double*,VCStr&,NBStr&);

void getrsig(double&,Vec3d&,Vec3d&);
void getsn(VLCHBOP&,int,int,double,double*);
inline void add2srtbl1(size_t,size_t,int,size_t,double,const Vec3d,double*,NBStr&);
inline void add2srtbl2(size_t,size_t,int,size_t,NBStr&);

void getinbji(size_t,size_t,size_t&,NBStr&);
void getndbap(VLCHBOP&,VCStr&,NBStr&,GhostTbls&);
void getndbc1(VLCHBOP&,NBStr&,size_t);
void getndbc2(VLCHBOP&,size_t,int,vector<uint8_t>&,NBStr&);
inline void getmrtbls1(VLCHBOP&,VCStr&,NBStr&,GhostTbls&);
inline void getndbij(VLCHBOP&,size_t,size_t,int,int,double,double&,NBStr&);
void submr(VLCHBOP&,int,size_t,size_t,int,int,int,double,const Vec3d,
double,double,int,int&,double,VmrLocDat&,vector<uint8_t>&,NBStr&);
inline void getmrnb0(VLCHBOP&,size_t,int,int,double,VCStr&,NBStr&);
inline void getmrnb1(VLCHBOP&,size_t,int,int,double,VCStr&,NBStr&);
inline void getmrnb2(VLCHBOP&,size_t,int,double,VCStr&,NBStr&,GhostTbls&);
inline void getmrnbij(VLCHBOP&,size_t,size_t,int,int,int,int&,
double,double,const Vec3d,VCStr&,NBStr&);

void getsumk(size_t,size_t,const Vec3d,double&,NBStr&);
void dndbij(VLCHBOP&,size_t,size_t,int,int,vector<uint8_t>&,double,NBStr&);
void dsumk(size_t,size_t,size_t,int,double,const Vec3d,NBStr&);
double fnelki(VLCHBOP&,double);
void fdfnelki(VLCHBOP&,double,double&,double&);
inline void loadnbssp(size_t,size_t,vector<uint8_t>&,LNBStr&,NBStr&);
void loadsrcmr(size_t,NBStr&,LNBStr&);
void vsmrij(VLCHBOP&,int,int,int,size_t,size_t,vector<uint8_t>&,
LNBStr*,FATPLoc&,NBStr&);
void subvsmrij(VLCHBOP&,int,int,int,size_t,size_t,double,double,
double&,double&,double&,double&,LNBStr*,NBStr&);
void subbij(VLCHBOP&,const uint8_t,int,int,size_t,double,
Vec3d,double&,LNBStr&,NBStr&);
void frcbijsym(int,int,size_t,double,LNBStr&,NBStr&);
void frcbijasym(int,int,size_t,size_t,size_t,double,LNBStr&,NBStr&);
void getrcmreff(LNBStr*);
void dsrcmreff(double,double,double,double,double&,LNBStr*);

inline void addvlr(size_t,size_t,double,Vec3d,Vec3d,NBStr&);

inline void chkaddvlr(int,size_t,size_t,double,Vec3d,double,VCStr&,NBStr&);

void wrtatpos(VCStr&);
//template<typename GridT,typename CellT>
//void wrtatpos(GridT&,CellT&,bool,size_t,size_t,vector<size_t>&,
//vector<size_t>&,vector<uint8_t>&,vector<Vec3d>&);

template<typename GridT,typename CellT>
void printCellDat(GridT&,CellT&,size_t,vector<size_t>&,vector<size_t>&);
     
template<typename GridT,typename CellT>
inline void chkcells(GridT&,CellT&,size_t,size_t,size_t,vector<size_t>&,
vector<size_t>&);
     
template<typename GridT>
inline void mklstgp2cp(GridT&,VCStr&,GhostTbls&);
     
void chknbgps(size_t,NBStr&,GhostTbls&);
void chknbigp(size_t,NBStr&,GhostTbls&);

inline void chkmrgps(VLCHBOP&,vector<uint8_t>&,vector<int8_t>&,NBStr&,GhostTbls&);
//void chkmrgps(VLCHBOP&,vector<int8_t>&,NBStr& nbs,GhostTbls& gtbls);
inline void chkmrigp(VLCHBOP&,size_t,vector<uint8_t>&,NBStr&,GhostTbls&);

void prtFrcCntr(int,size_t,Vec3d);

void printEsr(size_t,int,double,vector<double>&,VCStr&);
void printElr(size_t,int,double,vector<double>&,VCStr&);
void printEF(size_t,int,double,double,vector<Vec3d>&,VCStr&);
//void printEF(size_t,double,double,vector<Vec3d>&);
inline void printSRT(size_t,size_t,int,VCStr& vcs,NBStr&);
void printMRT0(size_t,size_t,NBStr&);
void printMRT(size_t,size_t,NBStr&,GhostTbls&,VCStr& vcs);
void printVIR(size_t,vector<Mat3d>&);
void orderlstjp(int,vector<size_t>&,vector<size_t>&);



namespace exaStamp
{

//#ifdef VIRIEL
//  template<class CellT,bool HasVirialTemp>
//  static inline void passVirial(CellT& cells,size_t icl,size_t ipc,Mat3d& virial,std::integral_constant<bool,HasVirialTemp> has_viriel)
//  {
//    cells[icl][field::virial][ipc]=virial; // big computation
//  }

//  template<class CellT>
//  static inline void getVirial(CellT& cells,size_t icl,size_t ipc,Mat3d& virial,std::integral_constant<bool,true>)
//  {
//     cells[icl][field::virial][ipc]=virial; // big computation
//  }

//  template<class CellT>
//  static inline void compute_viriel_lowlovel(CellT cells , size_t cell_i, size_t particle_j, std::integral_constant<bool,true> )
//  {
//     cells[cell_i][field::viriel][particle_j] = Mat3d{}; // big computation
//  }

//  template<class CellT>
//  static inline void compute_viriel_lowlovel(CellT, size_t , size_t , std::integral_constant<bool,false> )
//  {
//  }

//  template<class CellT, bool HasVirialTemp>
//  static inline void compute_viriel(CellT cells , size_t cell_i, size_t particle_j , std::integral_constant<bool,HasVirialTemp> has_viriel)
//  {
//
//    if( has_viriel )
//    {
//      // what if has viriel
//    }
//    else
//    {
//      // what if not
//    }
//
//    compute_viriel_lowlovel( cells,cell_i,particle_j,has_viriel  );
//  }
//#endif    

  template<
	typename GridT
        ,  class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz >
	>
  class LCHBOP_EF : public OperatorNode
  {
// compile time constant indicating if grid has viriel field
//#ifdef VIRIEL
//    using has_virial_field_t = std::integral_constant<bool,GridHasField<GridT,field::_virial>::value>;
//    static constexpr has_virial_field_t has_virial_field{};
//#endif    

// ========= I/O slots =======================
//    ADD_SLOT( MPI_Comm  , mpi         , INPUT );
    ADD_SLOT( Domain    , domain      , INPUT , REQUIRED );
//    ADD_SLOT( bool      , flag_vmr    , INPUT );
    ADD_SLOT( VLCHBOP   , parameters  , INPUT );
    ADD_SLOT( double    , nbstbl_fact , INPUT , 2. );
    ADD_SLOT( GridT     , grid        , INPUT_OUTPUT );
    ADD_SLOT( VCStr     , vcs         , INPUT_OUTPUT );
    ADD_SLOT( NBStr     , nbs         , INPUT_OUTPUT );
    ADD_SLOT( NbCells   , nbclist     , INPUT );

    public:
    
    void execute () override final
    {
      GridT& grid = *(this->grid);
      VLCHBOP& par = *(this->parameters);
      double ntfact = *(this->nbstbl_fact);
      VCStr& vcs = *(this->vcs);
      NBStr& nbs = *(this->nbs);
//      bool vmrflag = *(this->flag_vmr);
      assert(vcs.np=grid.number_of_particles());

//      compute_viriel(cells,0,0,has_viriel_field);

// Get neighbours tables, energy and forces
      auto cells=grid.cells();

//      cout << "vcs.mkvcs  =" << vcs.mkvcs  << endl;
//      cout << "vcs.NPT    =" << vcs.NPT    << endl;

      if(vcs.mkvcs)
      {  
        vcs.hmat=domain->xform();
        vcs.xyz0=grid.cell_position(vcs.ijkcb);
        getvcsgrid(vcs);
      }
      else if(vcs.NPT)
      {
        vcs.hmat=domain->xform();
        vcs.xyz0=grid.cell_position(vcs.ijkcb);
        vcs.mkvcs=chkvcsgrid(vcs);
//        if(vcs.mkvcs){printf("Case does occur \n");cin.get();}
      }

//      printf("getEF; hmat :%16.10lf%16.10lf%16.10lf\n",vcs.hmat.m11,vcs.hmat.m22,vcs.hmat.m33);
//      cin.get();

      if(vcs.mkvcs)
      {
        NbCells& nbc = *(this->nbclist);
        vcs.initvcs2(nbc.nnbc,nbc.nbclst);
        vclfill(cells,vcs);
        initip0vcl(grid,vcs);
        nbstblsizes(vcs,ntfact);
        nbs.initnbs1(vcs.np,vcs.npcbp2,vcs.nnbm,vcs.nnbmrcnm);
        vcs.mkvcs=false;
      }
      else
      {
        vcsupdate(grid,cells,vcs);
//        chkpmap(cells,vcs);
//        chkistat(vcs);
//        chkvcsupdate(cells,vcs);
//        printf("vcsupdate done \n");
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
// For testing only   
//        vcs.atposold=vcs.atpos;
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
      }  

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
// For testing only   
//      FILE *ATPOS;
//      int err,np;
//      Vec3d hxyz,atpos[41984];
//      char texte[7];
//      ATPOS=fopen("PolCEA.xyz","r");
//      err=fscanf(ATPOS,"%d\n",&np);
//      err=fscanf(ATPOS,"%lf%lf%lf\n",&hxyz.x,&hxyz.y,&hxyz.z);
//      printf("np   :%6d\n",np);
//      printf("hxyz :%16.10lf%16.10lf%16.10lf\n",hxyz.x,hxyz.y,hxyz.z);
//      for(int ip=0;ip<np;ip++)
//      {
//        err=fscanf(ATPOS,"%s%lf%lf%lf\n",texte,&atpos[ip].x,&atpos[ip].y,&atpos[ip].z);
////        printf("ip,atpos :%6d%16.10lf%16.10lf%16.10lf\n",ip,atpos[ip].x,atpos[ip].y,atpos[ip].z);
////        cin.get();
//        int iflag=0;
//        for(int ipp=0;ipp<np;ipp++)
//        {
////          printf("ip,ipp :%6d%6d\n",ip,ipp);
//          size_t ivcl=vcs.ivclp[ipp];
//          int ipvc=ipp-vcs.ip0vcl[ivcl];
//          Vec3d drip=vcs.atpos[ivcl][ipvc]+vcs.xyz0;
//          drip=drip-atpos[ip];
//          double dripsq=norm2(drip);
//          if(dripsq<1e-12)
//          {
//            iflag=1;
//            vcs.ipmap[ip]=ipp;
//            vcs.ipmapinv[ipp]=ip;
//            nbs.ipmap[ip]=ipp;
//            nbs.ipmapinv[ipp]=ip;
//            break;}
//        }
//        if(iflag==0){printf("BUG: ipp not found: ip =%6\n",ip);abort();}
//      }  
//      fclose(ATPOS);
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
//      cin.get();

      getEF(par,grid,vcs,nbs);
      passEFV(cells,vcs,nbs);

//      printf("Operator LCHBOP_EF done\n");
    }

  };

  template<class GridT> using LCHBOP_EF_Tmpl = LCHBOP_EF<GridT>;

// === register factories ===  
  __attribute__((constructor)) static void register_factories()
  {
    OperatorNodeFactory::instance()->register_factory( "LCHBOP_EF",
    make_grid_variant_operator< LCHBOP_EF_Tmpl > );
  }

} // end namespace exaStamp

//======================================================================//

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//
// General comments and meaning of crucial variables and arrays
// Variables:
//  - vcs.ncl    = total number of exaStamp cells
//  - vcs.np     = total number of atoms
//  - vcs.ngls0  = number of ghost cell layers in the exaStamp cell structure
//  - vcs.nvcl   = total number of verlet cells
//  - vcs.nvclcb = number of verlet cells within the verlet structure central box (VSCB)
//  - vcs.npcb   = number of atoms in VSCB
// Arrays:
//  - ip0cl[ncl]    = contains the number of the first atoms in each verlet cell
//  - icllst[nclcb] = contains the cell indices of all verlet inside central box 
//  - incb[ncl]     = contains labels for each cell; 
//      when incb[icl] > 0, cell icl is inside central box:
//        incb[icl] = 1 --> cell icl is in most outer layer of central box
//        incb[icl] = 2 --> cell icl is in second most outer layer of central box
//        incb[icl] = 3 --> cell icl is in third most outer layer of central box
//        incb[icl] = 4 --> cell icl is in the inner part of central box
//      when incb[icl] < 0, cell icl is outside central box:
//        incb[icl] =-1 --> cell icl is in active ghost zone
//        incb[icl] =-2 --> cell icl is in non-active ghost zone
//        (interactions of a central box atoms with ghost atoms in an active zone
//         have to computed, while those with ghost atoms in an inactive zone
//         should non be computed but are computed on a different processor)
//  - iflag[np] = label for (computational) status of each atom
//      istat[icl][ip] = 4 --> ip is in the ExaStamp Central Box (ESCB); 
//                             if flag[ip] is not equal to 4 then ip is a ghost atom  
//      istat[icl][ip] =-3 --> ip is ghost in the non-active zone and SR nor MR neighbor 
//                             tables are not (yet) computed/available
//      istat[icl][ip] =-2 --> ip is ghost in the active zone and SR nor MR neighbor 
//                             tables are not (yet) computed/available
//      istat[icl][ip] =-1 --> ip is atom in the non-active ghost zone and SR or candidate
//                             neighbor atom of one (or more) atoms(s) in the central box; 
//                             SR and MR candidate are listed to be computed                            
//      istat[icl][ip] = 0 --> ip is an atom in the active ghost zone and candidate 
//                             pure MR neighbor of a central box atom;
//                             at the moment that ip turns out to be an active pure MR neighbor 
//                             of a central box atom, then iflag[ip] is set to 1;
//                             SR and MR candidate are listed to be computed                            
//      istat[icl][ip] = 1 --> ip is an atom in the active ghost zone and SR and/or active
//                             MR neighbor of one (or more) atom(s) in the central box;
//                             SR and MR candidate are listed to be computed
//      istat[icl][ip] = 2 --> ip is an atom in the active ghost zone and neighbor of a
//                             ghost atom in the active zone that is a SR or MR neighbor; 
//                             of an atom in the central box;
//                             SR and MR candidate are listed to be computed
//      istat[icl][ip] = 3 --> ip is an atom in the non-active ghost zone and neighbor of 
//                             a ghost atom in the active zone that is a SR or MR neighbor 
//                             of an atom in the central box;
//                             SR and MR candidate are listed to be computed
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC//

  template<typename GridT>
  void getEF(VLCHBOP& par,GridT& grid,VCStr& vcs,NBStr& nbs)
  {
//    auto cells=grid.cells();    // should not be needed anymore after testing

    nbs.initnbs2(vcs.np);
    GhostTbls gtbls(vcs.np,vcs.nglst1m,vcs.nglst2m); 

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
// For testing
//    vcs.ivclt=vcs.pmap[vcs.iclt][vcs.ipct].first,vcs.ipvct=vcs.pmap[vcs.iclt][vcs.ipct].second;
//    vcs.jvclt=vcs.pmap[vcs.jclt][vcs.jpct].first,vcs.jpvct=vcs.pmap[vcs.jclt][vcs.jpct].second;
//    vcs.kvclt=vcs.pmap[vcs.kclt][vcs.kpct].first,vcs.kpvct=vcs.pmap[vcs.kclt][vcs.kpct].second;
//    vcs.lvclt=vcs.pmap[vcs.lclt][vcs.lpct].first,vcs.lpvct=vcs.pmap[vcs.lclt][vcs.lpct].second;
//    vcs.ipt=ipvcl(vcs.ivclt,vcs.ipvct,vcs);
//    vcs.jpt=ipvcl(vcs.jvclt,vcs.jpvct,vcs);
//    vcs.kpt=ipvcl(vcs.kvclt,vcs.kpvct,vcs);
//    vcs.lpt=ipvcl(vcs.lvclt,vcs.lpvct,vcs);
//    nbs.iclt=vcs.iclt; nbs.ipct=vcs.ipct; nbs.ipt=vcs.ipt; 
//    nbs.jclt=vcs.jclt; nbs.jpct=vcs.jpct; nbs.jpt=vcs.jpt; 
//    nbs.kclt=vcs.kclt; nbs.kpct=vcs.kpct; nbs.kpt=vcs.kpt; 
//    nbs.lclt=vcs.lclt; nbs.lpct=vcs.lpct; nbs.lpt=vcs.lpt; 
//    printf("iclt,ipct,ivclt,ipvct,ipt :%6zu%6zu%6zu%6d%6zu\n",vcs.iclt,vcs.ipct,vcs.ivclt,vcs.ipvct,vcs.ipt);
//    printf("jclt,jpct,jvclt,jpvct,jpt :%6zu%6zu%6zu%6d%6zu\n",vcs.jclt,vcs.jpct,vcs.jvclt,vcs.jpvct,vcs.jpt);
//    printf("kclt,kpct,kvclt,kpvct,kpt :%6zu%6zu%6zu%6d%6zu\n",vcs.kclt,vcs.kpct,vcs.kvclt,vcs.kpvct,vcs.kpt);
//    printf("lclt,lpct,lvclt,lpvct,lpt :%6zu%6zu%6zu%6d%6zu\n",vcs.lclt,vcs.lpct,vcs.lvclt,vcs.lpvct,vcs.lpt);
//    printf("istatip,istatjp,istatkp,istatlp :%6d%6d%6d%6d\n",
//    vcs.istat[vcs.ipt],vcs.istat[vcs.jpt],vcs.istat[vcs.kpt],vcs.istat[vcs.lpt]);
//    printf("incbip,incbjp,incbkp,incblp :%6d%6d%6d%6d\n",
//    vcs.incb[vcs.ivclt],vcs.incb[vcs.jvclt],vcs.incb[vcs.kvclt],vcs.incb[vcs.lvclt]);
//    IJK ijkipt=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(vcs.ivclt));
//    IJK ijkjpt=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(vcs.jvclt));
//    IJK ijkkpt=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(vcs.kvclt));
//    IJK ijklpt=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(vcs.lvclt));
//    printf("ipt,ivclt,ijkipt :%6zu%6zu%6d%6d%6d\n",vcs.ipt,vcs.ivclt,ijkipt.i,ijkipt.j,ijkipt.k);
//    printf("jpt,jvclt,ijkjpt :%6zu%6zu%6d%6d%6d\n",vcs.jpt,vcs.jvclt,ijkjpt.i,ijkjpt.j,ijkjpt.k);
//    printf("kpt,kvclt,ijkkpt :%6zu%6zu%6d%6d%6d\n",vcs.kpt,vcs.kvclt,ijkkpt.i,ijkkpt.j,ijkkpt.k);
//    printf("lpt,lvclt,ijklpt :%6zu%6zu%6d%6d%6d\n",vcs.lpt,vcs.lvclt,ijklpt.i,ijklpt.j,ijklpt.k);
//    size_t ipt=224;
//    size_t jpt=758;
//    printf("ipt,icl,ipc,istatip :%6zu%6zu%6zu%6d\n",ipt,iclipc(ipt,vcs).first,iclipc(ipt,vcs).second,vcs.istat[ipt]);
//    printf("jpt,jcl,jpc,istatjp :%6zu%6zu%6zu%6d\n",jpt,iclipc(jpt,vcs).first,iclipc(jpt,vcs).second,vcs.istat[jpt]);
////    printf("vcl,ipvc,npcli,npcbcli  :%6d%6d%6d\n",vcs.ivclt,vcs.ipvct,vcs.pmapinv[ivcl].size(),vcs.icllst[ilst].second);
//    cin.get();
//    nbs.ipt=vcs.ipt;
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

//C Getting SR and MR candidate neighbors of central (box) atoms
//C up to the first layer of verlet cells
    getnbcbp1(par,vcs,nbs,gtbls);
//    printf("getnbcbp1 done \n");

//    wrtptbl(vcs.np,6,6,vcs.ispc);
//    chkcells(grid,cells,nclcb,npcb,np,icllst,ip0cl);
//    printSRT(vcs.ipt,vcs.ipt+1,0,vcs,nbs);
//    printSRT(0,vcs.npcb,1,vcs,nbs);
//    printMRT(202,203,nbs,gtbls);
//    abort();

#ifdef TESTEnFrcVir
//C Getting central box images of ghost particles for comparison to fortran code
//C for a periodic system using just one mpi-process
    mklstgp2cp(grid,vcs,gtbls);
#endif    

//C Get SR and MR candidate neighbors of the atoms in ghostlist 1
    if(gtbls.nglst1>0){getnbgp1(par,vcs,nbs,gtbls);}
//    printf("getnbgp1 done \n");
//    printSRT(vcs.ipt,vcs.ipt+1,1,vcs,nbs);
//    printSRT(0,vcs.np,1,vcs,nbs);

//C Get SR and MR candidate neighbors of the atoms in ghostlist 2
    if(gtbls.nglst2>0){getnbgp2(par,vcs,nbs,gtbls);}
//    printf("getnbgp2 done \n");
//    printSRT(0,vcs.np,1,vcs,nbs);
//    printSRT(0,vcs.npcb,1,vcs,nbs);

//C Verify neighbor tables of ghost particles for periodic system 
//    chknbgps(npcb,nbs,gtbls);

//C Getting dangling bond numbers of all atoms in central box and ghostlist 1
    getndbap(par,vcs,nbs,gtbls);
//    printf("getndbap done \n");

//C Get active MR neighbor tables from the MR candidate neighbor lists
//C for all central atoms and first neighbor ghost atoms (ie those in glst1)
    getmrtbls1(par,vcs,nbs,gtbls);
//    printf("getmrtbls1 done \n");

//C Get pure MR and LR neighbors within the second verlet layer of 
//C all atoms in the central box
    getnbcbp2(par,vcs,nbs,gtbls);
//    printf("getnbcbp2 done \n");
//    printf("nmr %6zu\n",nbs.nmr);
//    printMRT(23,24,nbs,gtbls,vcs);
//    abort();

//C Verifying MR tables for ghostparticles in glst1 
//    chkmrgps(par,vcs.ispc,vcs.istat,nbs,gtbls);

//C Find LR neighbors within the third verlet layer of all atoms 
//C in the central box
    getnbcbp3(par,vcs,nbs);
//    printf("getnbcbp3 done \n");

//    double elr=0.;
//    for(size_t ip=0;ip<vcs.np;ip++){elr+=nbs.elrpp[ip];}
//    printf("elr = %18.10lf\n",elr/vcs.np);
//    printElr(vcs.np,1,elr,nbs.elrpp,vcs);
//    abort();
 
//C Calculate energy and forces
    size_t nnbmloc=6*vcs.nnbm;          // Check for dynamical optimizes adjustment
    LNBStr nbsij[2]={nnbmloc,nnbmloc};
 
    nbs.nnbsrm+=1;
    nbsij[0].nnbsrm=nbs.nnbsrm;   // for testing only
    nbsij[1].nnbsrm=nbs.nnbsrm;   // for testing only
    nbsij[0].nbrm=400;            // for testing only
    nbsij[1].nbrm=400;            // for testing only
    nbsij[0].nnbmloc=nnbmloc;     // for testing only
    nbsij[1].nnbmloc=nnbmloc;     // for testing only

    FATPLoc Ftbls(nbs.nnbsrm);

    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
    {
      size_t icl=vcs.icllst[ilst].first;
      int npcli=vcs.icllst[ilst].second;
      if(npcli==0)continue; 
      size_t ip0=vcs.ip0vcl[icl];
      for(int ipc=0;ipc<npcli;ipc++) 
      {
        size_t ip=ip0+ipc;
        size_t imrij=0,imrji=0;
        int ifracnb=0,ipmrnb=0,imrnb=0;
        size_t imr=nbs.indmr[ip];
        size_t imrip=imr;
       
        loadnbssp(ip,imr,vcs.ispc,nbsij[0],nbs);
        
        nbsij[0].inbj=1;
        for(;nbsij[0].inbj<nbsij[0].nnbfl;nbsij[0].inbj++) // full neighbors
        {
          size_t jp=nbsij[0].nbtbl[nbsij[0].inbj];
          int istatjp=vcs.istat[jp];
          if(istatjp==4){if(jp<ip)continue;}
          else if(istatjp<=0)continue;
          imr=nbs.indmr[jp];
          loadnbssp(jp,imr,vcs.ispc,nbsij[1],nbs);
          nbsij[1].inbj=nbs.indji[nbs.nnb0[ip]+nbsij[0].inbj-1]-nbs.nnb0[jp]+1;
          vsmrij(par,ifracnb,ipmrnb,imrnb,imrij,imrji,vcs.ispc,nbsij,Ftbls,nbs);
        }  
        
        ifracnb=1;
        for(;nbsij[0].inbj<nbsij[0].nnbfr;nbsij[0].inbj++) // fractional neighbors
        {
          size_t jp=nbsij[0].nbtbl[nbsij[0].inbj];
          int istatjp=vcs.istat[jp];
          if(istatjp==4){if(jp<ip)continue;}
          else if(istatjp<=0)continue;
          imr=nbs.indmr[jp];
       
          loadnbssp(jp,imr,vcs.ispc,nbsij[1],nbs);
       
          nbsij[1].inbj=nbs.indji[nbs.nnb0[ip]+nbsij[0].inbj-1]-nbs.nnb0[jp]+1;
          size_t imrij=nbsij[0].isr2mr[nbsij[0].inbj];
          size_t imrji=nbsij[1].isr2mr[nbsij[1].inbj];
          if(imrij==0||imrji==0)
          {
            imrnb=0;
            nbsij[0].rcmreff=par.pRC[nbsij[0].ispc[0]][nbsij[1].ispc[0]].r2sr;
          }  
          else
          {
            loadsrcmr(imrij,nbs,nbsij[0]);
            loadsrcmr(imrji,nbs,nbsij[1]);
            if(nbsij[0].smr<1e-15||nbsij[1].smr<1e-15)
            {
              imrnb=0;
              nbsij[0].rcmreff=par.pRC[nbsij[0].ispc[0]][nbsij[1].ispc[0]].r2sr;
              printf("Exceptional case: MR neighbor listed while smrij/ji<1e-15\n");
//              cin.get();
            }  
            else
            {
              imrnb=1;
              getrcmreff(nbsij);
            }  
          }  
          vsmrij(par,ifracnb,ipmrnb,imrnb,imrij,imrji,vcs.ispc,nbsij,Ftbls,nbs);
        }  
       
        ipmrnb=imrnb=1;
        imrij=imrip;
        for(;nbsij[0].inbj<nbsij[0].nnb;nbsij[0].inbj++) // pure MR neighbors
        {
          loadsrcmr(imrij,nbs,nbsij[0]);
          if(nbsij[0].smr<1e-15){imrij=nbs.nextind[imrij];continue;}
          size_t jp=nbsij[0].nbtbl[nbsij[0].inbj];
          imrji=nbs.indmrji[imrij];
          if(imrji==0){imrij=nbs.nextind[imrij];continue;}
          int istatjp=vcs.istat[jp];
          if(istatjp==4){if(jp<ip){imrij=nbs.nextind[imrij];continue;}}
          else if(istatjp<=0){imrij=nbs.nextind[imrij];continue;}
          loadnbssp(jp,imrji,vcs.ispc,nbsij[1],nbs);
          loadsrcmr(imrji,nbs,nbsij[1]);
          if(nbsij[1].smr<1e-15){imrij=nbs.nextind[imrij];continue;}
          nbsij[1].inbj=nbsij[1].nnbfr;
          getrcmreff(nbsij);
          if(nbsij[0].rstor[nbsij[0].inbj]>nbsij[0].rcmreff){imrij=nbs.nextind[imrij];continue;}
          vsmrij(par,ifracnb,ipmrnb,imrnb,imrij,imrji,vcs.ispc,nbsij,Ftbls,nbs);
          imrij=nbs.nextind[imrij];
        }
      }
    }
    
//    if(vcs.icnt==3)
//    {
//      double esr=0.,elr=0.;
//      for(size_t ip=0;ip<vcs.np;ip++){elr+=nbs.elrpp[ip];esr+=nbs.esrpp[ip];}
//      printEsr(vcs.np,0,esr,nbs.esrpp,vcs);
//      printElr(vcs.np,0,elr,nbs.elrpp,vcs);
//      printEF(vcs.np,1,esr,elr,nbs.force,vcs);
//    }
//    abort();

#ifdef TESTEnFrcVir
//C Loop for adding ghost contributions for testing periodic system 
//C in comparison to fortran code
    for(size_t igp=vcs.npcb;igp<vcs.np;igp++)   // would it be better to use separate loops ????? //
    {
      size_t ip=gtbls.gp2cp[igp];
      nbs.elrpp[ip]+=nbs.elrpp[igp];
      nbs.esrpp[ip]+=nbs.esrpp[igp];
      nbs.force[ip]=nbs.force[ip]+nbs.force[igp];
//      nbs.viriel[ip]=nbs.viriel[ip]+nbs.viriel[igp];
      nbs.elrpp[igp]=0.0;
      nbs.esrpp[igp]=0.0;
      nbs.force[igp]={0.0,0.0,0.0};
//      nbs.viriel[igp]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    }  
    double esr=0.,elr=0.;
    for(size_t ip=0;ip<vcs.npcb;ip++){elr+=nbs.elrpp[ip];esr+=nbs.esrpp[ip];}
    printEsr(vcs.npcb,0,esr,nbs.esrpp,vcs);
    printElr(vcs.npcb,0,elr,nbs.elrpp,vcs);
    printEF(vcs.npcb,0,esr,elr,nbs.force,vcs);
//    printVIR(npcb,nbs.viriel);
    printf("Energie, forces and viriel written \n");
    abort();
#endif    

//    cout << "Finished !!" << endl;
    fflush(stdout);
  }  // end of function getEF

//======================================================================//
//C Determine whether a bond is active or not

  inline bool activebond(int istatip, int istatjp)
  {
    if(istatip==4)
    {
      if(istatjp!=-1&&abs(istatjp)!=3){return true;}
      else{return false;}
    }  
    else if(istatjp==4)
    {
      if(istatip!=-1&&abs(istatip)!=3){return true;}
      else{return false;}
    }  
    else{return false;}

  }   // end function activebond


//======================================================================//
//C testing values of a table

//  template<typename CellT>
  inline void testvalue(size_t n, vector<uint8_t>& table)
  {
    for(size_t i=0;i<n;i++)
    {
      if(table[i]<0)
      {
        printf("BUG; i,table[i] %6zu%6d\n",i,table[i]);
        abort();
      }
    }
  }   // end function testvalue


//======================================================================//
// Get third verlet layer LR neighbors of all atoms in the central box

  template<typename GridT,typename CellT>
  inline void vcsupdate(GridT& grid,CellT& cells,VCStr& vcs)
  {
//C Updating vcs for particles in each exacell

    for(size_t icl=0;icl<vcs.ncl;icl++)
    {
//      double dripsqmax=0.0;
      size_t npcli=cells[icl].size();
      getvclpcl(cells,icl,npcli,vcs);
      for(size_t ipc=0;ipc<npcli;ipc++)
      {
        size_t ivcl=vcs.pmap[icl][ipc].first;
        int ipvc=vcs.pmap[icl][ipc].second;
        size_t ivclnew=vcs.ivclpcl[ipc];
        if(ivcl<0||ivcl>=vcs.nvcl)
        {
          printf("BUG; icl,ipc      :%6zu%6zu\n",icl,ipc);
          printf("BUG; ivcl,ivclnew,vcs.nvcl :%6zu%6zu%6zu\n",ivcl,ivclnew,vcs.nvcl);
          cin.get();
        }  
        if(ivclnew==ivcl)
        {
//          Vec3d drip=vcs.atposcl[ipc]-vcs.atpos[ivcl][ipvc];
          vcs.atpos[ivcl][ipvc]=vcs.atposcl[ipc];
//          Vec3d drip=vcs.atpos[ivcl][ipvc]-vcs.atposold[ivcl][ipvc];
//          double dripsq=norm2(drip);
//          if(dripsq>dripsqmax)dripsqmax=dripsq;
        }  
        else
        {
//          printf("particle migration; icl,ipc,ivcl,ipvc,incb :%6zu%6zu%6zu%6d%6d\n",icl,ipc,ivcl,ipvc,vcs.incb[ivcl]);
//          cin.get();
          vclmigr(icl,ipc,ivcl,ivclnew,ipvc,vcs);  // particle has migrated
        }  
      }  
//      printf("icl,dripmax %6zu%16.10lf\n",icl,sqrt(dripsqmax));
    }   // end loop icl

    initip0vcl(grid,vcs);
//    cin.get();

  }   // end function vcsupdate

//======================================================================//
     
  inline void vclmigr(size_t icl,size_t ipc,size_t ivcl,size_t ivclnew,
  int ipvc,VCStr& vcs)
  {

//C Remove particle from the old verlet cell structure
    if(vcs.incb[ivcl]==0){rmip1(ivcl,ipvc,vcs);}       // outside VSCB
    else                                               // inside VSCB
    {
      IJK ivclijk=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(ivcl));
      size_t ilst=grid_ijk_to_index(vcs.mvclcb,ivclijk-vcs.mvclg);
      assert(ilst>=0&&ilst<vcs.nvclcb);
      if(vcs.incb[ivcl]>1){vcs.icllst[ilst].second--;rmip1(ivcl,ipvc,vcs);} // ESCB atom in inner part of VSCB
      else if(vcs.istatcl[icl]!=4){rmip1(ivcl,ipvc,vcs);}                   // ghost atom
      else{vcs.icllst[ilst].second--;rmip2(ilst,ivcl,ipvc,vcs);}            // ESCB atom at edge of VSCB
    }
//    printf("particle removed from old cell; \n");
      
//C Add particle to the new verlet cell structure
    if(vcs.incb[ivclnew]==0){addip1(icl,ipc,ivclnew,vcs);}  // outside VSCB
    else                                                    // inside VSCB
    {
      IJK ivclijk=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(ivclnew));
//      printf("ivclnew,ivclijk :%6zu%6zd%6zd%6zd\n",ivclnew,ivclijk.i,ivclijk.j,ivclijk.k);
      size_t ilst=grid_ijk_to_index(vcs.mvclcb,ivclijk-vcs.mvclg);
      assert(ilst>=0&&ilst<vcs.nvclcb);
//      printf("ilst,incbnew,istatcl :%6zu%6d%6d\n",ilst,vcs.incb[ivclnew],vcs.istatcl[icl]);
      if(vcs.incb[ivclnew]>1){addip1(icl,ipc,ivclnew,vcs);vcs.icllst[ilst].second++;} // ESCB atom in inner part of VSCB
      else if(vcs.istatcl[icl]!=4){addip1(icl,ipc,ivclnew,vcs);}                      // ghost atom
      else{addip2(ilst,icl,ipc,ivclnew,vcs);vcs.icllst[ilst].second++;}               // edge of VSCB  
    }
//    printf("particle added to new cell; \n");

  } // end of routine vclmigr

//======================================================================//
     
  inline void rmip1(size_t ivcl,int ipvc,VCStr& vcs)
  {
    int ipvcm=vcs.pmapinv[ivcl].size()-1;
    if(ipvc<ipvcm)
    {
      vcs.atpos[ivcl][ipvc]=vcs.atpos[ivcl][ipvcm];
      vcs.pmapinv[ivcl][ipvc]=vcs.pmapinv[ivcl][ipvcm];
      getpmap(ivcl,ipvc,vcs.pmapinv[ivcl][ipvc],vcs.pmap);
    }  
    vcs.atpos[ivcl].pop_back();
    vcs.pmapinv[ivcl].pop_back();
  } // end of routine migrp

//======================================================================//
     
  inline void rmip2(size_t ilst,size_t ivcl,int ipvc,VCStr& vcs)
  {
    int ipvcm1=vcs.icllst[ilst].second;
    int ipvcm2=vcs.pmapinv[ivcl].size()-1;
    if(ipvc<ipvcm1)
    {
      vcs.atpos[ivcl][ipvc]=vcs.atpos[ivcl][ipvcm1];
      vcs.pmapinv[ivcl][ipvc]=vcs.pmapinv[ivcl][ipvcm1];
      getpmap(ivcl,ipvc,vcs.pmapinv[ivcl][ipvc],vcs.pmap);
      if(ipvcm1<ipvcm2)
      {
        vcs.atpos[ivcl][ipvcm1]=vcs.atpos[ivcl][ipvcm2];
        vcs.pmapinv[ivcl][ipvcm1]=vcs.pmapinv[ivcl][ipvcm2];
        getpmap(ivcl,ipvcm1,vcs.pmapinv[ivcl][ipvcm1],vcs.pmap);
      }  
    }
    else if(ipvc<ipvcm2)
    {
      vcs.atpos[ivcl][ipvc]=vcs.atpos[ivcl][ipvcm2];
      vcs.pmapinv[ivcl][ipvc]=vcs.pmapinv[ivcl][ipvcm2];
      getpmap(ivcl,ipvc,vcs.pmapinv[ivcl][ipvc],vcs.pmap);
    }  
    vcs.atpos[ivcl].pop_back();
    vcs.pmapinv[ivcl].pop_back();

//    getpmap(ivcl,ipvc,vcs.pmapinv[ivcl][ipvc],vcs.pmap);
//    vcs.pmapinv[ivcl][ipvcmcb]=vcs.pmapinv[ivcl][ipvcm];
//    getpmap(ivcl,ipvcmcb,vcs.pmapinv[ivcl][ipvcmcb],vcs.pmap);
//    vcs.pmapinv[ivcl].pop_back();
  } // end of routine migrp

//======================================================================//
     
  inline void addip1(size_t icl,size_t ipc,size_t ivclnew,VCStr& vcs)
  {
    vcs.atpos[ivclnew].push_back(vcs.atposcl[ipc]);
    vcs.pmapinv[ivclnew].push_back(pair<size_t,size_t>(icl,ipc));
    int ipvcnew=vcs.pmapinv[ivclnew].size()-1;
    vcs.pmap[icl][ipc]=pair<size_t,int>(ivclnew,ipvcnew);
  } // end of function addip1

//======================================================================//
     
  inline void addip2(size_t ilst,size_t icl,size_t ipc,size_t ivclnew,VCStr& vcs)
  {
    int ipvcm1=vcs.icllst[ilst].second;
    int npcli=vcs.pmapinv[ivclnew].size();

//    int npclit=vcs.atpos[ivclnew].size();                                 // FOR TESTING
//    printf("npcli,npclit,ipvcmcb :%6d%6d%6d\n",npcli,npclit,ipvcmcb);     // FOR TESTING

    if(npcli>ipvcm1)
    {
      vcs.atpos[ivclnew].push_back(vcs.atpos[ivclnew][ipvcm1]);
      vcs.atpos[ivclnew][ipvcm1]=vcs.atposcl[ipc];
      vcs.pmapinv[ivclnew].push_back(vcs.pmapinv[ivclnew][ipvcm1]);
      vcs.pmapinv[ivclnew][ipvcm1]=pair<size_t,size_t>(icl,ipc);
      vcs.pmap[icl][ipc]=pair<size_t,int>(ivclnew,ipvcm1);
      int ipvcnew=vcs.pmapinv[ivclnew].size()-1;
      getpmap(ivclnew,ipvcnew,vcs.pmapinv[ivclnew][ipvcnew],vcs.pmap);
    }
    else
    {
      vcs.atpos[ivclnew].push_back(vcs.atposcl[ipc]);
      vcs.pmapinv[ivclnew].push_back(pair<size_t,size_t>(icl,ipc));
      int ipvcnew=vcs.pmapinv[ivclnew].size()-1;
      getpmap(ivclnew,ipvcnew,vcs.pmapinv[ivclnew][ipvcnew],vcs.pmap);
    }

//    npclit=vcs.atpos[ivclnew].size();                                     // FOR TESTING
//    printf("npcli,npclit,ipvcmcb :%6d%6d%6d\n",npcli,npclit,ipvcmcb);     // FOR TESTING
//    printf("vcs.pmapinv[ivclnew][ipvcmcb] :%6zu%6zu\n",
//    vcs.pmapinv[ivclnew][ipvcmcb].first,vcs.pmapinv[ivclnew][ipvcmcb].second);     // FOR TESTING

  } // end of function addip2

//======================================================================//
     
  inline void getpmap(size_t ivcl,int ipvc,pair<size_t,size_t> pmapinv,
  vector<vector<pair<size_t,int>>>& pmap)
  {
//    size_t icl=pmapinv.first;
//    size_t ipc=pmapinv.second;
//    if(icl>vcs.ncl-1)
//    {
//      printf("BUG1; icl,ncl %6zu%6zu\n",icl,vcs.ncl);
//      abort();
//    }
//    size_t npcli=vcs.pmap[icl].size();
//    if(ipc>npcli-1)
//    {
//      printf("BUG2; ipc,npcli %6zu%6zu\n",ipc,npcli);
//      abort();
//    }
    pmap[pmapinv.first][pmapinv.second]=pair<size_t,int>(ivcl,ipvc);
//    size_t iclp=vcs.pmapinv[ivcl][ipvc].first;
//    size_t ipcp=vcs.pmapinv[ivcl][ipvc].second;
//    vcs.pmap[iclp][ipcp]=pair<size_t,int>(ivcl,ipvc);
  } // end of routine getpmap

//======================================================================//
     
  template<typename CellT>
  inline void chkpmap(CellT& cells,VCStr& vcs)
  {
    for(size_t icl=0;icl<vcs.ncl;icl++)
    {
      size_t npcli=cells[icl].size();
      for(size_t ipc=0;ipc<npcli;ipc++)
      {
        size_t ivcl=vcs.pmap[icl][ipc].first;
        int ipvc=vcs.pmap[icl][ipc].second;
        size_t iclt=vcs.pmapinv[ivcl][ipvc].first;
        size_t ipct=vcs.pmapinv[ivcl][ipvc].second;
        if(iclt!=icl)
        {  
          printf("BUG; icl,iclt %6zu%6zu\n",icl,iclt);
          printf("BUG; ipc,ipct %6zu%6zu\n",ipc,ipct);
          abort();
        }  
      }  
    }

//    size_t npt=0;
    for(size_t ivcl=0;ivcl<vcs.nvcl;ivcl++)
    {
      size_t npcli=vcs.pmapinv[ivcl].size();
//      npt+=npcli;
      for(int ipvc=0;ipvc<npcli;ipvc++)
      {
        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
        size_t ivclt=vcs.pmap[icl][ipc].first;
        int ipvct=vcs.pmap[icl][ipc].second;
        if(ivclt!=ivcl)
        {  
          printf("BUG; ivcl,ivclt %6zu%6zu\n",ivcl,ivclt);
          printf("BUG; ipvc,ipvct %6d%6d\n",ipvc,ipvct);
          abort();
        }  
      }  
    }
//    if(npt!=vcs.np)
//    {  
//      printf("BUG; np,npt %6zu%6zu\n",vcs.np,npt);
//      abort();
//    }  

    if(vcs.iprint>0)printf("pmap OK\n");

  } // end of routine chkpmap

//======================================================================//
     
  inline void chkistat(VCStr& vcs)
  {
    size_t ip=0;
    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++) // all central box verlet cells
    {
      int iwrong=0;
      size_t ivcl=vcs.icllst[ilst].first;
      int npcbcli=vcs.icllst[ilst].second;
      int npcli=vcs.pmapinv[ivcl].size();
      for(int ipvc=0;ipvc<npcbcli;ipvc++,ip++)
      {  
        if(vcs.istat[ip]!=4)
        {
          iwrong++;
          size_t icl=iclipc(ip,vcs).first;
          size_t ipc=iclipc(ip,vcs).second;
          printf("A) BUG istat wrong; icl,ipc,ivcl,ipvc,ip,istatip :%6zu%6zu%6zu%6d%6zu%6d\n",icl,ipc,ivcl,ipvc,ip,vcs.istat[ip]);
          printf("npcbcli,npcli :%6d%6d\n",npcbcli,npcli);
          abort();
        }
      }  
      for(int ipvc=npcbcli;ipvc<npcli;ipvc++,ip++)
      {  
        int istatip=vcs.istat[ip];
        if(abs(istatip)!=1)
        {
          iwrong++;
          size_t icl=iclipc(ip,vcs).first;
          size_t ipc=iclipc(ip,vcs).second;
          printf("B) BUG istat wrong; icl,ipc,ivcl,ipvc,ip,istatip :%6zu%6zu%6zu%6d%6zu%6d\n",icl,ipc,ivcl,ipvc,ip,istatip);
          printf("npcbcli,npcli :%6d%6d\n",npcbcli,npcli);
          abort();
        }
      }  
      if(iwrong>0)abort();
    }  

    if(vcs.iprint>0)printf("istat OK\n");

  } // end of routine chkistat

//======================================================================//
   
  template<typename CellT>
  inline void chkvcsupdate(CellT& cells,VCStr& vcs)
  {
    GridBlock cbox;
    cbox.start=vcs.mvclg;
    cbox.end=vcs.mvcl-vcs.mvclg;

    for(size_t ivcl=0;ivcl<vcs.nvcl;ivcl++) // all verlet cells
    {
//      int iwrong=0;
      size_t ip0=vcs.ip0vcl[ivcl];
      int npcli=vcs.pmapinv[ivcl].size();
//      IJK ijkcl=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(ivcl));
//      if(inside_block(cbox,ijkcl))
      if(vcs.incb[ivcl]>0)
      {
        IJK ijkcl=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(ivcl));
        ssize_t ilst=grid_ijk_to_index(vcs.mvclcb,ijkcl-vcs.mvclg);
//        if(ilst<0||ilst>=vcs.nvclcb)
//        {
//          printf("ivcl,ilst :%6zu%6zd\n",ivcl,ilst);
//          printf("ijkcl     :%8zd%8zd%8zd\n",ijkcl.i,ijkcl.j,ijkcl.k);
//          printf("mvcl      :%8zd%8zd%8zd\n",vcs.mvcl.i,vcs.mvcl.j,vcs.mvcl.k);
//          printf("mvclcb    :%8zd%8zd%8zd\n",vcs.mvclcb.i,vcs.mvclcb.j,vcs.mvclcb.k);
//          printf("mvclg     :%8zd%8zd%8zd\n",vcs.mvclg.i,vcs.mvclg.j,vcs.mvclg.k);
//        }
        assert(ilst>=0&&ilst<vcs.nvclcb);
        int npcbcli=vcs.icllst[ilst].second;
//        printf("icl,npcbcli,npcli :%6zu%6d%6d\n",icl,npcbcli,npcli);
        for(int ipvc=0;ipvc<npcbcli;ipvc++)
        {  
          size_t ip=ip0+ipvc;
          int istatdiff=4-vcs.istat[ip];
          if(abs(istatdiff)>0)
          {
            printf("01) BUG in vcsupdate; ivcl,ipvc :%6zu%6d\n",ivcl,ipvc);
            printf("ip,istat :%6zu%6d\n",ip,vcs.istat[ip]);
            abort();
          }
        }  
        for(int ipvc=npcbcli;ipvc<npcli;ipvc++)
        {  
          size_t ip=ip0+ipvc;
          int istatdiff=1-abs(vcs.istat[ip]);
          if(abs(istatdiff)>0)
          {
            printf("02) BUG in vcsupdate; ivcl,ipvc :%6zu%6d\n",ivcl,ipvc);
            printf("ip,istat :%6zu%6d\n",ip,vcs.istat[ip]);
//            printf("npcbcli,npcli :%6d%6d\n",npcbcli,npcli);
            abort();
          }
        }  
      }

      for(int ipvc=0;ipvc<npcli;ipvc++)
      {  
        size_t ip=ip0+ipvc;
        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
        GridBlock cbox;
        cbox.start={vcs.mvclg.i,vcs.mvclg.j,vcs.mvclg.k};
        cbox.end={vcs.mvcl.i-vcs.mvclg.i,vcs.mvcl.j-vcs.mvclg.j,vcs.mvcl.k-vcs.mvclg.k};
        IJK ijkcl=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(ivcl));
        bool incbox=inside_block(cbox,ijkcl);
        int istatdiff=100;
        if(incbox&&vcs.istat[ip]!=4)
        {
          if(vcs.istatcl[icl]==-2){istatdiff=vcs.istat[ip]-(vcs.istatcl[icl]+3);}
          else{istatdiff=vcs.istat[ip]-(vcs.istatcl[icl]+2);}
//          printf("Case 1\n");
//          printf("istatdiff,istat,istatcl :%6d%6d%6d\n",istatdiff,vcs.istat[ip],vcs.istatcl[icl]);
        }
        else
        {
          istatdiff=vcs.istat[ip]-vcs.istatcl[icl];
//          printf("Case 2\n");
//          printf("istatdiff,istat,istatcl :%6d%6d%6d\n",istatdiff,vcs.istat[ip],vcs.istatcl[icl]);
        }
        if(abs(istatdiff)>0)
        {
          printf("1) BUG in vcsupdate; ivcl,ipvc,icl,ipc :%6zu%6d%6zu%6zu\n",icl,ipc,ivcl,ipvc);
          printf("ip,istat,istatcl :%6zu%6d%6d\n",ip,vcs.istat[ip],vcs.istatcl[icl]);
          abort();
        }
        ssize_t iclz=floor(vcs.atpos[ivcl][ipvc].z*vcs.dvcli.z);
        ssize_t icly=floor(vcs.atpos[ivcl][ipvc].y*vcs.dvcli.y);
        ssize_t iclx=floor(vcs.atpos[ivcl][ipvc].x*vcs.dvcli.x);
        size_t ivclt=iclz*vcs.mvclxy+icly*vcs.mvcl.i+iclx;
        if(ivclt-ivcl>0)
        {
          printf("2) BUG in vcsupdate; ivcl,ivclt,ipvc,icl,ipc :%6zu%6zu%6d%6zu%6zu\n",ivcl,ivclt,ipvc,icl,ipc);
          abort();
        }

        Vec3d atposcl={cells[icl][field::rx][ipc],cells[icl][field::ry][ipc],cells[icl][field::rz][ipc]};
        Vec3d drip=vcs.atpos[ivcl][ipvc]-atposcl+vcs.xyz0;
        double dripsq=norm2(drip);
        if(dripsq>1e-10)
        {
          printf("3) BUG in vcsupdate; ivcl,ipvc,icl,ipc :%6zu%6d%6zu%6zu\n",ivcl,ipvc,icl,ipc);
          printf("atpos   :%16.10lf%16.10lf%16.10lf\n",vcs.atpos[ivcl][ipvc].x,vcs.atpos[ivcl][ipvc].y,vcs.atpos[ivcl][ipvc].z);
          printf("atposcl :%16.10lf%16.10lf%16.10lf\n",atposcl.x,atposcl.y,atposcl.z);
          abort();
        }
      }
//      if(iwrong>0)abort();
    }  
  } // end of routine chkvcsupdate

//======================================================================//
     
  template<typename CellT>
  inline void passEFV(CellT& cells,VCStr& vcs,NBStr& nbs)
  {

#ifdef EnIneV
    double cfe=1.0/1.03641882007443324881e-04;
#else
    double cfe=1.;
#endif    

    for(size_t ivcl=0;ivcl<vcs.nvcl;ivcl++)
    {
      size_t ip0=vcs.ip0vcl[ivcl];
      int npcli=vcs.pmapinv[ivcl].size();

      for(int ipvc=0;ipvc<npcli;ipvc++)
      {
        size_t ip=ip0+ipvc;
        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
        assert(icl>=0&&icl<vcs.ncl);
//#ifdef EnIneV
        cells[icl][field::ep][ipc]=cfe*(nbs.esrpp[ip]+nbs.elrpp[ip]);
        cells[icl][field::fx][ipc]=cfe*nbs.force[ip].x;
        cells[icl][field::fy][ipc]=cfe*nbs.force[ip].y;
        cells[icl][field::fz][ipc]=cfe*nbs.force[ip].z;
//        if( particle_viriels != nullptr )
//        {
//          particle_viriels[ipc]=cfe*viriel;
//        }
//
//#else
//
//        cells[icl][field::ep][ipc]=nbs.esrpp[ip]+nbs.elrpp[ip];
//        cells[icl][field::fx][ipc]=nbs.force[ip].x;
//        cells[icl][field::fy][ipc]=nbs.force[ip].y;
//        cells[icl][field::fz][ipc]=nbs.force[ip].z;
////        if( particle_viriels != nullptr )
////        {
//////          Mat3d viriel;
////          particle_viriels[ipc]=viriel;
////        }
//#endif        
      }  
    }

//C Loading virials
#ifdef VIRIEL
#ifdef VIRIELLOCAL
    for(size_t ivcl=0;ivcl<vcs.nvcl;ivcl++)
    {
     size_t ip0=vcs.ip0vcl[ivcl];
     int npcli=vcs.pmapinv[ivcl].size();
     for(int ipvc=0;ipvc<npcli;ipvc++)
     {
      size_t ip=ip0+ipvc;
      size_t icl=vcs.pmapinv[ivcl][ipvc].first;
      size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
      assert(icl>=0&&icl<vcs.ncl);
      Mat3d* particle_viriels=cells[icl].field_pointer_or_null(field::virial);
      particle_viriels[ipc]=-cfe*nbs.viriel[ip];
     }
    }
#else
    IJK iclijk{vcs.ngls0,vcs.ngls0,vcs.ngls0};
    size_t icl=grid_ijk_to_index(vcs.dims,iclijk);
    assert(icl>=0&&icl<vcs.ncl);
    Mat3d* particle_viriels=cells[icl].field_pointer_or_null(field::virial);
    particle_viriels[0]=-cfe*nbs.viriel[0];
#endif
#endif
//#ifdef VIRIEL
////    printVIR(0,nbs.viriel);
//    IJK iclijk{vcs.ngls0,vcs.ngls0,vcs.ngls0};
//    size_t icl=grid_ijk_to_index(vcs.dims,iclijk);
//    assert(icl>=0&&icl<vcs.ncl);
//    Mat3d* particle_viriels=cells[icl].field_pointer_or_null(field::virial);
//std::cout << icl << ' ' << -cfe*nbs.viriel[0] << std::endl;
//    particle_viriels[0]=-cfe*nbs.viriel[0];
////    double Ppot=nbs.viriel[0].m11*nbs.viriel[0].m22*nbs.viriel[0].m33/3.0;
////    Vec3d hcb={vcs.dimscb.i,vcs.dimscb.j,vcs.dimscb.k};
////    hcb=vcs.decl*hcb;
////    printf("hcb       :%14.8lf%14.8lf%14.8lf\n",hcb.x,hcb.y,hcb.z);
////    double vol=hcb.x*hcb.y*hcb.z;
////    Ppot=Ppot/vol;
////    printf("InitEF; hmat :%20.12lf%20.12lf%20.12lf\n",vcs.hmat.m11,vcs.hmat.m12,vcs.hmat.m13);
////    printf("InitEF; hmat :%20.12lf%20.12lf%20.12lf\n",vcs.hmat.m21,vcs.hmat.m22,vcs.hmat.m23);
////    printf("InitEF; hmat :%20.12lf%20.12lf%20.12lf\n",vcs.hmat.m31,vcs.hmat.m32,vcs.hmat.m33);
////    printf("Ppot =%20.8lf\n",Ppot);
////    nbs.viriel[0]=1.03641882007443324881e-04*nbs.viriel[0];
////    printf("getEF; viriel :%20.8lf%20.8lf%20.8lf\n",nbs.viriel[0].m11,nbs.viriel[0].m12,nbs.viriel[0].m13);
////    printf("getEF; viriel :%20.8lf%20.8lf%20.8lf\n",nbs.viriel[0].m21,nbs.viriel[0].m22,nbs.viriel[0].m23);
////    printf("getEF; viriel :%20.8lf%20.8lf%20.8lf\n",nbs.viriel[0].m31,nbs.viriel[0].m32,nbs.viriel[0].m33);
////    cin.get();
////    size_t npcli=cells[icl].size();
////    for(size_t ipc=1;ipc<npcli;ipc++)
////    {
////      if(abs(particle_viriels[ipc].m11)>1e-8)
////      {
////        printf("BUG, viriel non-zero; ipc,vir :%6zu%20.8lf\n",ipc,particle_viriels[ipc].m11);
////        abort();
////      }
////    }
//#endif    

  } // end of routine passEFV

//======================================================================//
// Finding the sr and first mr candidate neighbors the atoms in the central box

  void getnbcbp1(VLCHBOP& par,VCStr& vcs,NBStr& nbs,GhostTbls& gtbls)
  {
    size_t nnbfull0=0;
    size_t nnbfrac0=0;
    size_t nnbmrcn0=0;

    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
    {
      size_t icl=vcs.icllst[ilst].first;
      int npcli=vcs.pmapinv[icl].size();
      if(npcli==0)continue; 
      size_t ip0=vcs.ip0vcl[icl];
      for(int ipc=0;ipc<npcli;ipc++) 
      {
        size_t ip=ip0+ipc;
        int ispci=vcs.ispc[ip];
        int istatip=vcs.istat[ip];
        const Vec3d ri=vcs.atpos[icl][ipc];
        nbs.nnbfull[ip]=nnbfull0;
        nbs.nnbfrac[ip]=nnbfrac0;
        nbs.nnbmrcn[ip]=nnbmrcn0;
        nbs.nnbmrcn0[ip]=nnbmrcn0;

//C intracell neighbors first...
        if(npcli>1){getcbnb(par,ip,ispci,istatip,icl,ipc,ip0,ri,vcs,nbs);}

        for(int inbc=1;inbc<vcs.nnbc[1];inbc++)   // first layer neighbors
        {
          size_t jcl=icl+vcs.ijcl[inbc];
          assert(jcl>=0&&jcl<icl);
          int npclj=vcs.pmapinv[jcl].size();
          if(npclj==0)continue;
          size_t jp0=vcs.ip0vcl[jcl];
          getcbnb(par,ip,ispci,istatip,jcl,npclj,jp0,ri,vcs,nbs);
        }    // end loop inbc
  
        if(vcs.incb[icl]==1)    // cell icl at edge of central box
        {
          for(int inbc=vcs.nnbc[1]-1;inbc>0;inbc--)    // first layer neighbors
          {
            size_t jcl=icl-vcs.ijcl[inbc];
            assert(jcl>=0&&jcl<vcs.nvcl);
            if(vcs.incb[jcl]>0)continue;
            int npclj=vcs.pmapinv[jcl].size();
            if(npclj==0)continue;
            size_t jp0=vcs.ip0vcl[jcl];
            getcbnb(par,ip,ispci,istatip,jcl,npclj,jp0,ri,vcs,nbs);
          }    // end loop inbc
        }    // end if
        nnbfull0=nbs.nnbfull[ip];
        nnbfrac0=nbs.nnbfrac[ip];
        nnbmrcn0=nbs.nnbmrcn[ip];
      }    // end loop ipc
    }  // end loop ilst

// At this point the number of full and fractional neighbors of all central
// atoms are known. In the next we can now construct definitive neighbor tables.

    size_t npcb=vcs.npcb;
    size_t nnbtot=nbs.nnb0[npcb-1];
    size_t nnbipm1=nbs.nnb0[0];
    nbs.nnbsrm=nnbipm1;
    nbs.nnb0[0]=0;
    for(size_t ip=1;ip<npcb;ip++)
    {
      size_t nnbip=nbs.nnb0[ip];
      if(nnbip>nbs.nnbsrm){nbs.nnbsrm=nnbip;}
      nbs.nnb0[ip]=nbs.nnb0[ip-1]+nnbipm1;
      nnbipm1=nnbip;
    }
    nnbtot+=nbs.nnb0[npcb-1];
 
//C Add full SR neighbors of all central atoms to final neighbor list
    Vec3d vij,sigij;
    double sn[3]={0.,1.,0.};
    size_t ip=0,inb=0,inbji;
    size_t ngbsrfl2=0,ngbsrfr2=0;
    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
    {
      size_t icl=vcs.icllst[ilst].first;
      int npcli=vcs.pmapinv[icl].size();
      for(int ipc=0;ipc<npcli;ipc++,ip++)
      {
        int ispci=vcs.ispc[ip];
        int istatip=vcs.istat[ip];
        assert(istatip==-1||istatip==1||istatip==4);
        for(;inb<nbs.nnbfull[ip];inb++)
        {
          size_t jp=nbs.srfl[inb];
          int ispcj=vcs.ispc[jp];
          double rij=nbs.rijsrfl[inb]; 
          double rijinv=1.0/rij; 
          Vec3d sigij=rijinv*nbs.vijsrfl[inb]; 
          size_t inbij=nbs.nnb0[ip]++;
          if(jp<npcb)    // jcl,jp in central box
          {
            getinbji(jp,inbij,inbji,nbs);
            add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
            add2srtbl1(jp,ip,ispci,inbji,rij,-sigij,sn,nbs);
          }  
          else    // jcl,jp outside central box
          {
            add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
            int istatjp=vcs.istat[jp];
            assert(abs(istatjp)<4&&istatjp!=0);
            if(abs(istatjp)>1)
            {
              assert(gtbls.nglst1<gtbls.nglst1m-1);
              assert(gtbls.nglst2<gtbls.nglst2m-1);
              if(istatip==4)
              {
                if(abs(istatjp)==2){vcs.istat[jp]=1;}else{vcs.istat[jp]=-1;}
                assert(gtbls.nglst1<gtbls.nglst1m-1);
                assert(gtbls.ngbsrfl1<gtbls.nglst1m-1);
                gtbls.glst1[gtbls.nglst1++]=jp;
                gtbls.gbsrfl1[gtbls.ngbsrfl1++]=pair<size_t,size_t>(inbij,ip);
              }
              else
              {
                assert(gtbls.nglst2<gtbls.nglst2m-1);
                if(istatjp<0){gtbls.glst2[gtbls.nglst2++]=jp;vcs.istat[jp]=-istatjp;}
                assert(ngbsrfl2<gtbls.nglst2m-1);
                gtbls.gbsrfl2[ngbsrfl2++]=pair<size_t,size_t>(inbij,ip);
              }  
            }
            else  
            {
              assert(gtbls.ngbsrfl1<gtbls.nglst1m-1);
              gtbls.gbsrfl1[gtbls.ngbsrfl1++]=pair<size_t,size_t>(inbij,ip);
            }
          }  
        }  
      }  
    }  
//    nbs.srfl.clear();
    vector<size_t>::iterator itint=nbs.nnb0.begin();
    nbs.nnbfull.assign(itint,itint+npcb);

//C Add fractional SR neighbors of all central atoms to final neighbor list
    inb=0;
    ip=0;
    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
    {
      size_t icl=vcs.icllst[ilst].first;
      int npcli=vcs.pmapinv[icl].size();
      for(int ipc=0;ipc<npcli;ipc++,ip++)
      {
        int ispci=vcs.ispc[ip];
        int istatip=vcs.istat[ip];
        assert(istatip==-1||istatip==1||istatip==4);
        for(;inb<nbs.nnbfrac[ip];inb++)
        {
          size_t jp=nbs.srfr[inb];
          int ispcj=vcs.ispc[jp];
          double rij=nbs.rijsrfr[inb]; 
          double rijinv=1.0/rij; 
          Vec3d sigij=rijinv*nbs.vijsrfr[inb]; 
          getsn(par,ispci,ispcj,rij,sn);
          size_t inbij=nbs.nnb0[ip]++;
          if(jp<npcb)
          {
            getinbji(jp,inbij,inbji,nbs);
            add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
            add2srtbl1(jp,ip,ispci,inbji,rij,-sigij,sn,nbs);
          }  
          else
          {
            add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
            int istatjp=vcs.istat[jp];
            assert(abs(istatjp)<4&&istatjp!=0);
            if(abs(istatjp)>1)
            {
              if(istatip==4)
              {
                if(abs(istatjp)==2){vcs.istat[jp]=1;}else{vcs.istat[jp]=-1;}
                assert(gtbls.nglst1<gtbls.nglst1m-1);
                assert(gtbls.ngbsrfr1<gtbls.nglst1m-1);
                gtbls.glst1[gtbls.nglst1++]=jp;
                gtbls.gbsrfr1[gtbls.ngbsrfr1++]=pair<size_t,size_t>(inbij,ip);
              }
              else
              {
                assert(gtbls.nglst2<gtbls.nglst2m-1);
                if(istatjp<0){gtbls.glst2[gtbls.nglst2++]=jp;vcs.istat[jp]=-istatjp;}
                assert(ngbsrfr2<gtbls.nglst2m-1);
                gtbls.gbsrfr2[ngbsrfr2++]=pair<size_t,size_t>(inbij,ip);
              }
            }
            else  
            {
              assert(gtbls.ngbsrfr1<gtbls.nglst1m-1);
              gtbls.gbsrfr1[gtbls.ngbsrfr1++]=pair<size_t,size_t>(inbij,ip);
            }
          }  
        }  
      }  
    }  

//C Add candidate MR active ghost neighbors of all ESCB atoms to 
//C ghost list 1 (glst1)
    ip=0;
    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
    {
      size_t icl=vcs.icllst[ilst].first;
      int npcli=vcs.pmapinv[icl].size();
      for(int ipc=0;ipc<npcli;ipc++,ip++)
      {
        if(vcs.istat[ip]<0)continue;
        for(size_t inb=nbs.nnbmrcn0[ip];inb<nbs.nnbmrcn[ip];inb++)
        {
          size_t jp=nbs.mrcn[inb];
          if(jp<npcb)continue;
          int istatjp=vcs.istat[jp];
          assert(istatjp>-4&&istatjp<5);
          if(abs(istatjp)==2)
          {
            assert(gtbls.nglst1<gtbls.nglst1m-1);
            gtbls.glst1[gtbls.nglst1++]=jp;vcs.istat[jp]=0;
          }
        }
      }
    }

    for(size_t ilst=0;ilst<ngbsrfl2;ilst++)
    {
      size_t inbij=gtbls.gbsrfl2[ilst].first;
      size_t jp=nbs.nbtbl[inbij];
      if(abs(vcs.istat[jp])<2){gtbls.gbsrfl1[gtbls.ngbsrfl1++]=gtbls.gbsrfl2[ilst];}
      else{gtbls.gbsrfl2[gtbls.ngbsrfl2++]=gtbls.gbsrfl2[ilst];}
    }

    for(size_t ilst=0;ilst<ngbsrfr2;ilst++)
    {
      size_t inbij=gtbls.gbsrfr2[ilst].first;
      size_t jp=nbs.nbtbl[inbij];
      if(abs(vcs.istat[jp])<2){gtbls.gbsrfr1[gtbls.ngbsrfr1++]=gtbls.gbsrfr2[ilst];}
      else{gtbls.gbsrfr2[gtbls.ngbsrfr2++]=gtbls.gbsrfr2[ilst];}
    }

    itint=nbs.nnb0.begin();
    nbs.nnbfrac.assign(itint,itint+npcb);
    nbs.inbm=nbs.nnbfrac[npcb-1];
    nbs.inbmrm=nbs.nnbmrcn[npcb-1];

    nbs.nnb0[0]=0;
    for(size_t ip=1;ip<npcb;ip++){nbs.nnb0[ip]=nbs.nnbfrac[ip-1];}

//    printf("End of  getnbcbp1:%8zu\n",nnbtot);
//    cin.get();
  }  //  end getnbcbp1

//======================================================================//
// Getting the second mr candidate neighbors the atoms in the central box

  void getnbcbp2(VLCHBOP& par,VCStr& vcs,NBStr& nbs,GhostTbls& gtbls)
  {
    gtbls.nglst2=0;

    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++) // all central box verlet cells
    {
      size_t icl=vcs.icllst[ilst].first;
      int npcli=vcs.pmapinv[icl].size();
      if(npcli==0)continue; 
      size_t ip0=vcs.ip0vcl[icl];

      for(int ipc=0;ipc<npcli;ipc++) // all atoms in cell ilst
      {
        size_t ip=ip0+ipc;
        int ispci=vcs.ispc[ip];
        int istatip=vcs.istat[ip];
        const Vec3d ri=vcs.atpos[icl][ipc];
        double ndbij=nbs.ndb[ip];

        for(int inbc=vcs.nnbc[1];inbc<vcs.nnbc[2];inbc++)   // second layer neighbors
        {
          size_t jcl=icl+vcs.ijcl[inbc];
          assert(jcl>=0&&jcl<vcs.nvcl);
          int npclj=vcs.pmapinv[jcl].size();
          if(npclj==0)continue;
          size_t jp0=vcs.ip0vcl[jcl];
          for(int jpc=0;jpc<npclj;jpc++) 
          {
            size_t jp=jp0+jpc;
            int ispcj=vcs.ispc[jp];
            Vec3d vij=vcs.atpos[jcl][jpc]-ri;
            vij=vcs.hmat*vij;
            double rijsq=norm2(vij);
            if(rijsq<par.pVLR[ispci][ispcj].r2vlrsq)
            {
              double rij=0.;
              int istatjp=vcs.istat[jp];
              if(activebond(istatip,istatjp))
              {
                rij=sqrt(rijsq);
                double vlr,dvlr,ddvlr;
                par.pVLR[ispci][ispcj].vlr(rij,vlr,dvlr,ddvlr);
                double pfdfrc=dvlr/rij;
                Vec3d dfrc=pfdfrc*vij;
                addvlr(ip,jp,vlr,vij,dfrc,nbs);
                if(rijsq>=par.pRC[ispci][ispcj].rcmrmaxsq){continue;}
              }
              else if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq){rij=sqrt(rijsq);}
              else{continue;}

              int imrnb;
              if(istatjp==4)      // SR tables of jp are available
              {
                getmrnbij(par,ip,jp,ispci,ispcj,1,imrnb,ndbij,rij,vij,vcs,nbs);
              }
              else if(istatip>0)  // otherwise MR bond is inactive,not needed
              {
                assert(abs(istatjp)<4);
                if(istatjp==1)
                {
                  getmrnbij(par,ip,jp,ispci,ispcj,1,imrnb,ndbij,rij,vij,vcs,nbs);
                }
                else
                {
                  assert(abs(istatjp)<4||istatjp!=1);
                  if(abs(istatjp)==2||istatjp==0)
                  {
                    if(istatjp!=0)
                    {
                      if(istatjp==-2){getnbsgp(par,jcl,jpc,npclj,jp,jp0,vcs,nbs);}
                      vcs.istat[jp]=0;
                      for(size_t inbj=nbs.nnb0[jp];inbj<nbs.nnbfrac[jp];inbj++)
                      {
                        size_t kp=nbs.nbtbl[inbj];
                        size_t kcl=vcs.ivclp[kp];
                        assert(kcl>=0&&kcl<vcs.nvcl);
                        int istatkp=vcs.istat[kp];
                        if(istatkp<-1)
                        {
                          int npclk=vcs.pmapinv[kcl].size();
                          size_t kp0=vcs.ip0vcl[kcl];
                          int kpc=kp-kp0;
                          getnbsgp(par,kcl,kpc,npclk,kp,kp0,vcs,nbs);
                          vcs.istat[kp]=-istatkp;
                        }  
                      }
                    }
                    getndbc2(par,jp,ispcj,vcs.ispc,nbs);
                    getmrnbij(par,ip,jp,ispci,ispcj,1,imrnb,ndbij,rij,vij,vcs,nbs);
                    if(imrnb>0)
                    {
                      vcs.istat[jp]=1;
                      assert(nbs.nnbfrac[jp]>=nbs.nnb0[jp]);
                      size_t nnbsrjp=nbs.nnbfrac[jp]-nbs.nnb0[jp];
                      if(nnbsrjp>nbs.nnbsrm){nbs.nnbsrm=nnbsrjp;}
                      double ndbjk=nbs.ndb[jp];
                      getmrnb0(par,jp,ispcj,1,ndbjk,vcs,nbs);
                      getmrnb1(par,jp,ispcj,1,ndbjk,vcs,nbs);
                      getmrnb2(par,jp,ispcj,ndbjk,vcs,nbs,gtbls);
                      assert(gtbls.nglst2<gtbls.nglst2m-1);
                      gtbls.glst2[gtbls.nglst2++]=jp;
                    }
                  }  
                  else
                  {
                    vcs.istat[jp]=-1;
                    getmrnbij(par,ip,jp,ispci,ispcj,0,imrnb,ndbij,rij,vij,vcs,nbs);
                  }
                }
              }
            }
          }   // end loops jpc
        }   // end loop inbc

        if(istatip>0&&vcs.incb[icl]<3)  // cell icl within two layers from the edge of central box
        {
          for(int inbc=vcs.nnbc[2]-1;inbc>=vcs.nnbc[1];inbc--)   // second layer neighbors
          {
            size_t jcl=icl-vcs.ijcl[inbc];
            assert(jcl>=0&&jcl<vcs.nvcl);
            int npclj=vcs.pmapinv[jcl].size();
            if(npclj==0)continue;
            if(vcs.incb[jcl]>0)continue;      // only neighbors cells outside central box
            size_t jp0=vcs.ip0vcl[jcl];
            for(int jpc=0;jpc<npclj;jpc++) 
            {
              size_t jp=jp0+jpc;
              int istatjp=vcs.istat[jp];
              int ispcj=vcs.ispc[jp];
              Vec3d vij=vcs.atpos[jcl][jpc]-ri;
              vij=vcs.hmat*vij;
              double rijsq=norm2(vij);
              if(rijsq<par.pVLR[ispci][ispcj].r2vlrsq)
              {
                double rij=0.;
                if(istatip==4&&istatjp!=-1&&abs(istatjp)!=3)
                {
                  rij=sqrt(rijsq);
                  double vlr,dvlr,ddvlr;
                  par.pVLR[ispci][ispcj].vlr(rij,vlr,dvlr,ddvlr);
                  double pfdfrc=dvlr/rij;
                  Vec3d dfrc=pfdfrc*vij;
                  addvlr(ip,jp,vlr,vij,dfrc,nbs);
                  if(rijsq>=par.pRC[ispci][ispcj].rcmrmaxsq){continue;}
                }
                else if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq){rij=sqrt(rijsq);}
                else{continue;}

                int imrnb;
                assert(abs(istatjp)<4);
                if(istatjp==1)
                {
                  getmrnbij(par,ip,jp,ispci,ispcj,1,imrnb,ndbij,rij,vij,vcs,nbs);
                }
                else
                {
                  if(abs(istatjp)==2||istatjp==0)
                  {
                    if(istatjp!=0)
                    {
                      if(istatjp==-2){getnbsgp(par,jcl,jpc,npclj,jp,jp0,vcs,nbs);}
                      vcs.istat[jp]=0;
                      for(size_t inbj=nbs.nnb0[jp];inbj<nbs.nnbfrac[jp];inbj++)
                      {
                        size_t kp=nbs.nbtbl[inbj];
                        size_t kcl=vcs.ivclp[kp];
                        assert(kcl>=0&&kcl<vcs.nvcl);
                        int istatkp=vcs.istat[kp];
                        if(istatkp<-1)
                        {
                          int npclk=vcs.pmapinv[kcl].size();
                          size_t kp0=vcs.ip0vcl[kcl];
                          int kpc=kp-kp0;
                          getnbsgp(par,kcl,kpc,npclk,kp,kp0,vcs,nbs);
                          vcs.istat[kp]=-istatkp;
                        }  
                      }  
                    }
                    getndbc2(par,jp,ispcj,vcs.ispc,nbs);
                    getmrnbij(par,ip,jp,ispci,ispcj,1,imrnb,ndbij,rij,vij,vcs,nbs);
                    if(imrnb>0)
                    {
                      vcs.istat[jp]=1;
                      assert(nbs.nnbfrac[jp]>=nbs.nnb0[jp]);
                      size_t nnbsrjp=nbs.nnbfrac[jp]-nbs.nnb0[jp];
                      if(nnbsrjp>nbs.nnbsrm){nbs.nnbsrm=nnbsrjp;}
                      double ndbjk=nbs.ndb[jp];
                      getmrnb0(par,jp,ispcj,1,ndbjk,vcs,nbs);
                      getmrnb1(par,jp,ispcj,1,ndbjk,vcs,nbs);
                      getmrnb2(par,jp,ispcj,ndbjk,vcs,nbs,gtbls);
                      assert(gtbls.nglst2<gtbls.nglst2m-1);
                      gtbls.glst2[gtbls.nglst2++]=jp;
                    }
                  }
                  else
                  {
                    vcs.istat[jp]=-1;
                    getmrnbij(par,ip,jp,ispci,ispcj,0,imrnb,ndbij,rij,vij,vcs,nbs);
                  }
                }  
              }   // end if(rij<rvlr)
            }   // end loops jpc
          }   // end loop inbc
        }   // end if condition
      }   // end loops ipc
    }   // end loops ilst
  }   // end function getnbcbp2

//======================================================================//
// Get third verlet layer LR neighbors of all atoms in the central box

  void getnbcbp3(VLCHBOP& par,VCStr& vcs,NBStr& nbs)
  {
    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
    {
      size_t icl=vcs.icllst[ilst].first;
      int npcli=vcs.pmapinv[icl].size();
      if(npcli==0)continue; 
      size_t ip0=vcs.ip0vcl[icl];
      for(int ipc=0;ipc<npcli;ipc++) 
      {
        size_t ip=ip0+ipc;
        int ispci=vcs.ispc[ip];
        int istatip=vcs.istat[ip];
        Vec3d ri=vcs.atpos[icl][ipc];
        for(int inbc=vcs.nnbc[2];inbc<vcs.nnbc[3];inbc++)   // third layer neighbors
        {
          size_t jcl=icl+vcs.ijcl[inbc];
          assert(jcl>=0&&jcl<vcs.nvcl);
          int npclj=vcs.pmapinv[jcl].size();
          if(npclj==0)continue;
          size_t jp0=vcs.ip0vcl[jcl];
          for(int jpc=0;jpc<npclj;jpc++) 
          {
            size_t jp=jp0+jpc;
            int istatjp=vcs.istat[jp];
            if(activebond(istatip,istatjp)==false)continue;
            int ispcj=vcs.ispc[jp];
            Vec3d vij=vcs.atpos[jcl][jpc]-ri;
            vij=vcs.hmat*vij;
            double rijsq=norm2(vij);
            if(rijsq<par.pVLR[ispci][ispcj].r2vlrsq)
            {
              double rij=sqrt(rijsq),vlr,dvlr,ddvlr;
              par.pVLR[ispci][ispcj].vlr(rij,vlr,dvlr,ddvlr);
              double pfdfrc=dvlr/rij;
              Vec3d dfrc=pfdfrc*vij;
              addvlr(ip,jp,vlr,vij,dfrc,nbs);
            }
          }   // end loops jpc
        }   // end loop inbc

        if(istatip==4&&vcs.incb[icl]!=4)  // cell icl within three layers from the edge of central box
        {
          for(int inbc=vcs.nnbc[3]-1;inbc>=vcs.nnbc[2];inbc--)   // third layer neighbors
          {
            size_t jcl=icl-vcs.ijcl[inbc];
            assert(jcl>=0&&jcl<vcs.nvcl);
            if(vcs.incb[jcl]>0)continue;          // jp should be outside central box
            int npclj=vcs.pmapinv[jcl].size();
            if(npclj==0)continue;
            size_t jp0=vcs.ip0vcl[jcl];
            for(int jpc=0;jpc<npclj;jpc++) 
            {
              size_t jp=jp0+jpc;
              int istatjp=vcs.istat[jp];
              if(istatjp==-1||abs(istatjp)==3)continue;      // jp is an in-active ghost
              int ispcj=vcs.ispc[jp];
              Vec3d vij=vcs.atpos[jcl][jpc]-ri;
              vij=vcs.hmat*vij;
              double rijsq=norm2(vij);
              if(rijsq<par.pVLR[ispci][ispcj].r2vlrsq)
              {
                double rij=sqrt(rijsq),vlr,dvlr,ddvlr;
                par.pVLR[ispci][ispcj].vlr(rij,vlr,dvlr,ddvlr);
                double pfdfrc=dvlr/rij;
                Vec3d dfrc=pfdfrc*vij;
                addvlr(ip,jp,vlr,vij,dfrc,nbs);
              }
            }   // end loops jpc
          }   // end loop inbc
        }   // end if condition
      }   // end loops ipc
    }   // end loops ilst
  }   // end function getnbcbp3

//======================================================================//

  inline void getnbsgp(VLCHBOP& par,size_t icl,int ipc,
  int npcli,size_t ip,size_t ip0,VCStr& vcs,NBStr& nbs)
  {
    nbs.nnb0[ip]=nbs.inbm;
    nbs.nnbfull[ip]=nbs.inbm;
    nbs.nnbfrac[ip]=0;
    nbs.nnbmrcn0[ip]=nbs.inbmrm;
    nbs.nnbmrcn[ip]=nbs.inbmrm;

    int ispci=vcs.ispc[ip];
    Vec3d ri=vcs.atpos[icl][ipc];
    double sn[3]={0.,1.,0.};

    if(npcli>1)  // intracell neighbors
    { 
      for(int jpc=0;jpc<npcli;jpc++)
      {
        if(jpc==ipc){continue;}
        size_t jp=ip0+jpc;
        getgnb3(par,icl,jpc,ip,jp,ispci,ri,sn,vcs,nbs);
      }
    }
     
    for(int inbc=1;inbc<vcs.nnbc[1];inbc++)   // cells jcl<icl
    {
      size_t jcl=icl+vcs.ijcl[inbc];
      assert(jcl>=0&&jcl<vcs.nvcl);
      int npclj=vcs.pmapinv[jcl].size();
      if(npclj==0)continue; 
      size_t jp0=vcs.ip0vcl[jcl];
      for(int jpc=0;jpc<npclj;jpc++) 
      {
        size_t jp=jp0+jpc;
        getgnb3(par,jcl,jpc,ip,jp,ispci,ri,sn,vcs,nbs);
      }
    }
   
    for(int inbc=vcs.nnbc[1]-1;inbc>0;inbc--)    // cells jcl>icl
    {
      size_t jcl=icl-vcs.ijcl[inbc];
      assert(jcl>=0&&jcl<vcs.nvcl);
      int npclj=vcs.pmapinv[jcl].size();
      if(npclj==0)continue; 
      size_t jp0=vcs.ip0vcl[jcl];
      for(int jpc=0;jpc<npclj;jpc++) 
      {
        size_t jp=jp0+jpc;
        getgnb3(par,jcl,jpc,ip,jp,ispci,ri,sn,vcs,nbs);
      }
    }

    double rij;
    Vec3d sigij;
    nbs.inbm=nbs.nnbfull[ip];
    for(size_t inb=0;inb<nbs.nnbfrac[ip];inb++)
    {
      size_t inbij=nbs.inbm++;
      size_t jp=nbs.srfr[inb];
      int ispcj=vcs.ispc[jp];
      getrsig(rij,sigij,nbs.vijsrfr[inb]);
      getsn(par,ispci,ispcj,rij,sn);
      add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
    }
    nbs.nnbfrac[ip]=nbs.inbm;
    nbs.inbmrm=nbs.nnbmrcn[ip];

//    printf("End of getnbsgp \n");
  }  // end getnbsgp

//======================================================================//

  inline void getnbgp1(VLCHBOP& par,VCStr& vcs,NBStr& nbs,GhostTbls& gtbls)
  {
    size_t nnbfull0=0,nnbfrac0=0;

    for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
    {
      size_t ip=gtbls.glst1[ilst];
      size_t icl=vcs.ivclp[ip];
      size_t ip0=vcs.ip0vcl[icl];
      int ipc=ip-ip0;
      int ispci=vcs.ispc[ip];
      int istatip=vcs.istat[ip];
      assert(icl>=0&&icl<vcs.nvcl);
      int npcli=vcs.pmapinv[icl].size();
      const Vec3d ri=vcs.atpos[icl][ipc];

      nbs.nnbfull[ip]=nnbfull0;
      nbs.nnbfrac[ip]=nnbfrac0;
      nbs.nnbmrcn[ip]=nbs.inbmrm;
      nbs.nnbmrcn0[ip]=nbs.inbmrm;

      if(npcli>1) //   intracell neighbors
      { 
        size_t ip0=ip-ipc;
        for(int jpc=0;jpc<npcli;jpc++)
        {
          if(jpc==ipc){continue;}
          size_t jp=ip0+jpc;
          getgnb1(par,icl,jpc,ip,jp,ispci,istatip,ri,vcs,nbs);
        }
      }
     
      for(int inbc=1;inbc<vcs.nnbc[1];inbc++)   // cells jcl<icl
      {
        size_t jcl=icl+vcs.ijcl[inbc];
        assert(jcl>=0&&jcl<vcs.nvcl);
        if(vcs.incb[jcl]>0){continue;}
        int npclj=vcs.pmapinv[jcl].size();
        if(npclj==0)continue;
        size_t jp0=vcs.ip0vcl[jcl];
        for(int jpc=0;jpc<npclj;jpc++) 
        {
          size_t jp=jp0+jpc;
          getgnb1(par,jcl,jpc,ip,jp,ispci,istatip,ri,vcs,nbs); 
        }
      }
   
      for(int inbc=vcs.nnbc[1]-1;inbc>0;inbc--)    // cells jcl>icl
      {
        size_t jcl=icl-vcs.ijcl[inbc];
        assert(jcl>=0&&jcl<vcs.nvcl);
        if(vcs.incb[jcl]>0){continue;}
        int npclj=vcs.pmapinv[jcl].size();
        if(npclj==0)continue;
        size_t jp0=vcs.ip0vcl[jcl];
        for(int jpc=0;jpc<npclj;jpc++) 
        {
          size_t jp=jp0+jpc;
          getgnb1(par,jcl,jpc,ip,jp,ispci,istatip,ri,vcs,nbs); 
        }
      }
      nnbfull0=nbs.nnbfull[ip];
      nnbfrac0=nbs.nnbfrac[ip];
      nbs.inbmrm=nbs.nnbmrcn[ip];
    }

    size_t ip0=gtbls.glst1[0];
    size_t nnbip0=nbs.nnb0[ip0];
    if(nnbip0>nbs.nnbsrm){nbs.nnbsrm=nnbip0;}
    nbs.nnb0[ip0]=nbs.inbm;
    for(size_t ilst=1;ilst<gtbls.nglst1;ilst++)
    {
      size_t ip=gtbls.glst1[ilst];
      size_t nnbip=nbs.nnb0[ip];
      if(nnbip>nbs.nnbsrm){nbs.nnbsrm=nnbip;}
      nbs.nnb0[ip]=nbs.nnb0[ip0]+nnbip0;
      ip0=ip;
      nnbip0=nnbip;
    }  

    size_t inb=0,inbji;
    double rij;
    Vec3d sigij;
    double sn[3]={0.,1.,0.};
    for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
    {
      size_t ip=gtbls.glst1[ilst];
      int ispci=vcs.ispc[ip];
      int istatip=vcs.istat[ip];
      for(;inb<nbs.nnbfull[ip];inb++)
      {
        size_t jp=nbs.srfl[inb];
        int ispcj=vcs.ispc[jp];
        int istatjp=vcs.istat[jp];
        getrsig(rij,sigij,nbs.vijsrfl[inb]);
        size_t inbij=nbs.nnb0[ip]++;
        assert(abs(istatjp)<4);
        if(abs(istatjp)>1)           // istatjp=-3,-2,2,3
        {
          add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
          gtbls.gbsrfl2[gtbls.ngbsrfl2++]=pair<size_t,size_t>(inbij,ip);
          if(istatjp<0)
          {
            if(istatip>=0)
            {
              assert(gtbls.nglst2<gtbls.nglst2m-1);
              gtbls.glst2[gtbls.nglst2++]=jp;
              vcs.istat[jp]=-istatjp;
            }
          }
        }
        else
        {
          getinbji(jp,inbij,inbji,nbs);
          add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
          add2srtbl1(jp,ip,ispci,inbji,rij,-sigij,sn,nbs);
        }
      }
    }

    for(size_t ilst=0;ilst<gtbls.ngbsrfl1;ilst++)
    {
      size_t inbij=gtbls.gbsrfl1[ilst].first;
      size_t ip=gtbls.gbsrfl1[ilst].second;
      size_t jp=nbs.nbtbl[inbij];
      int ispci=vcs.ispc[ip];
      add2srtbl2(ip,jp,ispci,inbij,nbs);
    }

    for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
    {
      size_t ip=gtbls.glst1[ilst];
      nbs.nnbfull[ip]=nbs.nnb0[ip];
    }  

    inb=0;
    for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
    {
      size_t ip=gtbls.glst1[ilst];
      int ispci=vcs.ispc[ip];
      int istatip=vcs.istat[ip];
      for(;inb<nbs.nnbfrac[ip];inb++)
      {
        size_t jp=nbs.srfr[inb];
        int ispcj=vcs.ispc[jp];
        int istatjp=vcs.istat[jp];
        getrsig(rij,sigij,nbs.vijsrfr[inb]);
        getsn(par,ispci,ispcj,rij,sn);
        size_t inbij=nbs.nnb0[ip]++;
        assert(abs(istatjp)<4);
        if(abs(istatjp)>1)  // istatjp=-3,-2,2,3
        {
          add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
          gtbls.gbsrfr2[gtbls.ngbsrfr2++]=pair<size_t,size_t>(inbij,ip);
          if(istatjp<0)
          {
            if(istatip>=0)
            {
              assert(gtbls.nglst2<gtbls.nglst2m-1);
              gtbls.glst2[gtbls.nglst2++]=jp;
              vcs.istat[jp]=-istatjp;
            }  
          }  
        }
        else
        {
          getinbji(jp,inbij,inbji,nbs);
          add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
          add2srtbl1(jp,ip,ispci,inbji,rij,-sigij,sn,nbs);
        }
      }
    }

    for(size_t ilst=0;ilst<gtbls.ngbsrfr1;ilst++)
    {
      size_t inbij=gtbls.gbsrfr1[ilst].first;
      size_t ip=gtbls.gbsrfr1[ilst].second;
      size_t jp=nbs.nbtbl[inbij];
      int ispci=vcs.ispc[ip];
      add2srtbl2(ip,jp,ispci,inbij,nbs);
    }

    ip0=gtbls.glst1[0];
    nbs.nnbfrac[ip0]=nbs.nnb0[ip0];
    nbs.nnb0[ip0]=nbs.inbm;
    for(size_t ilst=1;ilst<gtbls.nglst1;ilst++)
    {
      size_t ip=gtbls.glst1[ilst];
      nbs.nnbfrac[ip]=nbs.nnb0[ip];
      nbs.nnb0[ip]=nbs.nnbfrac[ip0];
      ip0=ip;
    }
    nbs.inbm=nbs.nnbfrac[ip0];

  } // end function getnbgp1

////======================================================================//
  
  inline void getnbgp2(VLCHBOP& par,VCStr& vcs,NBStr& nbs,GhostTbls& gtbls)
  {
    size_t nnbfull0=0;
    size_t nnbfrac0=0;
    size_t nglst2=0;
    for(size_t ilst=0;ilst<gtbls.nglst2;ilst++)
    {
      size_t ip=gtbls.glst2[ilst];
      int istatip=vcs.istat[ip];
      if(abs(istatip)<2)continue;
      gtbls.glst2[nglst2++]=ip;
      int ispci=vcs.ispc[ip];
      size_t icl=vcs.ivclp[ip];
      size_t ip0=vcs.ip0vcl[icl];
      int ipc=ip-ip0;
      assert(icl>=0&&icl<vcs.nvcl);
      int npcli=vcs.pmapinv[icl].size();
      const Vec3d ri=vcs.atpos[icl][ipc];

      nbs.nnbfull[ip]=nnbfull0;
      nbs.nnbfrac[ip]=nnbfrac0;
      nbs.nnbmrcn[ip]=nbs.inbmrm;
      nbs.nnbmrcn0[ip]=nbs.inbmrm;

      if(npcli>1) //   intracell neighbors
      { 
        size_t ip0=ip-ipc;
        for(int jpc=0;jpc<npcli;jpc++)
        {
          if(jpc==ipc){continue;}
          size_t jp=ip0+jpc;
          getgnb2(par,icl,jpc,ip,jp,ispci,istatip,ri,vcs,nbs);
        }
      }

      for(int inbc=1;inbc<vcs.nnbc[1];inbc++)    // cells jcl>icl
      {
        size_t jcl=icl+vcs.ijcl[inbc];
        assert(jcl>=0&&jcl<vcs.nvcl);
        if(vcs.incb[jcl]>0){continue;}
        int npclj=vcs.pmapinv[jcl].size();
        if(npclj==0)continue;
        size_t jp0=vcs.ip0vcl[jcl];
        for(int jpc=0;jpc<npclj;jpc++) 
        {
          size_t jp=jp0+jpc;
          getgnb2(par,jcl,jpc,ip,jp,ispci,istatip,ri,vcs,nbs);
        }
      }
  
      for(int inbc=vcs.nnbc[1]-1;inbc>0;inbc--)    // cells jcl>icl
      {
        size_t jcl=icl-vcs.ijcl[inbc];
        assert(jcl>=0&&jcl<vcs.nvcl);
        if(vcs.incb[jcl]>0){continue;}
        int npclj=vcs.pmapinv[jcl].size();
        if(npclj==0)continue;
        size_t jp0=vcs.ip0vcl[jcl];
        for(int jpc=0;jpc<npclj;jpc++) 
        {
          size_t jp=jp0+jpc;
          getgnb2(par,jcl,jpc,ip,jp,ispci,istatip,ri,vcs,nbs);
        }
      }
      nnbfull0=nbs.nnbfull[ip];
      nnbfrac0=nbs.nnbfrac[ip];
      nbs.inbmrm=nbs.nnbmrcn[ip];
    }
    gtbls.nglst2=nglst2;

    size_t ip0=gtbls.glst2[0];
    size_t nnbip0=nbs.nnb0[ip0];
    if(nnbip0>nbs.nnbsrm){nbs.nnbsrm=nnbip0;}
    nbs.nnb0[ip0]=nbs.inbm;
    for(size_t ilst=1;ilst<gtbls.nglst2;ilst++)
    {
      size_t ip=gtbls.glst2[ilst];
      size_t nnbip=nbs.nnb0[ip];
      if(nnbip>nbs.nnbsrm){nbs.nnbsrm=nnbip;}
      nbs.nnb0[ip]=nbs.nnb0[ip0]+nnbip0;
      ip0=ip;
      nnbip0=nnbip;
    }  

    size_t inb=0,inbji;
    double rij;
    Vec3d vij,sigij;
    double sn[3]={0.,1.,0.};
    for(size_t ilst=0;ilst<gtbls.nglst2;ilst++)
    {
      size_t ip=gtbls.glst2[ilst];
      int ispci=vcs.ispc[ip];
      for(;inb<nbs.nnbfull[ip];inb++)
      {
        size_t jp=nbs.srfl[inb];
        int ispcj=vcs.ispc[jp];
        int istatjp=vcs.istat[jp];
        getrsig(rij,sigij,nbs.vijsrfl[inb]);
        size_t inbij=nbs.nnb0[ip]++;
        assert(abs(istatjp)<4);
        if(istatjp<-1)  
        {
          add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
        }
        else if(istatjp>1)
        {
          getinbji(jp,inbij,inbji,nbs);
          add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
          add2srtbl1(jp,ip,ispci,inbji,rij,-sigij,sn,nbs);
        }
      }
    }

    for(size_t ilst=0;ilst<gtbls.ngbsrfl2;ilst++)
    {
      size_t inbij=gtbls.gbsrfl2[ilst].first;
      size_t jp=nbs.nbtbl[inbij];
      if(vcs.istat[jp]<-1)continue;
      size_t ip=gtbls.gbsrfl2[ilst].second;
      int ispci=vcs.ispc[ip];
      add2srtbl2(ip,jp,ispci,inbij,nbs);
    }

    for(size_t ilst=0;ilst<gtbls.nglst2;ilst++)
    {
      size_t ip=gtbls.glst2[ilst];
      nbs.nnbfull[ip]=nbs.nnb0[ip];
    }  

    inb=0;
    for(size_t ilst=0;ilst<gtbls.nglst2;ilst++)
    {
      size_t ip=gtbls.glst2[ilst];
      assert(ip<vcs.np);
      int ispci=vcs.ispc[ip];
      for(;inb<nbs.nnbfrac[ip];inb++)
      {
        size_t jp=nbs.srfr[inb];
        int ispcj=vcs.ispc[jp];
        vij=nbs.vijsrfr[inb];
        getrsig(rij,sigij,vij);
        getsn(par,ispci,ispcj,rij,sn);
        size_t inbij=nbs.nnb0[ip]++;
        int istatjp=vcs.istat[jp];
        assert(abs(istatjp)<4);
        if(istatjp<-1)  
        {
          add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
        }
        else if(istatjp>1)
        {
          getinbji(jp,inbij,inbji,nbs);
          add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
          add2srtbl1(jp,ip,ispci,inbji,rij,-sigij,sn,nbs);
        }
      }
    }

    for(size_t ilst=0;ilst<gtbls.ngbsrfr2;ilst++)
    {
      size_t inbij=gtbls.gbsrfr2[ilst].first;
      size_t jp=nbs.nbtbl[inbij];
      if(vcs.istat[jp]<-1)continue;
      size_t ip=gtbls.gbsrfr2[ilst].second;
      int ispci=vcs.ispc[ip];
      add2srtbl2(ip,jp,ispci,inbij,nbs);
    }

    ip0=gtbls.glst2[0];
    nbs.nnbfrac[ip0]=nbs.nnb0[ip0];
    nbs.nnb0[ip0]=nbs.inbm;
    for(size_t ilst=1;ilst<gtbls.nglst2;ilst++)
    {
      size_t ip=gtbls.glst2[ilst];
      nbs.nnbfrac[ip]=nbs.nnb0[ip];
      nbs.nnb0[ip]=nbs.nnbfrac[ip0];
      ip0=ip;
    }
    nbs.inbm=nbs.nnbfrac[ip0];

  }  // end getnbgp2

//======================================================================//
//
//  template<typename GridT>
//  inline void applypbc(ipbc,ri)
//  {
//    for(size_t idir=0;idir<3;idir++)
//    {
//      if(ipbc[idir]>0)
//      {
//        ri.x+=cbox.x;}
//    {  
//  }
//
//======================================================================//

  inline void getcbnb(VLCHBOP& par,size_t ip,int ispci,int istatip,size_t jcl,
  size_t jpcmax,size_t jp0,const Vec3d& ri,VCStr& vcs,NBStr& nbs)
  {
    for(size_t jpc=0;jpc<jpcmax;jpc++) 
    {
      size_t jp=jp0+jpc;
      int ispcj=vcs.ispc[jp];
      Vec3d vij=vcs.atpos[jcl][jpc]-ri;
      vij=vcs.hmat*vij;
      double rijsq=norm2(vij);
      if(rijsq<par.pRC[ispci][ispcj].r1srsq)
      {
        size_t inb=nbs.nnbfull[ip]++;
        assert(inb<nbs.nbndm);
        nbs.srfl[inb]=jp;
        nbs.rijsrfl[inb]=sqrt(rijsq);
        nbs.vijsrfl[inb]=vij;
        nbs.nnb0[ip]++;
        nbs.nnb0[jp]++;
      }  
      else if(rijsq<par.pVLR[ispci][ispcj].r2vlrsq)
      {
        double rij=0.;
        if(activebond(istatip,vcs.istat[jp]))
        {
          double vlr,dvlr,ddvlr;
          rij=sqrt(rijsq);
          par.pVLR[ispci][ispcj].vlr(rij,vlr,dvlr,ddvlr);
          double pfdfrc=dvlr/rij;
          Vec3d dfrc=pfdfrc*vij;
          addvlr(ip,jp,vlr,vij,dfrc,nbs);
          if(rijsq>=par.pRC[ispci][ispcj].rcmrmaxsq)continue;
        } 
        else if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq){rij=sqrt(rijsq);}
        else{continue;}

        if(rijsq<par.pRC[ispci][ispcj].r2srsq)
        {
          size_t inb=nbs.nnbfrac[ip]++;
          assert(inb<nbs.nbndm);
          nbs.srfr[inb]=jp;
          nbs.rijsrfr[inb]=rij;
          nbs.vijsrfr[inb]=vij;
          nbs.nnb0[ip]++;
          nbs.nnb0[jp]++;
        }
        else   // HIER: but if istatip=-1 and jp>npcb we don't need this MR candidate (?)
        {
          size_t inb=nbs.nnbmrcn[ip]++;
          assert(inb<nbs.nbndmrcnm);
          nbs.mrcn[inb]=jp;
          nbs.rijmrcn[inb]=rij;
          nbs.vijmrcn[inb]=vij;
        }
      }
    }
  }

//======================================================================//

  inline void getgnb1(VLCHBOP& par,size_t jcl,int jpc,size_t ip,size_t jp,
  int ispci,int istatip,const Vec3d ri,VCStr& vcs,NBStr& nbs)
  {
    assert(abs(istatip)<2);
    int istatjp=vcs.istat[jp];
    assert(abs(istatjp)<4);
    if(abs(istatjp)>1||jp<ip) // (jp not in glst1) or (jp in glst1 and jp<ip (avoiding double counting))
    {                                       // (jp in glst1 and jp<ip (avoiding double counting)) 
      int ispcj=vcs.ispc[jp];
      Vec3d vij=vcs.atpos[jcl][jpc]-ri;
      vij=vcs.hmat*vij;
      double rijsq=norm2(vij);
      if(rijsq<par.pRC[ispci][ispcj].r1srsq)
      {
        size_t inb=nbs.nnbfull[ip]++;
        assert(inb<nbs.nbndm);
        nbs.srfl[inb]=jp;
        nbs.vijsrfl[inb]=vij;
        nbs.nnb0[ip]++;
        nbs.nnb0[jp]++;
      }
      else if(rijsq<par.pRC[ispci][ispcj].r2srsq)
      {
        size_t inb=nbs.nnbfrac[ip]++;
        assert(inb<nbs.nbndm);
        nbs.srfr[inb]=jp;
        nbs.vijsrfr[inb]=vij;
        nbs.nnb0[ip]++;
        nbs.nnb0[jp]++;
      }
      else if(istatip>-1)
      {
        if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq)
        {
          size_t inb=nbs.nnbmrcn[ip]++;
          assert(inb<nbs.nbndmrcnm);
          nbs.mrcn[inb]=jp;
          nbs.rijmrcn[inb]=sqrt(rijsq);
          nbs.vijmrcn[inb]=vij;
        }
      }
    }  
    else if(istatip>-1)  // (jp in glst1 and jp>ip and ip active ghost; getting mrcn neighbors) 
    {
      int ispcj=vcs.ispc[jp];
      Vec3d vij=vcs.atpos[jcl][jpc]-ri;
      vij=vcs.hmat*vij;
      double rijsq=norm2(vij);
      if(rijsq<par.pRC[ispci][ispcj].r2srsq)return;
      else if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq)
      {
        size_t inb=nbs.nnbmrcn[ip]++;
        assert(inb<nbs.nbndmrcnm);
        nbs.mrcn[inb]=jp;
        nbs.rijmrcn[inb]=sqrt(rijsq);
        nbs.vijmrcn[inb]=vij;
      }
    }
  }

//======================================================================//

  inline void getgnb2(VLCHBOP& par,size_t jcl,int jpc,size_t ip,size_t jp,
  int ispci,int istatip,const Vec3d& ri,VCStr& vcs,NBStr& nbs)
  {
    int ispcj=vcs.ispc[jp];
    int istatjp=vcs.istat[jp];
    assert(abs(istatjp)<4);
    Vec3d vij=vcs.atpos[jcl][jpc]-ri;
    vij=vcs.hmat*vij;
    double rijsq=norm2(vij);
    if(istatjp<-1||(istatjp>1&&jp<ip)) // (jp not in glst1 or glst2) or
    {                                  // (jp in glst2 and jp<ip (avoiding double counting)) 
      if(rijsq<par.pRC[ispci][ispcj].r1srsq)
      {
        size_t inb=nbs.nnbfull[ip]++;
        assert(inb<nbs.nbndm);
        nbs.srfl[inb]=jp;
        nbs.vijsrfl[inb]=vij;
        nbs.nnb0[ip]++;
        nbs.nnb0[jp]++;
      }
      else if(rijsq<par.pRC[ispci][ispcj].r2srsq)
      {
        size_t inb=nbs.nnbfrac[ip]++;
        assert(inb<nbs.nbndm);
        nbs.srfr[inb]=jp;
        nbs.vijsrfr[inb]=vij;
        nbs.nnb0[ip]++;
        nbs.nnb0[jp]++;
      }
      else if(istatip==2)  // (ip is an active ghost; getting its mrcn neighbors) 
      {
        if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq)
        {
          size_t inb=nbs.nnbmrcn[ip]++;
          assert(inb<nbs.nbndmrcnm);
          nbs.mrcn[inb]=jp;
          nbs.rijmrcn[inb]=sqrt(rijsq);
          nbs.vijmrcn[inb]=vij;
        }
      }
    }
    else if(istatip==2)  // (ip is an active ghost; getting its mrcn neighbors) 
    {
      if(rijsq<par.pRC[ispci][ispcj].r2srsq)return;
      else if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq)
      {
        size_t inb=nbs.nnbmrcn[ip]++;
        assert(inb<nbs.nbndmrcnm);
        nbs.mrcn[inb]=jp;
        nbs.rijmrcn[inb]=sqrt(rijsq);
        nbs.vijmrcn[inb]=vij;
      }
    }
  }

//======================================================================//

  inline void getgnb3(VLCHBOP& par,size_t jcl,int jpc,size_t ip,size_t jp,
  int ispci,const Vec3d ri,double* sn,VCStr& vcs,NBStr& nbs)
  {
    int ispcj=vcs.ispc[jp];
    Vec3d vij=vcs.atpos[jcl][jpc]-ri;
    vij=vcs.hmat*vij;
    double rijsq=norm2(vij);
    if(rijsq<par.pRC[ispci][ispcj].r1srsq)
    {
      size_t inbij=nbs.nnbfull[ip]++;
      double rij;
      Vec3d sigij;
      getrsig(rij,sigij,vij);
      add2srtbl1(ip,jp,ispcj,inbij,rij,sigij,sn,nbs);
    }
    else if(rijsq<par.pRC[ispci][ispcj].r2srsq)
    {
      size_t inb=nbs.nnbfrac[ip]++;
      assert(inb<nbs.nbndm);
      nbs.srfr[inb]=jp;
      nbs.vijsrfr[inb]=vij;
    }
    else if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq)
    {
      size_t inb=nbs.nnbmrcn[ip]++;
      assert(inb<nbs.nbndmrcnm);
      nbs.mrcn[inb]=jp;
      nbs.rijmrcn[inb]=sqrt(norm2(vij));
      nbs.vijmrcn[inb]=vij;
    }
  }

//=====================================================================//

  inline void getrsig(double& rij,Vec3d& sigij,Vec3d& vij)
  {
    rij=sqrt(norm2(vij));
    double rijinv=1.0/rij;   //  ??????????  //
    sigij=rijinv*vij;    
//    sigij=(1.0/rij)*vij;   // Est-ce-que ce serait mieux ?????? //
//    sig=sig/rij;       
  }

//=====================================================================//

  inline void getsn(VLCHBOP& par,int ispci,int ispcj,double rij,double* sn)
  {
    double r1=par.pRC[ispci][ispcj].r1sr;
    double r2=par.pRC[ispci][ispcj].r2sr;
    double drinv=par.pRC[ispci][ispcj].r21inv;
    double p2=par.pRC[ispci][ispcj].p2z;
    double Sq,dSq;
    fcpdown(rij,r1,r2,drinv,p2,Sq,dSq);
    if(Sq<0)
    {
      cout<<"rij   ="<<rij<<endl;
      cout<<"r1    ="<<r1<<endl;
      cout<<"r2    ="<<r2<<endl;
      cout<<"drinv ="<<drinv<<endl;
      cout<<"p2    ="<<p2<<endl;
      cout<<"ERROR: Sq<0; Sq ="<<Sq<<endl;
      abort();
    }  
    sn[0]=Sq-1.;
    sn[1]=Sq;
    sn[2]=dSq;
  }  

//=====================================================================//

  inline void add2srtbl1(size_t ip,size_t jp,int ispcj,size_t inbij,
  double r,const Vec3d sig,double* sn,NBStr& nbs)
  {
    assert(inbij<nbs.nnnbm);
    nbs.nbtbl[inbij]=jp;
    nbs.rstor[inbij]=r;
    nbs.sigstor[inbij]=sig;
    nbs.snstor[inbij][0]=sn[0];
    nbs.snstor[inbij][1]=sn[1];
    nbs.snstor[inbij][2]=sn[2];
    nbs.crd[ip]+=sn[1];
    nbs.crdspc[ip][ispcj]+=sn[1];
  }  

//=====================================================================//

  inline void add2srtbl2(size_t ip,size_t jp,int ispci,size_t inbij,NBStr& nbs)
  {
    assert(inbij<nbs.nnnbm);
    size_t inbji;
    getinbji(jp,inbij,inbji,nbs);
    assert(inbji<nbs.nnnbm);
    double sn1=nbs.snstor[inbij][1];
    nbs.nbtbl[inbji]=ip;
    nbs.rstor[inbji]=nbs.rstor[inbij];
    nbs.sigstor[inbji]=-nbs.sigstor[inbij];
    nbs.snstor[inbji][0]=nbs.snstor[inbij][0];
    nbs.snstor[inbji][1]=sn1;
    nbs.snstor[inbji][2]=nbs.snstor[inbij][2];
    nbs.crd[jp]+=sn1;
    nbs.crdspc[jp][ispci]+=sn1;
  }  

//======================================================================//

  inline void getinbji(size_t jp,size_t inbij,size_t& inbji,NBStr& nbs)
  {
    inbji=nbs.nnb0[jp]++;
    nbs.indji[inbij]=inbji;
    nbs.indji[inbji]=inbij;
  }  

//======================================================================//
//C Getting first verlet layer MR neighbors of all atoms in the central box

  inline void getmrtbls1(VLCHBOP& par,VCStr& vcs,NBStr& nbs,GhostTbls& gtbls)
  {
//C Getting all MR neighbors of all central box atoms
    size_t ip=0;
    double ndbij;
    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++) // all central box verlet cells
    {
      size_t icl=vcs.icllst[ilst].first;
      int npcbcli=vcs.icllst[ilst].second;
      int npcli=vcs.pmapinv[icl].size();
      for(int ipc=0;ipc<npcbcli;ipc++,ip++)
      {  
        int ispci=vcs.ispc[ip];
        double ndbij=nbs.ndb[ip];
        getmrnb0(par,ip,ispci,4,ndbij,vcs,nbs);
        getmrnb1(par,ip,ispci,4,ndbij,vcs,nbs);
      }  
      for(int ipc=npcbcli;ipc<npcli;ipc++,ip++)
      {  
        int ispci=vcs.ispc[ip];
        int istatip=vcs.istat[ip];
        double ndbij=nbs.ndb[ip];
        getmrnb0(par,ip,ispci,istatip,ndbij,vcs,nbs);
        getmrnb1(par,ip,ispci,istatip,ndbij,vcs,nbs);
      }
    }  

//C Getting all MR neighbors of ghost atoms that are active first
//C verlet layer neighbors of atoms in the central box
    for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
    {
      size_t ip=gtbls.glst1[ilst];
      assert(ip>=0&&ip<vcs.np);
      int istatip=vcs.istat[ip];
      if(istatip<=0)continue;       // HIER: is this correct ?
      assert(istatip==1);
      int ispci=vcs.ispc[ip];
      ndbij=nbs.ndb[ip];
      getmrnb0(par,ip,ispci,istatip,ndbij,vcs,nbs);
      getmrnb1(par,ip,ispci,istatip,ndbij,vcs,nbs);
      getmrnb2(par,ip,ispci,ndbij,vcs,nbs,gtbls);
    }
    nbs.inbmrm=0;
  } // end function getmrtbls1

//======================================================================//
// Getting MR contribution for fractional SR neighbors of atom ip

  inline void getmrnb0(VLCHBOP& par,size_t ip,int ispci,int istatip,
  double ndbip,VCStr& vcs,NBStr& nbs)
  {
    int ispcc=par.ispcc,isrnb=1,imrnb;
    VmrLocDat vmrloc;

    if(istatip>0)   // MR neighbors of ip required   (HIER)
    {
      for(size_t inbij=nbs.nnbfull[ip];inbij<nbs.nnbfrac[ip];inbij++)
      {
        size_t jp=nbs.nbtbl[inbij];
        if(jp<vcs.npcb&&jp>ip){continue;}
        int istatjp=vcs.istat[jp];
        if(istatip==4)assert(abs(istatjp)<2||istatjp==4);
        if(istatip==1)assert(istatjp>-2);
        int ispcj=vcs.ispc[jp];
        double rij=nbs.rstor[inbij];
        Vec3d sigij=nbs.sigstor[inbij];
        double snij=nbs.snstor[inbij][1];
        double dsnij=nbs.snstor[inbij][2];
        double ndbij=ndbip;
        getndbij(par,ip,jp,ispci,ispcj,snij,ndbij,nbs);
        par.pRC[ispci][ispcj].rcmrij(ispci,ispcj,ispcc,ndbij,
        vmrloc.rcmr[0],vmrloc.drcmrdndb[0]);
//        assert(istatjp==4||abs(istatjp)<2);
        if(istatjp==4||istatjp==1||istatjp==0) 
        {
          double ndbji=nbs.ndb[jp];
          getndbij(par,jp,ip,ispcj,ispci,snij,ndbji,nbs);
          par.pRC[ispci][ispcj].rcmrij(ispci,ispcj,ispcc,ndbji,
          vmrloc.rcmr[1],vmrloc.drcmrdndb[1]);
          submr(par,0,ip,jp,ispci,ispcj,ispcc,rij,sigij,snij,dsnij,
          isrnb,imrnb,ndbij,vmrloc,vcs.ispc,nbs);
          if(imrnb>0){nbs.isr2mr[inbij]=nbs.nmr;}else{nbs.isr2mr[inbij]=0;}

          submr(par,1,jp,ip,ispcj,ispci,ispcc,rij,-sigij,snij,dsnij,
          isrnb,imrnb,ndbji,vmrloc,vcs.ispc,nbs);
          size_t inbji=nbs.indji[inbij];
          if(imrnb>0){nbs.isr2mr[inbji]=nbs.nmr;}else{nbs.isr2mr[inbji]=0;}
        }
        else
        {
          submr(par,0,ip,jp,ispci,ispcj,ispcc,rij,sigij,snij,dsnij,
          isrnb,imrnb,ndbij,vmrloc,vcs.ispc,nbs);
          if(imrnb>0){nbs.isr2mr[inbij]=nbs.nmr;}else{nbs.isr2mr[inbij]=0;}
        }
      }
    }  
    else         // istatip=-1 --> ip is an inactive ghost in the VSCB
    {
      for(size_t inbij=nbs.nnbfull[ip];inbij<nbs.nnbfrac[ip];inbij++)
      {
        size_t jp=nbs.nbtbl[inbij];
        if(jp<vcs.npcb&&jp>ip){continue;}
        int istatjp=vcs.istat[jp];
        if(istatjp<0||istatjp==2||istatjp==3){continue;} 
        int ispcj=vcs.ispc[jp];
        double ndbji=nbs.ndb[jp];
        double rij=nbs.rstor[inbij];
        double snij=nbs.snstor[inbij][1];
        getndbij(par,jp,ip,ispcj,ispci,snij,ndbji,nbs);
        par.pRC[ispci][ispcj].rcmrij(ispci,ispcj,ispcc,ndbji,
        vmrloc.rcmr[1],vmrloc.drcmrdndb[1]);
        Vec3d sigij=nbs.sigstor[inbij];
        double dsnij=nbs.snstor[inbij][2];
        submr(par,1,jp,ip,ispcj,ispci,ispcc,rij,-sigij,snij,dsnij,
        isrnb,imrnb,ndbji,vmrloc,vcs.ispc,nbs);
        size_t inbji=nbs.indji[inbij];
        if(imrnb>0){nbs.isr2mr[inbji]=nbs.nmr;}else{nbs.isr2mr[inbji]=0;}
      }
    }
  }     // end function getmrnb0

//======================================================================//
// Getting first verlet layer pure MR neighbor of atom ip from candidate list

  inline void getmrnb1(VLCHBOP& par,size_t ip,int ispci,int istatip,
  double ndbij,VCStr& vcs,NBStr& nbs)
  {
    for(size_t inb=nbs.nnbmrcn0[ip];inb<nbs.nnbmrcn[ip];inb++)
    {
      int imrnb=0;
      size_t jp=nbs.mrcn[inb];
      if(jp<vcs.npcb){assert(jp<ip);}      // HIER 
      int ispcj=vcs.ispc[jp];
      int istatjp=vcs.istat[jp];
//      if(ip==vcs.jpt)
//      {
//        printf("getmrnb1; ip,jp,istatip,istatjp :%6zu%6zu%6d%6d\n",ip,jp,istatip,istatjp);
////        printf("ip,jp,istatjp,nmr :%6zu%6d%6d%6zu\n",ip,jp,istatjp,nbs.nmr);
////        assert(istatjp==4||istatjp==0||abs(istatjp)==1||abs(istatjp)==3);
//      }
      double rij=nbs.rijmrcn[inb];
      Vec3d vij=nbs.vijmrcn[inb];
// (ip outside cb) or (ip in cb and (jp inactive ghost or not a first nb))   
//      if(ip==180)
////      if(ip==11&&jp==360||jp==11&&ip==360)
//      {
//        printf("getmrnb1; ip,jp,iab,istatjp :%6zu%6zu%6d%6d\n",ip,jp,iab,istatjp);
//        cin.get();
//      }  

      if(istatip>0)     // istatip=1 or 4
      {
//        if(ip==43&&jp==377)
//        if((ip==44&&jp==389)||(ip==389&&jp==44))
//        {
//          printf("YES; ip,jp,istatjp :%6zu%6zu%6d\n",ip,jp,istatjp);
//          abort();
//        }
//        assert(abs(istatjp)!=2);      // HIER
        if(ip<vcs.npcb&&(istatjp==0||istatjp==1||istatjp==4))
        {
          getmrnbij(par,ip,jp,ispci,ispcj,1,imrnb,ndbij,rij,vij,vcs,nbs);
//          if(ip==268&&jp==251)
////          if(ip==23&&jp==260)
//          {
//            printf("ip,jp,imrnb :%6zu%6zu%6d\n",ip,jp,imrnb);
//            cin.get();
//          }
          if(imrnb>0){if(istatip==4&&istatjp==0)vcs.istat[jp]=1;}    // HIER CHECKEN
//          if(imrnb>0&&jp>=vcs.npcb)vcs.istat[jp]=1;    // HIER CHECKEN
        }
        else
        {
          getmrnbij(par,ip,jp,ispci,ispcj,0,imrnb,ndbij,rij,vij,vcs,nbs);
        }  
      }  
      else    // istatip=-1
      {
        if(istatjp==0||istatjp==1||istatjp==4)
        {
          double ndbji=nbs.ndb[jp];
          getmrnbij(par,jp,ip,ispcj,ispci,0,imrnb,ndbji,rij,-vij,vcs,nbs);
//          getmrnbij(par,ip,jp,ispci,ispcj,0,imrnb,ndbij,rij,vij,vcs,nbs);
        }  
      }

    }
  }   // end getmrnb1

//======================================================================//

  inline void getmrnb2(VLCHBOP& par,size_t ip,int ispci,
  double ndbij,VCStr& vcs,NBStr& nbs,GhostTbls& gtbls)
  {
//    size_t icl,jcl,npclj,ipc,jp0;
//    Vec3d dfrc;
//    vector<double*>pri,prj;
//    pri.resize(3);  
//    prj.resize(3);  
  
    int imrnb;
    size_t icl=vcs.ivclp[ip];
    size_t ip0=vcs.ip0vcl[icl];
    int ipc=ip-ip0;
    const Vec3d ri=vcs.atpos[icl][ipc];
    for(int inbc=vcs.nnbc[1];inbc<vcs.nnbc[2];inbc++)  // second layer neighbors
    {
      size_t jcl=icl+vcs.ijcl[inbc];
      assert(jcl>=0&&jcl<vcs.nvcl);
      if(vcs.incb[jcl]>0)continue;
      int npclj=vcs.pmapinv[jcl].size();
      if(npclj==0)continue;
      size_t jp0=vcs.ip0vcl[jcl];
      for(int jpc=0;jpc<npclj;jpc++) 
      {
        size_t jp=jp0+jpc;
        int ispcj=vcs.ispc[jp];
        Vec3d vij=vcs.atpos[jcl][jpc]-ri;
        vij=vcs.hmat*vij;
        double rijsq=norm2(vij);
        if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq)
        {
          double rij=sqrt(rijsq);
          getmrnbij(par,ip,jp,ispci,ispcj,0,imrnb,ndbij,rij,vij,vcs,nbs);
//          getmrnbij(par,ip,jp,npcb,ispci,ispcj,0,imrnb,ndbij,rij,vij,nbs);
        }
      }   // end loops jpc
    }   // end loop inbc

    for(int inbc=vcs.nnbc[2]-1;inbc>=vcs.nnbc[1];inbc--)   // second layer neighbors
    {
      size_t jcl=icl-vcs.ijcl[inbc];
      assert(jcl>=0&&jcl<vcs.nvcl);
      if(vcs.incb[jcl]>0)continue;
      int npclj=vcs.pmapinv[jcl].size();
      if(npclj==0)continue;
      size_t jp0=vcs.ip0vcl[jcl];
      for(int jpc=0;jpc<npclj;jpc++) 
      {
        size_t jp=jp0+jpc;
        int ispcj=vcs.ispc[jp];
        Vec3d vij=vcs.atpos[jcl][jpc]-ri;
        vij=vcs.hmat*vij;
        double rijsq=norm2(vij);
        if(rijsq<par.pRC[ispci][ispcj].rcmrmaxsq)
        {
          double rij=sqrt(rijsq);
          getmrnbij(par,ip,jp,ispci,ispcj,0,imrnb,ndbij,rij,vij,vcs,nbs);
//          getmrnbij(par,ip,jp,npcb,ispci,ispcj,0,imrnb,ndbij,rij,vij,nbs);
        }
      }   // end loops jpc
    }   // end loop inbc
  }   // end getmrnb2

//======================================================================//
// Determine whether a MR pair ij is activated, from both sides (ij and ji).

  inline void getmrnbij(VLCHBOP& par,size_t ip,size_t jp,
  int ispci,int ispcj,int iab,int& imrnb,double ndbij,
  double rij,const Vec3d vij,VCStr& vcs,NBStr& nbs)
  {
    int ispcc=par.ispcc;
    int isrnb=0,imrnbij=0,imrnbji=0;

    if(iab>0)
    {
      VmrLocDat vmrloc;
      par.pRC[ispci][ispcj].rcmrij(ispci,ispcj,ispcc,ndbij,
      vmrloc.rcmr[0],vmrloc.drcmrdndb[0]);

      double ndbji=nbs.ndb[jp];
      par.pRC[ispci][ispcj].rcmrij(ispci,ispcj,ispcc,ndbji,
      vmrloc.rcmr[1],vmrloc.drcmrdndb[1]);
      
      double rcmax=std::max(vmrloc.rcmr[0],vmrloc.rcmr[1]);
      
      if(rij<rcmax)
      {
//        double rijinv=1.0/rij;
//        vij=rijinv*vij;
        Vec3d sigij=vij/rij;
        submr(par,0,ip,jp,ispci,ispcj,ispcc,rij,sigij,0.,0.,
        isrnb,imrnbij,ndbij,vmrloc,vcs.ispc,nbs);
        submr(par,1,jp,ip,ispcj,ispci,ispcc,rij,-sigij,0.,0.,
        isrnb,imrnbji,ndbji,vmrloc,vcs.ispc,nbs);
        imrnb=imrnbij*imrnbji;
//        if(ip==268)
////        if(ip==268&&jp==251)
////        if(ip==23&&jp==260)
//        {
//          printf("ip,jp,imrnbij,imrnbji,rij :%6zu%6zu%6d%6d%16.10lf\n",ip,jp,imrnbij,imrnbji,rij);
//          cin.get();
//        }
        if(imrnb>0)
        {
          nbs.indmrji[nbs.nmr-1]=nbs.nmr;
          nbs.indmrji[nbs.nmr]=nbs.nmr-1;
        }  
      }
      else{imrnb=0;}  
    }  
    else
    {
      VmrLocDat vmrloc;
      par.pRC[ispci][ispcj].rcmrij(ispci,ispcj,ispcc,ndbij,
      vmrloc.rcmr[0],vmrloc.drcmrdndb[0]);
      if(rij<vmrloc.rcmr[0])
      {
//        double rijinv=1.0/rij;
//        vij=rijinv*vij;
//        vij=vij/rij;
        Vec3d sigij=vij/rij;
        submr(par,0,ip,jp,ispci,ispcj,ispcc,rij,sigij,0.,0.,
        isrnb,imrnb,ndbij,vmrloc,vcs.ispc,nbs);
      }  
      else{imrnb=0;}  
    }  
  }

//======================================================================//
// Getting dangling bond number of all atoms

  void getndbap(VLCHBOP& par,VCStr& vcs,NBStr& nbs,GhostTbls& gtbls)
//  inline void getndbap(VLCHBOP& par,size_t npcb,vector<int8_t>& iflag,
//  NBStr& nbs,GhostTbls& gtbls)
  {
//    size_t ip,ilst;

    if(par.nspc==1)
    {
      if(par.ispcc==0)
      {
        for(size_t ip=0;ip<vcs.npcb;ip++)
        {
          if(vcs.istat[ip]<0)continue;    // HIER
          getndbc1(par,nbs,ip);
        }
        for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
        {
          size_t ip=gtbls.glst1[ilst];
          if(vcs.istat[ip]<0)continue;
//          if(iflag[ip]<0)continue;
          getndbc1(par,nbs,ip);
        }  
      }  
      else
      {  
        for(size_t ip=0;ip<vcs.npcb;ip++){nbs.ndb[ip]=1.-nbs.crd[ip];}
        for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
        {
          size_t ip=gtbls.glst1[ilst];
          if(vcs.istat[ip]<0)continue;
//          if(iflag[ip]<0)continue;
          nbs.ndb[ip]=1.0-nbs.crd[ip];
        }  
      }
    }  
    else
    {
      for(size_t ip=0;ip<vcs.npcb;ip++)
      {
        if(vcs.istat[ip]<0)continue;    // HIER
        int ispci=vcs.ispc[ip];
        getndbc2(par,ip,ispci,vcs.ispc,nbs);
      }  
//      for(size_t ip=0;ip<vcs.npcb;ip++){getndbc2(par,nbs,ip);}
      for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
      {
        size_t ip=gtbls.glst1[ilst];
//        if(ip==389)
////        if(ip==377)
//        {
//          printf("G; ip,ispci,istat :%6zu%6d%6d\n",ip,vcs.ispc[ip],vcs.istat[ip]);
//        }  
        if(vcs.istat[ip]<0)continue;
//        if(iflag[ip]<0)continue;
        int ispci=vcs.ispc[ip];
        getndbc2(par,ip,ispci,vcs.ispc,nbs);
//        getndbc2(par,nbs,ip);
      }  
    }
  }  // end getndbap

//======================================================================//

  inline void getndbc1(VLCHBOP& par,NBStr& nbs,size_t ip)
  {
    size_t inb=nbs.nnb0[ip];
    double xndb=4.;
    for(;inb<nbs.nnbfull[ip];inb++)
    {
      size_t jp=nbs.nbtbl[inb];
      double nji=nbs.crd[jp]-1.;
      double nelki=4./(nji+1.);
      xndb=xndb-fnelki(par,nelki);
    }  
    for(;inb<nbs.nnbfrac[ip];inb++)
    {
      size_t jp=nbs.nbtbl[inb];
      double fij=nbs.snstor[inb][1];
      double nji=nbs.crd[jp]-fij;
      double nelki=4./(nji+1.);
      xndb=xndb-fij*fnelki(par,nelki);
    }  
    nbs.ndb[ip]=xndb;
  }  

//======================================================================//

  inline void getndbc2(VLCHBOP& par,size_t ip,int ispci,
  vector<uint8_t>& ispc,NBStr& nbs)
  {
//    int ispci=vcs.ispc[ip];
    if(ispci==par.ispcc)
    {
      double xndb=4.-nbs.crdspc[ip][par.ispch];
//      if(ip==377)
////      if(ip==389)
//      {
//        printf("ip,xndb :%6zu%16.10lf\n",ip,xndb);
//      }  
      size_t inb=nbs.nnb0[ip];
      for(;inb<nbs.nnbfull[ip];inb++)
      {  
        size_t jp=nbs.nbtbl[inb];
        int ispcj=ispc[jp];
        if(ispcj==par.ispcc)
        {
          double ncji=nbs.crdspc[jp][par.ispcc]-1.;
          double nelki=(4.-nbs.crdspc[jp][par.ispch])/(ncji+1.);
          xndb=xndb-fnelki(par,nelki);
////          if(ip==377)
//          if(ip==389)
//          {
//            printf("Full; ip,jp,ncji,nelki,xndb :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",ip,jp,ncji,nelki,xndb);
//            printf("Full; ispci,ispcj,crdspcjp_H:%6d%6d%16.10lf\n",ispci,ispcj,nbs.crdspc[jp][par.ispch]);
//          }  
        }  
      }
      for(;inb<nbs.nnbfrac[ip];inb++)
      {  
        size_t jp=nbs.nbtbl[inb];
        int ispcj=ispc[jp];
        double fij=nbs.snstor[inb][1];
        if(ispcj==par.ispcc)
        {
          double ncji=nbs.crdspc[jp][par.ispcc]-fij;
          double nelki=(4.-nbs.crdspc[jp][par.ispch])/(ncji+1.);
          xndb=xndb-fij*fnelki(par,nelki);
////          if(ip==377)
//          if(ip==389)
//          {
//            printf("Frac; ip,jp,ncji,nelki,xndb :%6zu%6zu%16.10lf%16.10lf%16.10lf\n",ip,jp,ncji,nelki,xndb);
//            printf("Frac; ispci,ispcj,crdspcjp_H:%6d%6d%16.10lf\n",ispci,ispcj,nbs.crdspc[jp][par.ispch]);
//          }  
        }
      }
      nbs.ndb[ip]=xndb;
    }  
    else
    {
      nbs.ndb[ip]=1.-nbs.crd[ip];
    }  
  }  

//======================================================================//

  inline void getndbij(VLCHBOP& par,size_t ip,size_t jp,int ispci,int ispcj,
  double sn,double& ndb,NBStr& nbs)
  {
////    if(ip==377)
//    if(ip==389)
//    {
//      printf("ip,jp,ispci,ispcj :%6zu%6zu%6d%6d\n",ip,jp,ispci,ispcj);
//      printf("ndb,sn :%16.10lf%16.10lf\n",ndb,sn);
//      cin.get();
//    }
    if(ispci==par.ispcc)
//    if(vcs.ispc[ip]==par.ispcc)
    {
      if(ispcj==par.ispch){ndb+=sn;}
//      if(vcs.ispc[jp]==par.ispch){ndb+=sn;}
      else
      {
        if(par.nspc==1)
        {
          double nji=nbs.crd[jp]-sn;
          double nelji=4./(nji+1.);
          ndb+=sn*fnelki(par,nelji);
        }  
        else
        {
          double nhji=nbs.crdspc[jp][par.ispch];
          double ncji=nbs.crdspc[jp][par.ispcc]-sn;
          double nelji=(4.-nhji)/(ncji+1.);
          ndb+=sn*fnelki(par,nelji);
        }
      }
    }
    else{ndb+=sn;}
    ndb=std::max(0.,ndb);
  }

//======================================================================//

  inline void submr(VLCHBOP& par,int ij,size_t ip,size_t jp,
  int ispci,int ispcj,int ispcc,double rij,const Vec3d sigij,
  double snij,double dsnij,int isrnb,int& imrnb,double ndb,
  VmrLocDat& vmrloc,vector<uint8_t>& ispc,NBStr& nbs)
  {
//C Steric entanglement
    double sumk;
    getsumk(ip,jp,sigij,sumk,nbs);
      
//C Screening parameter
    double nij=nbs.crd[ip]-snij;
    size_t nnbij=nbs.nnbfrac[ip]-nbs.nnb0[ip]-isrnb;
    par.pVMR[ispci].gamma(nnbij,nij,sumk,ndb,imrnb,vmrloc);
//    if(ij==1&&ip==nbs.jpt)
//    {
//      printf("YES3; ip,jp,isrnb,imrnb,ndb :%6zu%6zu%6d%6d%16.10lf\n",ip,jp,isrnb,imrnb,ndb);
//      printf("YES3; nnbij,nij,sumk,vmrloc.gij :%6zu%16.10lf%16.10lf%16.10lf\n",nnbij,nij,sumk,vmrloc.gij);
//      printf("YES3; nnbfrac,nnb0 :%6zu%6zu\n",nbs.nnbfrac[ip],nbs.nnb0[ip]);
//      cin.get();
//    }
//    if(nbs.nmr>517)
//    if(ip==11&&jp==360||jp==11&&ip==360)
//    {
//      FILE *MRDAT;
//      MRDAT=fopen("TMPDAT/MR.dat","a");
//      fprintf(MRDAT,"submr; ip,jp,imrnb,nbs.nmr :%6zu%6zu%6d%6zu\n",ip,jp,imrnb,nbs.nmr);
//      fclose(MRDAT);
//    }  
//    if(nbs.nmr==520)abort();

    if(imrnb>0)
    {
      nbs.nmr++;
      assert(nbs.nmr<nbs.nbndmrm);
      if(isrnb==0)
      {
        if(nbs.nnbmr[ip]==0){nbs.nextind[nbs.nmr]=nbs.nmr;}  
        else
        {
          size_t lastind=nbs.indmr[ip];
          nbs.nextind[nbs.nmr]=nbs.nextind[lastind];
          nbs.nextind[lastind]=nbs.nmr;
        }
        nbs.nnbmr[ip]++;
        nbs.indmr[ip]=nbs.nmr;
        nbs.nbtblmr[nbs.nmr]=jp;
        nbs.rstormr[nbs.nmr]=rij;
        nbs.sigstormr[nbs.nmr]=sigij;
//        if(ip==268&&jp==251)
////        if(ip==23&&jp==260)
//        {
//          printf("ip,jp,nnbmr :%6zu%6zu%6zu\n",ip,jp,nbs.nnbmr[ip]);
//          cin.get();
//        }
      }

      nbs.gij[nbs.nmr]=vmrloc.gij;
      nbs.rcmr[nbs.nmr]=vmrloc.rcmr[ij];
      nbs.drcmrdndb[nbs.nmr]=vmrloc.drcmrdndb[ij];

//C derivatives w.r.t. MR quantities ndb
      if(ndb>0.){dndbij(par,ip,jp,ispci,isrnb,ispc,ndb,nbs);}
      else{nbs.nderndb[nbs.nmr]=nbs.nderndb[nbs.nmr-1];}
      double r1=par.pRC[ispci][ispcj].r1sr;
      double p2=par.pRC[ispci][ispcj].p2z;
      double drcmri=1./(vmrloc.rcmr[ij]-r1);
      double snmr,dsnmrdr,dsnmrdrc;
      fcpdownrc(rij,r1,vmrloc.rcmr[ij],drcmri,p2,snmr,dsnmrdr,dsnmrdrc);
     
//C loading MR tables
      if(imrnb==2)
      {
        nbs.smr[nbs.nmr]=1.;
        nbs.snsmr[nbs.nmr]=snmr-snij;
        nbs.dsnsmrdr[nbs.nmr]=dsnmrdr-dsnij;
        nbs.dsnsmrdndb[nbs.nmr]=dsnmrdrc*vmrloc.drcmrdndb[ij];
        nbs.ndersumk[nbs.nmr]=nbs.ndersumk[nbs.nmr-1];
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
// Additional initialisation for new code
        nbs.dsnsmrdnik[nbs.nmr]=0.0;
        nbs.dsmrdgij[nbs.nmr]=0.0;
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      }  
      else
      {
//C derivatives wrt MR quantities nij, sumk    
        double smr,dsmrdgij;
        dsumk(nnbij,ip,jp,ispci,rij,sigij,nbs);
        sqdsqup(vmrloc.gij,0.,1.,1.,smr,dsmrdgij);
        nbs.smr[nbs.nmr]=smr;
        nbs.snsmr[nbs.nmr]=smr*(snmr-snij);
        nbs.dsnsmrdr[nbs.nmr]=smr*(dsnmrdr-dsnij);
        double dsnsmrdgik=(snmr-snij)*dsmrdgij;
        nbs.dsnsmrdndb[nbs.nmr]=dsnsmrdgik*vmrloc.dgijdndb+
        smr*dsnmrdrc*vmrloc.drcmrdndb[ij];
        nbs.dsmrdgij[nbs.nmr]=dsmrdgij;
        nbs.dgijdndb[nbs.nmr]=vmrloc.dgijdndb;
        nbs.dgijdsumk[nbs.nmr]=vmrloc.dgijdsumk;
        nbs.dgijdnij[nbs.nmr]=vmrloc.dgijdnij;
        nbs.dsnsmrdnik[nbs.nmr]=dsnsmrdgik*vmrloc.dgijdnij;
        nbs.dsnsmrdsumk[nbs.nmr]=dsnsmrdgik*vmrloc.dgijdsumk;
      }
    }
  }

//======================================================================//

#ifdef PowRepl
  template <typename T>
  constexpr auto pow4(T const x) noexcept -> T 
  {
    auto const y = x * x;
    return y * y;
  }
#endif

  inline void getsumk(size_t ip,size_t jp,const Vec3d sigij,double& sumk,
  NBStr& nbs)
  {
    size_t kp;
    double cosijk,onepcosijk4;

    sumk=0.;
    size_t inb=nbs.nnb0[ip];
    for(;inb<nbs.nnbfull[ip];inb++)
    {
      kp=nbs.nbtbl[inb];
      if(kp==jp)continue;
      cosijk=dot(sigij,nbs.sigstor[inb]);
#ifdef PowRepl
      sumk+=pow4(1+cosijk);
#else
      onepcosijk4=pow(1.+cosijk,4);
      sumk+=onepcosijk4;
#endif      
    }  

    for(;inb<nbs.nnbfrac[ip];inb++)
    {
      kp=nbs.nbtbl[inb];
      if(kp==jp)continue;
      cosijk=dot(sigij,nbs.sigstor[inb]);
#ifdef PowRepl
      sumk+=nbs.snstor[inb][1]*pow4(1+cosijk);
#else      
      onepcosijk4=pow(1.+cosijk,4);
      sumk+=nbs.snstor[inb][1]*onepcosijk4;
#endif      
    }  
  }

//======================================================================//

  inline void dndbij(VLCHBOP& par,size_t ip,size_t jp,int ispci,
  int isrnb,vector<uint8_t>& ispc,double ndb,NBStr& nbs)
  {    
    size_t kp,lp,ider,ideri,iderk,iderk0,iderl;
    int ispck,ispcl;
    
    double fik,dfik,dfkl;
    double pref0;
    double nelkip,nelkiq,nelki0,nelki,dnelki;
      
//C  calculate ndb,frcndb,vrlndb 
    ider=nbs.nderndb[nbs.nmr-1]-1;
    if(ispci==par.ispcc)
    {
      ideri=ider+1;
      iderk=ideri;
      iderk0=0;
      nbs.nderatndb[ideri]=ip;
      nbs.frcndb[ideri]={0.,0.,0.};

      size_t nnbfullip=nbs.nnbfull[ip];
      for(size_t inbi=nbs.nnb0[ip];inbi<nbs.nnbfrac[ip];inbi++)
      {
        kp=nbs.nbtbl[inbi];
        if(kp==jp){continue;}
        ispck=ispc[kp];
        fik=nbs.snstor[inbi][1];
        dfik=nbs.snstor[inbi][2];
    
        if(ispck==par.ispch)
        {
          if(dfik!=0.)
          {
            iderk++;
            nbs.nderatndb[iderk]=kp;
            nbs.frcndb[ideri]=nbs.frcndb[ideri]-dfik*nbs.sigstor[inbi];
            nbs.frcndb[iderk]=dfik*nbs.sigstor[inbi];
#ifdef VIRIELLOCAL
            Vec3d vecik=dfik*nbs.sigstor[inbi];
            vecik=-0.5*nbs.rstor[inbi]*vecik;
            nbs.vrlndb[ideri]=nbs.vrlndb[ideri]+tensor(vecik,nbs.sigstor[inbi]);
            nbs.vrlndb[iderk]=nbs.vrlndb[iderk]+tensor(vecik,nbs.sigstor[inbi]);
#else 
#ifdef VIRIEL
            Vec3d vecik=dfik*nbs.sigstor[inbi];
            vecik=-nbs.rstor[inbi]*vecik;
            nbs.vrlndb[nbs.nmr]=nbs.vrlndb[nbs.nmr]+tensor(vecik,nbs.sigstor[inbi]);
#endif
#endif
          } 
        } 
        else
        {
          nelkip=4.-nbs.crdspc[kp][par.ispch];
          nelkiq=nbs.crdspc[kp][par.ispcc]-fik+1.;
          nelki0=nelkip/nelkiq;
    
//  correction to nelki to reach asymptotic limits at nelki0=1 (1.) and nelki0=3.5 (3.)         
          fdfnelki(par,nelki0,nelki,dnelki);

          if(inbi>=nnbfullip)
          {
            iderk0=++iderk;
            nbs.nderatndb[iderk]=kp;
            nbs.frcndb[ideri]=nbs.frcndb[ideri]-dfik*nelki*nbs.sigstor[inbi];
            nbs.frcndb[iderk]=dfik*nelki*nbs.sigstor[inbi];
#ifdef VIRIELLOCAL
            Vec3d vecik=dfik*nelki*nbs.sigstor[inbi];
            vecik=-0.5*nbs.rstor[inbi]*vecik;
            nbs.vrlndb[ideri]=nbs.vrlndb[ideri]+tensor(vecik,nbs.sigstor[inbi]);
            nbs.vrlndb[iderk]=nbs.vrlndb[iderk]+tensor(vecik,nbs.sigstor[inbi]);
#else
#ifdef VIRIEL
            Vec3d vecik=dfik*nelki*nbs.sigstor[inbi];
            vecik=-nbs.rstor[inbi]*vecik;
            nbs.vrlndb[nbs.nmr]=nbs.vrlndb[nbs.nmr]+tensor(vecik,nbs.sigstor[inbi]);
#endif            
#endif            
          }
          if(dnelki!=0.)
          {
            if(iderk>iderk0)
            {
              iderk0=++iderk;
              nbs.nderatndb[iderk]=kp;
              nbs.frcndb[iderk]={0.0,0.,0.};
            }
            else if(inbi<nnbfullip){nbs.nderatndb[iderk]=kp;}
            for(size_t inbk=nbs.nnbfull[kp];inbk<nbs.nnbfrac[kp];inbk++)
            {
              lp=nbs.nbtbl[inbk];
              if(lp==ip){continue;}
              iderl=++iderk;
              nbs.nderatndb[iderl]=lp;
              ispcl=ispc[lp];
              dfkl=fik*dnelki*nbs.snstor[inbk][2];
              if(ispcl==par.ispcc){pref0=-dfkl*nelki0/nelkiq;} 
              else{pref0=-dfkl/nelkiq;}
              nbs.frcndb[iderk0]=nbs.frcndb[iderk0]-pref0*nbs.sigstor[inbk];
              nbs.frcndb[iderl]=pref0*nbs.sigstor[inbk];
#ifdef VIRIELLOCAL
              Vec3d vecil=pref0*nbs.sigstor[inbk];
              vecil=0.5*nbs.rstor[inbk]*vecil;
              nbs.vrlndb[iderk0]=nbs.vrlndb[iderk0]-tensor(vecil,nbs.sigstor[inbk]);
              nbs.vrlndb[iderl]=nbs.vrlndb[iderl]-tensor(vecil,nbs.sigstor[inbk]);
#else 
#ifdef VIRIEL

              Vec3d vecil=pref0*nbs.sigstor[inbk];
              vecil=nbs.rstor[inbk]*vecil;
              nbs.vrlndb[nbs.nmr]=nbs.vrlndb[nbs.nmr]-tensor(vecil,nbs.sigstor[inbk]);
#endif
#endif
            }
          }
        }
      }  
      if(iderk>ideri){ider=iderk;}
    }
    else                 // H atom
    {
      ideri=ider+1; 
      iderk=ideri; 
      nbs.nderatndb[ideri]=ip;
      nbs.frcndb[ideri]={0.,0.,0.};
      for(size_t inbi=nbs.nnbfull[ip];inbi<nbs.nnbfrac[ip];inbi++)
      {
        kp=nbs.nbtbl[inbi];
        if(kp==jp){continue;}
    
        dfik=nbs.snstor[inbi][2];
        iderk++; 
        nbs.nderatndb[iderk]=kp;
        nbs.frcndb[ideri]=nbs.frcndb[ideri]-dfik*nbs.sigstor[inbi];
        nbs.frcndb[iderk]=dfik*nbs.sigstor[inbi];
#ifdef VIRIELLOCAL
        Vec3d vecik=dfik*nbs.sigstor[inbi];
        vecik=-0.5*nbs.rstor[inbi]*vecik;
        nbs.vrlndb[ideri]=nbs.vrlndb[ideri]+tensor(vecik,nbs.sigstor[inbi]);
        nbs.vrlndb[iderk]=nbs.vrlndb[iderk]+tensor(vecik,nbs.sigstor[inbi]);
#else 
#ifdef VIRIEL
        Vec3d vecik=dfik*nbs.sigstor[inbi];
        vecik=-nbs.rstor[inbi]*vecik;
        nbs.vrlndb[nbs.nmr]=nbs.vrlndb[nbs.nmr]+tensor(vecik,nbs.sigstor[inbi]);
#endif
#endif
      }
      if(iderk>ideri){ider=iderk;}
    }
    nbs.nderndb[nbs.nmr]=++ider;
  }

//========================================================================//

#ifdef PowRepl
  template <typename T>
  constexpr auto pow3(T const x) noexcept -> T 
  {
    auto const y = x * x;
    return x * y;
  }
#endif  

  inline void dsumk(size_t nnbij,size_t ip,size_t jp,int ispci,
  double rij,const Vec3d sigij,NBStr& nbs)
  {
    size_t kp,ider,ideri,iderj,iderk;
    double riji,rik,riki;
    double cosijk,onepcosijk,onepcosijk3,onepcosijk4;
    double pref1,pref2;
    Vec3d frcij,frcik;
      
//C sum over cosijk
    ider=nbs.ndersumk[nbs.nmr-1]-1;
  
    if(nnbij>0)
    {
      ideri=++ider;
      iderj=++ider;
      nbs.nderatsumk[ideri]=ip;
      nbs.nderatsumk[iderj]=jp;
      nbs.frcsumk[ideri]={0.,0.,0.};
      nbs.frcsumk[iderj]={0.,0.,0.};
     
      riji=1./rij;
      for(size_t inbi=nbs.nnb0[ip];inbi<nbs.nnbfrac[ip];inbi++)
      {
        kp=nbs.nbtbl[inbi];
        if(kp==jp){continue;}
   
//C  sum_S*(1+cos)**4
        cosijk=dot(sigij,nbs.sigstor[inbi]);
        onepcosijk=1.+cosijk;
#ifdef PowRepl
        onepcosijk3=pow3(onepcosijk);
        onepcosijk4=onepcosijk3*onepcosijk;
#else
        onepcosijk3=pow(onepcosijk,3);
        onepcosijk4=onepcosijk3*onepcosijk;
#endif
      
//C  forces and viriel from sum_S*(1+cos)**4
        iderk=++ider;
//        assert(ider<nnbmrm);
        nbs.nderatsumk[iderk]=kp;
  
        rik=nbs.rstor[inbi];
        riki=1./rik;
        pref1=4.*nbs.snstor[inbi][1]*onepcosijk3;
        pref2=nbs.snstor[inbi][2]*onepcosijk4;
       
        frcij=pref1*(nbs.sigstor[inbi]-cosijk*sigij)*riji;          // cosine part
        frcik=pref1*(sigij-cosijk*nbs.sigstor[inbi])*riki;
        nbs.frcsumk[ideri]=nbs.frcsumk[ideri]+(frcij+frcik);
        nbs.frcsumk[iderj]=nbs.frcsumk[iderj]-frcij;
        nbs.frcsumk[iderk]=-frcik;
#ifdef VIRIELLOCAL
        Vec3d vecij=0.5*rij*frcij;
        Vec3d vecik=0.5*rik*frcik;
        Mat3d vir1=tensor(vecij,sigij);
        Mat3d vir2=tensor(vecik,nbs.sigstor[inbi]);
        nbs.vrlsumk[ideri]=nbs.vrlsumk[ideri]+vir1+vir2;
        nbs.vrlsumk[iderj]=nbs.vrlsumk[iderj]+vir1;
        nbs.vrlsumk[iderk]=nbs.vrlsumk[iderk]+vir2;
#else 
#ifdef VIRIEL
        Vec3d vecij=rij*frcij;
        Vec3d vecik=rik*frcik;
        Mat3d vir1=tensor(vecij,sigij);
        Mat3d vir2=tensor(vecik,nbs.sigstor[inbi]);
        nbs.vrlsumk[nbs.nmr]=nbs.vrlsumk[nbs.nmr]+vir1+vir2;
#endif
#endif
  
        frcij=pref2*nbs.sigstor[inbi];                   // cutoff part
        nbs.frcsumk[ideri]=nbs.frcsumk[ideri]+frcij;
        nbs.frcsumk[iderk]=nbs.frcsumk[iderk]-frcij;
#ifdef VIRIELLOCAL
        vecik=0.5*rik*frcij;
        nbs.vrlsumk[ideri]=nbs.vrlsumk[ideri]+tensor(vecik,nbs.sigstor[inbi]);
        nbs.vrlsumk[iderk]=nbs.vrlsumk[iderk]+tensor(vecik,nbs.sigstor[inbi]);
#else 
#ifdef VIRIEL
        vecik=rik*frcij;
        nbs.vrlsumk[nbs.nmr]=nbs.vrlsumk[nbs.nmr]+tensor(vecik,nbs.sigstor[inbi]);
#endif
#endif
      }
    }
    nbs.ndersumk[nbs.nmr]=++ider;
  }

//========================================================================//
//  template <typename T>
//  constexpr auto pow2(T const x) noexcept -> T 
//  {
//    return x * x;
//  }

  double fnelki(VLCHBOP& par,double nelki)
  {
//    double gnelki;
    
    if(nelki<1.0)
    {
      if(nelki>par.nelkibnd)
      {
        double gnelki=nelki-par.nelkibnd;
        gnelki=gnelki*gnelki;
        return (par.nelkimin+3.0*gnelki);
      }  
      else{return par.nelkimin;}
    }
    else  
    {
      if(nelki>2.5)
      {
//        if(nelki<3.5){return (3.0-0.50*(pow2(3.50-nelki)));}
        if(nelki<3.5){return (3.0-0.50*(pow(3.50-nelki,2)));}
        else{return 3.0;}  
      }  
      else{return nelki;}
    }
  }  

//========================================================================//

  void fdfnelki(VLCHBOP& par,double nelki0,double& nelki,double& dnelki)
  {
    if(nelki0<1.0)
    {
      if(nelki0>par.nelkibnd)
      {
        dnelki=nelki0-par.nelkibnd;
        nelki=par.nelkimin+3.0*dnelki*dnelki;
        dnelki=6.0*dnelki;
      }  
      else
      {
        nelki=par.nelkimin;
        dnelki=0;
      }
    }
    else  
    {
      if(nelki0>2.5)
      {
        if(nelki0<3.5)
        {
          dnelki=3.5-nelki0;
          nelki=3.-0.5*dnelki*dnelki;
        }  
        else
        {
          nelki=3.;
          dnelki=0.;
        }
      }  
      else
      {
        nelki=nelki0;
        dnelki=1.;
      }
    }
  }  

//======================================================================//
// SR and MR interaction for pair ij
  void vsmrij(VLCHBOP& par,int ifracnb,int ipmrnb,int imrnb,size_t imrij,
  size_t imrji,vector<uint8_t>& ispc,LNBStr* nbsij,FATPLoc& Ftbls,NBStr& nbs)
  {
//    FILE *VIRBIJDAT;    // for testing only
    double vsrr=0,dvsrrdr=0,vsmra=0,dvsmradr=0;
    double fcnij=0,pchij=0,phcij=0,phhij=0;
    double pccij[2]={0,0};
    double bij,bji,bijtot;
    double vsmr,dvsmr;

//    VIRBIJDAT=fopen("TMPDAT/virbij.dat","a");    // for testing only

#ifdef VIRIELLOCAL
    for(int ij=0;ij<2;ij++)
    {
      for(int inb=0;inb<nbsij[0].nnb;inb++)
      {
        nbsij[ij].virielbij[inb]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
      }  
    }  
#else 
#ifdef VIRIEL
    for(int ij=0;ij<2;ij++)
    {
      nbsij[ij].virielbij[0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
    }  
#endif
#endif
      
    size_t ip=nbsij[0].nbtbl[0];
    size_t jp=nbsij[1].nbtbl[0];
    int ispci=nbsij[0].ispc[0];
    int ispcj=nbsij[1].ispc[0];
    double rij=nbsij[0].rstor[nbsij[0].inbj];
    Vec3d sigij=nbsij[0].sigstor[nbsij[0].inbj];

//C Get repulsive and attractive pair potentials
    subvsmrij(par,ispci,ispcj,imrnb,imrij,imrji,rij,
    nbsij[0].rcmreff,vsrr,dvsrrdr,vsmra,dvsmradr,nbsij,nbs);

//    if(ip==nbs.ipt||jp==nbs.ipt)
//    {
////      size_t ip=nbsij[0].nbtbl[0]; 
////      size_t jp=nbsij[0].nbtbl[nbsij[0].inbj]; 
//      printf("ip,jp :%6d%6d\n",nbsij[0].nbtbl[0],nbsij[0].nbtbl[nbsij[0].inbj]);
//      printf("jp,ip :%6d%6d\n",nbsij[1].nbtbl[0],nbsij[1].nbtbl[nbsij[1].inbj]);
//      cin.get();
//    }  

//    printf("OK1; bij,bji              : %16.10lf%16.10lf\n",bij,bji);
    subbij(par,0,ispci,ispcj,imrij,rij,sigij,bij,nbsij[0],nbs);
    subbij(par,1,ispcj,ispci,imrji,rij,-sigij,bji,nbsij[1],nbs);
//    printf("OK2; bij,bji              : %16.10lf%16.10lf\n",bij,bji);
    
//C Get specie pair dependent bond order corrections
    if(ispci==par.ispcc&&ispcj==par.ispcc)
    {
      par.pFATP.subFcnj_CH(par.ispcc,par.ispch,ifracnb,ipmrnb,rij,-sigij,
      vsmra,ispc,nbsij,fcnij,pccij,Ftbls,nbs);
    }    
    else if(ispci==par.ispcc&&ispcj==par.ispch)
    {
      par.pFATP.subPCH(0,par.ispcc,par.ispch,ifracnb,ipmrnb,vsmra,pchij,nbsij,Ftbls,nbs);
      par.pFATP.subPHC(1,ifracnb,vsmra,phcij,nbsij,nbs);
    }   
    else if(ispci==par.ispch&&ispcj==par.ispcc)
    {
      par.pFATP.subPHC(0,ifracnb,vsmra,phcij,nbsij,nbs);
      par.pFATP.subPCH(1,par.ispcc,par.ispch,ifracnb,ipmrnb,vsmra,pchij,nbsij,Ftbls,nbs);
    }
    else if(ispci==par.ispch&&ispcj==par.ispch)
    {
      par.pFATP.subPHH(ifracnb,vsmra,phhij,nbsij,nbs);
      par.pFATP.subFcnj_HH(ifracnb,vsmra,fcnij,nbsij,Ftbls,nbs);
    }  
      
//C Collect terms  
    bijtot=0.5*(bij+bji);
    bijtot+=pccij[0]+pccij[1]+phcij+pchij+phhij+fcnij;
    vsmr=vsrr-bijtot*vsmra;
    nbs.esrpp[ip]+=0.5*vsmr;
    nbs.esrpp[jp]+=0.5*vsmr;
//    if(nbs.ipmapinv[ip]==0||nbs.ipmapinv[jp]==0)
//    {
////      printf("icl,ipc,jcl,jpc,vsmr : %6zu%6zu%6zu%6zu%20.10lf\n",icl,ipc,jcl,jpc,vsmr);
//      printf("ip,jp,imrnb,vsmr     : %6d%6d%6d%20.10lf\n",nbs.ipmapinv[ip],nbs.ipmapinv[jp],imrnb,vsmr);
//      printf("ispci,ispcj          : %6d%6d\n",ispci,ispcj);
//      printf("rij,bijtot           : %16.10lf%16.10lf\n",rij,bijtot);
//      printf("bij,bji              : %16.10lf%16.10lf\n",bij,bji);
//      printf("nbsij[0].rcmreff     : %16.10lf\n",nbsij[0].rcmreff);
//      printf("vsrr,vsmra           : %16.10lf%20.10lf\n",vsrr,vsmra);
//      printf("pccij,fcnij          : %16.10lf%16.10lf%16.10lf\n",pccij[0],pccij[1],fcnij);
//      printf("pchij,phcij          : %16.10lf%16.10lf\n",pchij,phcij);
//      printf("\n");
//      cin.get();
//    }  
 
//C size_t iptest=0;
      
//C SR pair contributions to forces op atoms ip and jp
    dvsmr=dvsrrdr-bijtot*dvsmradr;
    nbs.force[ip]=nbs.force[ip]+dvsmr*sigij;
    nbs.force[jp]=nbs.force[jp]-dvsmr*sigij;
//    prtFrcCntr(1,ip,dvsmr*sigij);
//    prtFrcCntr(1,jp,-dvsmr*sigij);
#ifdef VIRIELLOCAL
    double hdvsmrrij=0.5*dvsmr*rij;
    nbs.viriel[ip]=nbs.viriel[ip]+hdvsmrrij*tensor(sigij,sigij);
    nbs.viriel[jp]=nbs.viriel[jp]+hdvsmrrij*tensor(sigij,sigij);
#else 
#ifdef VIRIEL
    double dvsmrrij=dvsmr*rij;
    nbs.viriel[0]=nbs.viriel[0]+dvsmrrij*tensor(sigij,sigij);
//    if(ip==0&&jp==2)
//    {
//      FILE *VIRDAT;
//      VIRDAT=fopen("TMPDAT/vir_ip.dat","a");
//      double cfe=1.03641882007443324881e-04;
//      fprintf(VIRDAT,"rij,dvsrrdr,dvsmradr,bijtot,dvsmr =%16.8lf%20.8lf%20.8lf%16.8lf%20.8lf\n",rij,dvsrrdr,dvsmradr,bijtot,dvsmr);
//      fprintf(VIRDAT,"%16.8lf%16.8lf%16.8lf\n",cfe*nbs.viriel[0].m11,cfe*nbs.viriel[0].m21,cfe*nbs.viriel[0].m31);
//      fprintf(VIRDAT,"%16.8lf%16.8lf%16.8lf\n",cfe*nbs.viriel[0].m12,cfe*nbs.viriel[0].m22,cfe*nbs.viriel[0].m32);
//      fprintf(VIRDAT,"%16.8lf%16.8lf%16.8lf\n",cfe*nbs.viriel[0].m13,cfe*nbs.viriel[0].m23,cfe*nbs.viriel[0].m33);
//      fclose(VIRDAT);
//    }  
#endif
#endif

//C Non-pair symmetric contributions from derivatives of bij
#ifdef PowRepl
    double pref0=-0.25*vsmra*pow3(bij);  //  bij contributions  
#else
    double pref0=-0.25*vsmra*pow(bij,3);  //  bij contributions  
#endif
    frcbijsym(0,nbsij[0].inbj,imrij,pref0,nbsij[0],nbs);

#ifdef PowRepl
    pref0=-0.25*vsmra*pow3(bji);  //  bij contributions  
#else
    pref0=-0.25*vsmra*pow(bji,3);  //  bij contributions  
#endif
    frcbijsym(1,nbsij[1].inbj,imrji,pref0,nbsij[1],nbs);
      
//C Non-pair a-symmetric contributions from derivatives of bij
    if(nbsij[0].inbj>=nbsij[0].nnbfl)
    {
      frcbijasym(0,nbsij[0].inbj,ip,jp,imrij,bijtot,nbsij[0],nbs);
      frcbijasym(1,nbsij[1].inbj,jp,ip,imrji,bijtot,nbsij[1],nbs);
    }  

//    fclose(VIRBIJDAT);
      
 }  //    end vsmrij

//========================================================================//

  void subvsmrij(VLCHBOP& par,int ispci,int ispcj,int imrnb,
  size_t imrij,size_t imrji,double rij,double rcmreff,double& vsrr,double& dvsrrdr,
  double& vsmra,double& dvsmradr,LNBStr* nbsij,NBStr& nbs)
  {
//C initialise output data
    double smreff,dvsmradrc,dvsmradsmr;
    vsrr=dvsrrdr=0.;
    vsmra=dvsmradr=0.;
  
    if(imrnb==0)
    {
      par.pVSR[ispci][ispcj].vdvsrra(rij,rcmreff,0.,vsrr,dvsrrdr,
      vsmra,dvsmradr,dvsmradrc,dvsmradsmr);
      nbsij[0].dsmreffdgij=nbsij[1].dsmreffdgij=0.;
      nbsij[0].drcmreffdndb=nbsij[1].drcmreffdndb=0.;
      nbsij[0].drcmreffdgij=nbsij[1].drcmreffdgij=0.;
      nbsij[0].dvsmradndb=nbsij[1].dvsmradndb=0.; 
      nbsij[0].dvsmradsumk=nbsij[1].dvsmradsumk=0.;
      nbsij[0].dvsmradnij=nbsij[1].dvsmradnij=0.; 
    }
    else
    {
//C For testing only
//     size_t ip=nbsij[0].nbtbl[0]; 
//     size_t jp=nbsij[0].nbtbl[nbsij[0].inbj]; 
      
      dsrcmreff(nbs.drcmrdndb[imrij],nbs.drcmrdndb[imrji],
      nbs.dsmrdgij[imrij],nbs.dsmrdgij[imrji],smreff,nbsij);

//      if((icl==13&&ipc==49)&&(jcl==12&&jpc==57))
////      if((ip==43&&jp==377)||(ip==377&&jp==43))
////      if((ip==44&&jp==389)||(ip==389&&jp==44))
//      {
//        printf("subvsmrij; ip,jp,smreff :%6zu%6zu%18.10lf\n",ip,jp,smreff);
//        printf("subvsmrij; ispci,ispcj  :%6d%6d\n",ispci,ispcj);
//        printf("subvsmrij; smrij,smrji  :%18.10lf%18.10lf\n",nbsij[0].smr,nbsij[1].smr);
//        cin.get();
//      }
    
      par.pVSR[ispci][ispcj].vdvsrra(rij,rcmreff,smreff,vsrr,dvsrrdr,
      vsmra,dvsmradr,dvsmradrc,dvsmradsmr);

      double prefgij=dvsmradsmr*nbsij[0].dsmreffdgij+dvsmradrc*nbsij[0].drcmreffdgij;
      nbsij[0].dvsmradndb=prefgij*nbs.dgijdndb[imrij]+dvsmradrc*nbsij[0].drcmreffdndb;
      nbsij[0].dvsmradsumk=prefgij*nbs.dgijdsumk[imrij];
      nbsij[0].dvsmradnij=prefgij*nbs.dgijdnij[imrij];
    
      prefgij=dvsmradsmr*nbsij[1].dsmreffdgij+dvsmradrc*nbsij[1].drcmreffdgij;
      nbsij[1].dvsmradndb=prefgij*nbs.dgijdndb[imrji]+dvsmradrc*nbsij[1].drcmreffdndb;
      nbsij[1].dvsmradsumk=prefgij*nbs.dgijdsumk[imrji];
      nbsij[1].dvsmradnij=prefgij*nbs.dgijdnij[imrji];
    } 
  }  // end subvsmrij

//======================================================================//

inline void loadsrcmr(size_t imrij,NBStr& nbs,LNBStr& nbsij)
{
  nbsij.smr=nbs.smr[imrij];
  nbsij.rcmr=nbs.rcmr[imrij];
}

//======================================================================//

inline void getrcmreff(LNBStr* nbsij)
{
  nbsij[0].pfeff=1./(nbsij[0].smr+nbsij[1].smr);
  nbsij[0].rcmreff=nbsij[0].pfeff*(nbsij[0].smr*nbsij[0].rcmr+nbsij[1].smr*nbsij[1].rcmr);
}

//======================================================================//

void dsrcmreff(double drcmrdndbij,double drcmrdndbji,
double dsmrdgij,double dsmrdgji,double& smreff,LNBStr* nbsij)
{
//  printf("dsrcmreff; smrij :%18.10lf\n",nbsij[0].smr);
//  printf("dsrcmreff; smrji :%18.10lf\n",nbsij[1].smr);
//rcmreff
  double pf=nbsij[0].pfeff;
  smreff=2.*pf*(nbsij[0].smr*nbsij[1].smr);
//  printf("dsrcmreff; smrji*smrji :%18.10lf\n",nbsij[0].smr*nbsij[1].smr);
//  printf("dsrcmreff; smreff      :%18.10lf\n",smreff);
  
//drc/dndb = drc/drcij*drcij/dndb
  nbsij[0].drcmreffdndb=pf*nbsij[0].smr*drcmrdndbij;
  nbsij[1].drcmreffdndb=pf*nbsij[1].smr*drcmrdndbji;
  
//drc/dgij = drc/dsmrij*dsmrij/dgij
  double pfsq=nbsij[0].pfeff*nbsij[0].pfeff;
  nbsij[0].drcmreffdgij=pfsq*nbsij[1].smr*(nbsij[0].rcmr-nbsij[1].rcmr)*dsmrdgij;
  nbsij[1].drcmreffdgij=pfsq*nbsij[0].smr*(nbsij[1].rcmr-nbsij[0].rcmr)*dsmrdgji;
  nbsij[0].dsmreffdgij=2.0*pfsq*nbsij[1].smr*nbsij[1].smr*dsmrdgij;
  nbsij[1].dsmreffdgij=2.0*pfsq*nbsij[0].smr*nbsij[0].smr*dsmrdgji;
}

//======================================================================//

 void subbij(VLCHBOP& par,const uint8_t ij,int ispci,int ispcj,
 size_t imrj,double rij,Vec3d sigij,double& bij,LNBStr& nbsij,NBStr& nbs)
 {
//C local data
   int inb,inbk,inbm,inbl,ispck; 
   size_t imrk,imrl; 
   double pref0,snsmrik;
   double y,z,u,rik;
   Vec3d sigik,sigil,sigim,ddbij[2];
   double SGH,dSGHdy,dSGHdu,dSGHdz;
   double dSGHdsnsmrik,dsnsmrdrtot,dSGHdril,dSGHdnil,dSGHdsnim;
   double Gyz,dGyzdy,dGyzdz;
   double Gyzden,Gyznum,dGyzden,dGyznum,dGyz;
   double Hduz,dHdu;

//C For testing only
//   size_t ip=nbsij.nbtbl[0]; 
//   size_t jp=nbsij.nbtbl[nbsij.inbj]; 

//   if(ip==nbs.ipt&&jp=nbs,jpt)
//   {
//     printf("subbij; ij,ip,jp          :%6d%6zu%6zu\n",ij,ip,jp);
//     printf("subbij; ispci,ispcj,ispcc :%6d%6d%6d\n",ispci,ispcj,par.ispcc);
//     cin.get();
//   }
      
//C special case with only one neighbour for ip
   if(nbsij.nnb==2){bij=1.;return;}
      
//C  initialise prefactors for forces
   for(inb=0;inb<nbsij.nnb;inb++)
   {
     nbsij.dbijdndb[inb]=0.;
     nbsij.dbijdsumk[inb]=0.;
   }

   double snsmrij=nbsij.snstor[nbsij.inbj][1]+nbs.snsmr[imrj];
    
//C initialise bij-derivative list
   for(inb=0;inb<nbsij.nnb;inb++){nbsij.isdbij[inb]=0;}
   nbsij.isdbij[0]=1;
   nbsij.isdbij[nbsij.inbj]=1;
      
//C initialise sumijk
   double sumijk=1.;

//C loop over all neighbors kp of ip
   for(inbk=1;inbk<nbsij.nnb;inbk++)
   {
     if(inbk==nbsij.inbj)continue;
     if(inbk<nbsij.nnbfl)
     {
       snsmrik=1.;
     }  
     else
     {
       imrk=nbsij.isr2mr[inbk];
       snsmrik=nbsij.snstor[inbk][1]+nbs.snsmr[imrk];
       if(inbk<nbsij.nnbfr){snsmrik=nbsij.snstor[inbk][1]+nbs.snsmr[imrk];}
       else{snsmrik=nbs.snsmr[imrk];if(snsmrik==0.)continue;}
     }
     ispck=nbsij.ispc[inbk];
     rik=nbsij.rstor[inbk];
     sigik=nbsij.sigstor[inbk];
     y=dot(sigij,sigik);
     z=nbsij.crdmr-snsmrij-snsmrik;   // crdmr not defined ???? //
   
//C Gyz and Huz functions
     if(ispci==par.ispcc)
     {
       par.pGHc.GdGyz(y,z,Gyz,dGyzdy,dGyzdz);
//       if(ip==nbs.jpt&&jp==nbs.ipt)
//       {
////         size_t kcl=iclipc(nbsij.nbtbl[inbk],vcs).first;
////         size_t kpc=iclipc(nbsij.nbtbl[inbk],vcs).second;
////         printf("icl,ipc,jcl,jpc               :%6zu%6zu%6zu%6zu\n",icl,ipc,jcl,jpc);
////         printf("subbij; kcl,kpc,rij,rik       :%6zu%6zu%18.10lf%18.10lf\n",kcl,kpc,rij,rik);
//         printf("subbij; ip,jp,kp,rij,rik      :%6d%6d%6d%18.10lf%18.10lf\n",ip,jp,nbsij.nbtbl[inbk],rij,rik);
//         printf("subbij; ispci,ispcj,ispck     :%6d%6d%6d\n",ispci,ispcj,ispck);
//         printf("subbij; crdmr,snsmrij,snsmrik :%18.10lf%18.10lf%18.10lf\n",nbsij.crdmr,snsmrij,snsmrik);
//         printf("subbij; ispck,y,z             :%6d%18.10lf%18.10lf\n",ispck,y,z);
//         printf("subbij; Gyz,dGyzdy,dGyzdz     :%18.10lf%18.10lf%18.10lf\n",Gyz,dGyzdy,dGyzdz);
//         printf("\n");
//       }

       if(ispcj==par.ispch||ispck==par.ispch)
       {
         Gyznum=(par.pGHc.gCbc0[ispcj][ispck]+par.pGHc.gCbc1[ispcj][ispck]*Gyz)*Gyz;
         dGyznum=par.pGHc.gCbc0[ispcj][ispck]+2.*par.pGHc.gCbc1[ispcj][ispck]*Gyz;
         Gyzden=1.+par.pGHc.gC*Gyz;
         dGyzden=par.pGHc.gC;
     
         Gyz=Gyznum/Gyzden;
         dGyz=(dGyznum*Gyzden-Gyznum*dGyzden)/(Gyzden*Gyzden);
         dGyzdy=dGyz*dGyzdy;
         dGyzdz=dGyz*dGyzdz;

         Hduz=1.;
         dHdu=0.;
       }  
       else
       {
         double u=rij-rik;
         par.pGHc.HdHu(u,Hduz,dHdu);
       }
     }
     else
     {
       par.pGHh.GdGyz(y,z,Gyz,dGyzdy,dGyzdz);
//       if(ip==nbs.jpt&&jp==nbs.ipt)
//       {
////         size_t kcl=iclipc(nbsij.nbtbl[inbk],vcs).first;
////         size_t kpc=iclipc(nbsij.nbtbl[inbk],vcs).second;
////         printf("icl,ipc,jcl,jpc               :%6zu%6zu%6zu%6zu\n",icl,ipc,jcl,jpc);
////         printf("subbij; kcl,kpc,rij,rik       :%6zu%6zu%18.10lf%18.10lf\n",kcl,kpc,rij,rik);
//         printf("subbij; ip,jp,kp,rij,rik      :%6d%6d%6d%18.10lf%18.10lf\n",ip,jp,nbsij.nbtbl[inbk],rij,rik);
//         printf("subbij; crdmr,snsmrij,snsmrik :%18.10lf%18.10lf%18.10lf\n",nbsij.crdmr,snsmrij,snsmrik);
//         printf("subbij; ispck,y,z             :%6d%18.10lf%18.10lf\n",ispck,y,z);
//         printf("subbij; Gyz,dGyzdy,dGyzdz     :%18.10lf%18.10lf%18.10lf\n",Gyz,dGyzdy,dGyzdz);
//         printf("\n");
//       }

       if(ispcj==par.ispcc||ispck==par.ispcc)
       {
         Hduz=1.;
         dHdu=0.;
       }
       else
       {
         u=rij-rik;
         par.pGHh.HdHu(u,Hduz,dHdu);
       }
     }
      
     SGH=snsmrik*Gyz*Hduz;
     sumijk+=SGH;
      
//C Start calculating forces
     nbsij.isdbij[inbk]=1;
      
//C MR independent-part !!!
      
//C  Contribution from y=COSijk
     dSGHdy=snsmrik*Hduz*dGyzdy;
     ddbij[0]=dSGHdy*(sigik-y*sigij)/rij;
     ddbij[1]=dSGHdy*(sigij-y*sigik)/rik;
     nbsij.dbij[0]=nbsij.dbij[0]-ddbij[0]-ddbij[1];
     nbsij.dbij[nbsij.inbj]=nbsij.dbij[nbsij.inbj]+ddbij[0];
     nbsij.dbij[inbk]=nbsij.dbij[inbk]+ddbij[1];

#ifdef VIRIELLOCAL
     double hdSGHdy=0.5*dSGHdy;
     Vec3d vec1=hdSGHdy*(sigik-y*sigij);
     Vec3d vec2=hdSGHdy*(sigij-y*sigik);
     Mat3d vir1=tensor(sigij,vec1);
     Mat3d vir2=tensor(sigik,vec2);
     nbsij.virielbij[0]=nbsij.virielbij[0]-vir1-vir2;
     nbsij.virielbij[nbsij.inbj]=nbsij.virielbij[nbsij.inbj]-vir1;
     nbsij.virielbij[inbk]=nbsij.virielbij[inbk]-vir2;
#else 
#ifdef VIRIEL
     Vec3d vec1=dSGHdy*(sigik-y*sigij);
     Vec3d vec2=dSGHdy*(sigij-y*sigik);
     Mat3d vir1=tensor(sigij,vec1);
     Mat3d vir2=tensor(sigik,vec2);
     nbsij.virielbij[0]=nbsij.virielbij[0]-vir1-vir2;
#endif
#endif
      
//C Contribution from u=rij-rik-(rij0-rik0)
     dSGHdu=snsmrik*Gyz*dHdu;
     if(dSGHdu!=0.)
     {
       ddbij[0]=dSGHdu*sigij;
       ddbij[1]=dSGHdu*sigik;
       nbsij.dbij[0]=nbsij.dbij[0]-ddbij[0]+ddbij[1];
       nbsij.dbij[nbsij.inbj]=nbsij.dbij[nbsij.inbj]+ddbij[0];
       nbsij.dbij[inbk]=nbsij.dbij[inbk]-ddbij[1];
#ifdef VIRIELLOCAL
       double hdSGHdu=0.5*dSGHdu;
       vir1=hdSGHdu*rij*tensor(sigij,sigij);
       vir2=hdSGHdu*rik*tensor(sigik,sigik);
       nbsij.virielbij[0]=nbsij.virielbij[0]-vir1+vir2;
       nbsij.virielbij[nbsij.inbj]=nbsij.virielbij[nbsij.inbj]-vir1;
       nbsij.virielbij[inbk]=nbsij.virielbij[inbk]+vir2;
#else 
#ifdef VIRIEL
       vir1=dSGHdu*rij*tensor(sigij,sigij);
       vir2=dSGHdu*rik*tensor(sigik,sigik);
       nbsij.virielbij[0]=nbsij.virielbij[0]-vir1+vir2;
#endif
#endif
     }
//#ifdef VERBOSE
//      write(90,'(a,3f18.12)')"Precision ddbij 2 : ",dbij(1,0,ij),dbij(2,0,ij),dbij(3,0,ij)
//      write(90,'(a,3f18.12)')"Precision ddbij 2 : ",dbij(1,jnbi,ij),dbij(2,jnbi,ij),dbij(3,jnbi,ij)
//      write(90,'(a,3f18.12)')"Precision ddbij 2 : ",dbij(1,knbi,ij),dbij(2,knbi,ij),dbij(3,knbi,ij)
//#endif
      
//C MR dependent-part
//C  contributions from snsmril (due to snsmrik for l=k, to zijkmr=crdmr(ip)-snsmrij-snsmrik otherwise)
     dSGHdsnsmrik=Gyz*Hduz;
     dSGHdz=snsmrik*Hduz*dGyzdz;
     for(inbl=nbsij.nnbfl;inbl<nbsij.nnb;inbl++)
     {
       if(inbl==nbsij.inbj)continue;
       imrl=nbsij.isr2mr[inbl];
       dsnsmrdrtot=nbsij.snstor[inbl][2]+nbs.dsnsmrdr[imrl];
       if(inbl==inbk){pref0=dSGHdsnsmrik;}
       else{pref0=dSGHdz;}
       if(abs(pref0)<1e-12)continue;   // should not be needed (???) //
          
//C ril contribution from snsmril (can be purely SR, with isMR=0)
       assert(nbsij.isdbij[inbl]=1);
       nbsij.isdbij[inbl]=1;      //  Why ??? isdbij[inbl] should 1 already. //
       dSGHdril=pref0*dsnsmrdrtot;
       sigil=nbsij.sigstor[inbl];
       ddbij[0]=dSGHdril*sigil;
       nbsij.dbij[0]=nbsij.dbij[0]-ddbij[0];
       nbsij.dbij[inbl]=nbsij.dbij[inbl]+ddbij[0];
#ifdef VIRIELLOCAL
       double hril=0.5*nbsij.rstor[inbl];
       nbsij.virielbij[0]=nbsij.virielbij[0]-hril*tensor(ddbij[0],sigil);
       nbsij.virielbij[inbl]=nbsij.virielbij[inbl]-hril*tensor(ddbij[0],sigil);
#else 
#ifdef VIRIEL
       double ril=nbsij.rstor[inbl];
       nbsij.virielbij[0]=nbsij.virielbij[0]-ril*tensor(ddbij[0],sigil);
#endif
#endif
      
       if(imrl>0)
       {
//C other contributions from snsmril (purely MR)   // why purely MR ???? //
         nbsij.dbijdndb[inbl]+=pref0*nbs.dsnsmrdndb[imrl];
         nbsij.dbijdsumk[inbl]+=pref0*nbs.dsnsmrdsumk[imrl];
      
//  contribution from nil
         dSGHdnil=pref0*nbs.dsnsmrdnik[imrl];
         for(inbm=1;inbm<nbsij.nnbfr;inbm++) // only purely SR frac neighbours contribute
         {
           if(inbm==inbl)continue;
           assert(nbsij.isdbij[inbm]=1);
           nbsij.isdbij[inbm]=1;            // should not be needed (????) //
           dSGHdsnim=dSGHdnil*nbsij.snstor[inbm][2];
           sigim=nbsij.sigstor[inbm];
           nbsij.dbij[0]=nbsij.dbij[0]-dSGHdsnim*sigim;
           nbsij.dbij[inbm]=nbsij.dbij[inbm]+dSGHdsnim*sigim;
#ifdef VIRIELLOCAL
           Vec3d vecim=dSGHdsnim*sigim;
           vecim=-0.5*nbsij.rstor[inbm]*vecim;
           nbsij.virielbij[0]=nbsij.virielbij[0]+tensor(vecim,sigim);
           nbsij.virielbij[inbm]=nbsij.virielbij[inbm]+tensor(vecim,sigim);
#else 
#ifdef VIRIEL
           Vec3d vecim=dSGHdsnim*sigim;
           vecim=-nbsij.rstor[inbm]*vecim;
           nbsij.virielbij[0]=nbsij.virielbij[0]+tensor(vecim,sigim);
#endif
#endif
         }    // loop inbm
       }    // if(imrl>0)   
     }    // loop inbl
   }    // loop inbk

   bij=sqrt(1./sumijk);
 }    //  end subbij

//========================================================================//
// Symmetric contributions from derivatives of bij (r<rceff)

  void frcbijsym(int ij,int inbj, size_t imrj,double pref0,
  LNBStr& nbsij,NBStr& nbs)
  {
    int inbk;
    size_t imrk,ider,kp;
    double prefndb,prefsumk;

//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
// For testing
//    size_t iptest=0;
//    size_t ip=nbsij.nbtbl[0];
//    size_t jp=nbsij.nbtbl[inbj];
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC    
    
//C Bij_SR contribs
    for(inbk=0;inbk<nbsij.nnb;inbk++)
    {
//      assert(nbsij.isdbij[inbk]!=0);
      if(nbsij.isdbij[inbk]==0){continue;}
      kp=nbsij.nbtbl[inbk];
      nbs.force[kp]=nbs.force[kp]+pref0*nbsij.dbij[inbk];
//      prtFrcCntr(2,kp,pref0*nbsij.dbij[inbk]);
      nbsij.dbij[inbk]={0.,0.,0.};
#ifdef VIRIELLOCAL
        nbs.viriel[kp]=nbs.viriel[kp]+pref0*nbsij.virielbij[inbk];
        nbsij.virielbij[inbk]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
    }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
    nbs.viriel[0]=nbs.viriel[0]+pref0*nbsij.virielbij[0];
    nbsij.virielbij[0]={0.,0.,0.,0.,0.,0.,0.,0.,0.};
#endif
#endif
      
//C Bij_MR contribs
    for(inbk=nbsij.nnbfl;inbk<nbsij.nnb;inbk++)
    {
      if(inbk==inbj)continue;
      imrk=nbsij.isr2mr[inbk];
      if(imrk==0)continue;
      if(nbsij.dbijdndb[inbk]==0.)continue;  // should not be needed ??? //
      prefndb=pref0*nbsij.dbijdndb[inbk];
      for(ider=nbs.nderndb[imrk-1];ider<nbs.nderndb[imrk];ider++)
      {
        kp=nbs.nderatndb[ider];
        nbs.force[kp]=nbs.force[kp]-prefndb*nbs.frcndb[ider];
//        prtFrcCntr(3,kp,-prefndb*nbs.frcndb[ider]);
#ifdef VIRIELLOCAL
        nbs.viriel[kp]=nbs.viriel[kp]-prefndb*nbs.vrlndb[ider];
#endif
      } 
#ifdef VIRIEL
#ifndef VIRIELLOCAL
      nbs.viriel[0]=nbs.viriel[0]-prefndb*nbs.vrlndb[imrk];
#endif
#endif
    }
      
//C Sumk
    for(inbk=nbsij.nnbfl;inbk<nbsij.nnb;inbk++)
    {
      if(inbk==inbj)continue;
      imrk=nbsij.isr2mr[inbk];
      if(imrk==0)continue;
      if(nbsij.dbijdsumk[inbk]==0.)continue;  // should not be needed ??? //
      prefsumk=pref0*nbsij.dbijdsumk[inbk];
      for(ider=nbs.ndersumk[imrk-1];ider<nbs.ndersumk[imrk];ider++)
      {
        kp=nbs.nderatsumk[ider];
        nbs.force[kp]=nbs.force[kp]-prefsumk*nbs.frcsumk[ider];
//        prtFrcCntr(4,kp,-prefsumk*nbs.frcsumk[ider]);
#ifdef VIRIELLOCAL
        nbs.viriel[kp]=nbs.viriel[kp]-prefsumk*nbs.vrlsumk[ider];
#endif
      }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
      nbs.viriel[0]=nbs.viriel[0]-prefsumk*nbs.vrlsumk[imrk];
#endif
#endif
    }
  }  // end frcbijsym

//========================================================================//
//  asymmetric contributions to forces (due to LRcorr contributions): 
//  need to compute for r>rceff; only from MR contributions to Vsmra

  void frcbijasym(int ij,int inbj,size_t ip,size_t jp,size_t imrj,
  double bijtot,LNBStr& nbsij,NBStr& nbs)
  {
    size_t ider,kp;
//    size_t iptest=-1;
      
// VsmrA_Ndb
    if(imrj>0)
    {
      if(nbsij.dvsmradndb!=0.)  // Is this needed ??? //
      {
        nbsij.dvsmradndb=-nbsij.dvsmradndb*bijtot;  // SR contrib
        {
          for(ider=nbs.nderndb[imrj-1];ider<nbs.nderndb[imrj];ider++)
          {
            kp=nbs.nderatndb[ider];
            nbs.force[kp]=nbs.force[kp]+nbsij.dvsmradndb*nbs.frcndb[ider];
//            prtFrcCntr(5,kp,nbsij.dvsmradndb*nbs.frcndb[ider]);
#ifdef VIRIELLOCAL
            nbs.viriel[kp]=nbs.viriel[kp]+nbsij.dvsmradndb*nbs.vrlndb[ider];
#endif
          }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
          nbs.viriel[0]=nbs.viriel[0]+nbsij.dvsmradndb*nbs.vrlndb[imrj];
#endif
#endif
        }
      }
      
// VsmrA_Sumk
      if(nbsij.dvsmradsumk!=0.)
      {
        nbsij.dvsmradsumk=-nbsij.dvsmradsumk*bijtot;
        for(ider=nbs.ndersumk[imrj-1];ider<nbs.ndersumk[imrj];ider++)
        {
          kp=nbs.nderatsumk[ider];
          nbs.force[kp]=nbs.force[kp]+nbsij.dvsmradsumk*nbs.frcsumk[ider];
//          prtFrcCntr(6,kp,nbsij.dvsmradsumk*nbs.frcsumk[ider]);
#ifdef VIRIELLOCAL
          nbs.viriel[kp]=nbs.viriel[kp]+nbsij.dvsmradsumk*nbs.vrlsumk[ider];
#endif
        }
#ifdef VIRIEL
#ifndef VIRIELLOCAL
        nbs.viriel[0]=nbs.viriel[0]+nbsij.dvsmradsumk*nbs.vrlsumk[imrj];
#endif
#endif
      }
    }
      
//C VsmrA_Nij
    if(abs(nbsij.dvsmradnij)>1e-12)
    {
      nbsij.dvsmradnij=-nbsij.dvsmradnij*bijtot;             // SR contrib (pair pot.)
      for(int inbk=nbsij.nnbfl;inbk<nbsij.nnb;inbk++)
      {
        if(inbk==inbj)continue;
        if(abs(nbsij.snstor[inbk][2])<1e-12)continue;   // if would be better to limit loop over frac nbs ???? //
        kp=nbsij.nbtbl[inbk];
        assert(kp!=jp);
        if(kp!=jp)        // should not be needed ?????? //
        {
          Vec3d vecik=nbsij.dvsmradnij*nbsij.snstor[inbk][2]*nbsij.sigstor[inbk];
          nbs.force[ip]=nbs.force[ip]+vecik;
          nbs.force[kp]=nbs.force[kp]-vecik;
//          prtFrcCntr(7,ip,vecik);
//          prtFrcCntr(7,kp,-vecik);
#ifdef VIRIELLOCAL
          vecik=0.5*nbsij.rstor[inbk]*vecik;
          nbs.viriel[ip]=nbs.viriel[ip]+tensor(vecik,nbsij.sigstor[inbk]);
          nbs.viriel[kp]=nbs.viriel[kp]+tensor(vecik,nbsij.sigstor[inbk]);
#else 
#ifdef VIRIEL
          vecik=nbsij.rstor[inbk]*vecik;
          nbs.viriel[0]=nbs.viriel[0]+tensor(vecik,nbsij.sigstor[inbk]);
#endif
#endif
        }
      }
    }
  }  // end frcbijasym

//========================================================================//

  inline void loadnbssp(size_t ip,size_t imr,vector<uint8_t>& ispc,
  LNBStr& nbsij,NBStr& nbs)
  {
    int inb=0;
    size_t jp;
    nbsij.nbtbl[inb]=ip;
    nbsij.ispc[inb++]=ispc[ip];
    size_t inb0=nbs.nnb0[ip];
    nbsij.nnbfl=nbs.nnbfull[ip]-inb0+1;
    nbsij.nnbfr=nbs.nnbfrac[ip]-inb0+1;
    nbsij.nnb=nbsij.nnbfr+nbs.nnbmr[ip];
    nbsij.crdmr=nbs.crd[ip];
    for(;inb0<nbs.nnbfull[ip];inb0++,inb++)
    {
      assert(inb<nbsij.nnbmloc);
      jp=nbs.nbtbl[inb0];
      nbsij.nbtbl[inb]=jp;
      nbsij.ispc[inb]=ispc[jp];
      nbsij.isr2mr[inb]=nbs.isr2mr[inb0];
      nbsij.rstor[inb]=nbs.rstor[inb0];
      nbsij.sigstor[inb]=nbs.sigstor[inb0];
      nbsij.snstor[inb][0]=nbs.snstor[inb0][0];
      nbsij.snstor[inb][1]=nbs.snstor[inb0][1];
      nbsij.snstor[inb][2]=nbs.snstor[inb0][2];
    }
    for(;inb0<nbs.nnbfrac[ip];inb0++,inb++)
    {
      jp=nbs.nbtbl[inb0];
      assert(inb<nbsij.nnbmloc);
      nbsij.nbtbl[inb]=jp;
      nbsij.ispc[inb]=ispc[jp];
      nbsij.isr2mr[inb]=nbs.isr2mr[inb0];
      nbsij.rstor[inb]=nbs.rstor[inb0];
      nbsij.sigstor[inb]=nbs.sigstor[inb0];
      nbsij.snstor[inb][0]=nbs.snstor[inb0][0];
      nbsij.snstor[inb][1]=nbs.snstor[inb0][1];
      nbsij.snstor[inb][2]=nbs.snstor[inb0][2];
      nbsij.crdmr+=nbs.snsmr[nbsij.isr2mr[inb]];
    }

    for(inb0=0;inb0<nbs.nnbmr[ip];inb0++,inb++)
    {
      jp=nbs.nbtblmr[imr];
      assert(inb<nbsij.nnbmloc);
      nbsij.nbtbl[inb]=jp;
      nbsij.ispc[inb]=ispc[jp];
      nbsij.isr2mr[inb]=imr;
      nbsij.rstor[inb]=nbs.rstormr[imr];
      nbsij.sigstor[inb]=nbs.sigstormr[imr];
      nbsij.snstor[inb][0]=-1.;
      nbsij.snstor[inb][1]=0.;
      nbsij.snstor[inb][2]=0.;
      nbsij.crdmr+=nbs.snsmr[imr];
      imr=nbs.nextind[imr];
    }  
  }  

//=====================================================================//

  inline void addvlr(size_t ip,size_t jp,double vlr,Vec3d vij,Vec3d dfrc,NBStr& nbs)
  {
//    return;
    nbs.elrpp[ip]=nbs.elrpp[ip]+0.5*vlr;
    nbs.elrpp[jp]=nbs.elrpp[jp]+0.5*vlr;
    nbs.force[ip]=nbs.force[ip]+dfrc;
    nbs.force[jp]=nbs.force[jp]-dfrc;
//    if(ip==nbs.ipt||jp==nbs.ipt)
//    {
//      printf("addvlr; ip,jp,vlr :%6zu%6zu%16.10lf\n",ip,jp,0.5*vlr);
//      cin.get();
//    }

#ifdef VIRIELLOCAL
      nbs.viriel[ip]=nbs.viriel[ip]+tensor(0.5*vij,dfrc);
      nbs.viriel[jp]=nbs.viriel[jp]+tensor(0.5*vij,dfrc);
#else
#ifdef VIRIEL
      nbs.viriel[0]=nbs.viriel[0]+tensor(vij,dfrc);
#endif
#endif
  }  

//======================================================================//

  inline void chkaddvlr(int ichkvlr,size_t ip,size_t jp,double rij,
  Vec3d vij,double vlr,VCStr& vcs,NBStr& nbs)
  {
    size_t ivcl=vcs.ivclp[ip];
    int ipvc=ip-vcs.ip0vcl[ivcl];
    size_t icl=vcs.pmapinv[ivcl][ipvc].first;
    size_t ipc=vcs.pmapinv[ivcl][ipvc].second;

    size_t jvcl=vcs.ivclp[jp];
    int jpvc=jp-vcs.ip0vcl[jvcl];
    size_t jcl=vcs.pmapinv[jvcl][jpvc].first;
    size_t jpc=vcs.pmapinv[jvcl][jpvc].second;

    if(icl==1&&ipc==142)
    {
      size_t ivcl=vcs.ivclp[jp];
      int ipvc=jp-vcs.ip0vcl[ivcl];
      size_t icl=vcs.pmapinv[ivcl][ipvc].first;
      size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
      printf("1,ichkvlr,icl,ipc,istatjp,rij,vlr :%6d%6zu%6zu%6d%20.12lf%20.12lf\n",ichkvlr,icl,ipc,vcs.istat[jp],rij,vlr);
    }
    else if(jcl==1&&jpc==142)
    {
      size_t ivcl=vcs.ivclp[ip];
      int ipvc=ip-vcs.ip0vcl[ivcl];
      size_t icl=vcs.pmapinv[ivcl][ipvc].first;
      size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
      printf("2,ichkvlr,icl,ipc,istatjp,rij,vlr :%6d%6zu%6zu%6d%20.12lf%20.12lf\n",ichkvlr,icl,ipc,vcs.istat[jp],rij,vlr);
    }
  }  

//======================================================================//

  void prtFrcCntr(int icntr,size_t ip,Vec3d dfrc)
  {
    FILE *FRCCNTR;
//    FRCCNTR=fopen("TMPDAT/frccntr.dat","w");
    FRCCNTR=fopen("TMPDAT/frccntr.dat","a");
    fprintf(FRCCNTR,"%5d%8zu%18.10lf%18.10lf%18.10lf\n",icntr,ip,dfrc.x,dfrc.y,dfrc.z);
    fclose(FRCCNTR);
  }

//======================================================================//

  void printEsr(size_t n,int iout,double esr,vector<double>& esrpp,VCStr& vcs)
  {
    FILE *ESRDAT;
    if(iout==0)
    {
      ESRDAT=fopen("TMPDAT/esrpp.dat","w");
      fprintf(ESRDAT,"esr = %18.10lf\n",esr);
      for(size_t ip=0;ip<n;ip++)
      {
        fprintf(ESRDAT,"%8zu%18.10lf\n",ip,esrpp[ip]);
//        size_t ivcl=vcs.ivclp[ip];
//        int ipvc=ip-vcs.ip0vcl[ivcl];
//        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
//        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
//        fprintf(ESRDAT,"%8zu%18.10lf\n",ipc,esrpp[ipc]);
//        size_t ipp=vcs.ipmap[ip];
//        fprintf(ESRDAT,"%8zu%8zu%18.10lf\n",ip,ipp,esrpp[ipp]);
      }  
      fclose(ESRDAT);
    }
    else if(iout==1)
    {
      ESRDAT=fopen("TMPDAT/esrpcl.dat","w");
      fprintf(ESRDAT,"esr = %18.10lf\n",esr);
      for(size_t ip=0;ip<n;ip++)
      {
        size_t ivcl=vcs.ivclp[ip];
        int ipvc=ip-vcs.ip0vcl[ivcl];
        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
//        if(icl==13&&ip>vcs.npcb-1)
//        {
//          size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
//          printf("%8zu%8zu%8zu%18.10lf\n",icl,ipc,ip,esrpp[ip]);
////          printf("%8zu%8d%8zu%18.10lf\n",ivcl,ipvc,ip,esrpp[ip]);
//          abort();
//        }
        if(icl!=13)continue;
//        if(ip>vcs.npcb-1)continue;
        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
        fprintf(ESRDAT,"%8zu%8zu%8zu%18.10lf\n",icl,ipc,ip,esrpp[ip]);
//        fprintf(ESRDAT,"%8zu%8d%8zu%18.10lf\n",ivcl,ipvc,ip,esrpp[ip]);
      }  
      fclose(ESRDAT);
    }
//    else if(iout==2)
//    {
//      ESRDAT=fopen("TMPDAT/esrpcl.dat","w");
//      fprintf(ESRDAT,"esr = %18.10lf\n",esr);
//      for(ssize_t iclk=vcs.ngls0;iclk<vcs.dims.k-vcs.ngls0;iclk++)
//      for(ssize_t iclj=vcs.ngls0;iclj<vcs.dims.j-vcs.ngls0;iclj++)
//      for(ssize_t icli=vcs.ngls0;icli<vcs.dims.i-vcs.ngls0;icli++)
//      {
//        IJK iclijk{icli,iclj,iclk};
//        size_t icl=grid_ijk_to_index(vcs.dims,iclijk);
//        for(size_t ipc=0;ipc<cells[icl].size();ipc++)
//        {
//          size_t ivcl=vcs.pmap[icl][ipc].first;
//          int ipvc=vcs.pmap[icl][ipc].second;
//          size_t ip=vcs.ip0vcl[ivcl]+ipvc;
//          fprintf(ESRDAT,"%8zu%8zu%8zu%18.10lf\n",icl,ipc,ip,esrpp[ip]);
//        }  
//      }  
//    }
//    fclose(ESRDAT);
  }

//======================================================================//

  void printElr(size_t n,int iout,double elr,vector<double>& elrpp,VCStr& vcs)
  {
    FILE *ELRDAT;
    if(iout==0)
    {
      ELRDAT=fopen("TMPDAT/elrpp.dat","w");
      fprintf(ELRDAT,"elr = %18.10lf\n",elr);
      for(size_t ip=0;ip<n;ip++)
      {
        fprintf(ELRDAT,"%8zu%18.10lf\n",ip,elrpp[ip]);
//        size_t ivcl=vcs.ivclp[ip];
//        int ipvc=ip-vcs.ip0vcl[ivcl];
//        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
//        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
//        Vec3d rip=vcs.atpos[ivcl][ipvc]+vcs.xyz0;
//        fprintf(ELRDAT,"%8zu%18.10lf%18.10lf%18.10lf%18.10lf\n",ipc,elrpp[ipc],rip.x,rip.y,rip.z);
//        size_t ipp=vcs.ipmap[ip];
//        fprintf(ELRDAT,"%8zu%8zu%18.10lf\n",ip,ipp,elrpp[ipp]);
      }  
      fclose(ELRDAT);
    }
    else if(iout==1)
    {
      ELRDAT=fopen("TMPDAT/elrpcl.dat","w");
      fprintf(ELRDAT,"elr = %18.10lf\n",elr);
      for(size_t ip=0;ip<n;ip++)
      {
        size_t ivcl=vcs.ivclp[ip];
        int ipvc=ip-vcs.ip0vcl[ivcl];
        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
//        if(icl!=13)continue;
        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
        fprintf(ELRDAT,"%8zu%8zu%18.10lf\n",icl,ipc,elrpp[ip]);
      }  
      fclose(ELRDAT);
    }
    else if(iout==2)
    {
      ELRDAT=fopen("TMPDAT/elrpcl.dat","w");
      fprintf(ELRDAT,"elr = %18.10lf\n",elr);
      for(size_t ip=0;ip<n;ip++)
      {
        size_t ivcl=vcs.ivclp[ip];
        int ipvc=ip-vcs.ip0vcl[ivcl];
        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
//        if(icl!=13)continue;
        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
        fprintf(ELRDAT,"%8zu%8zu%8zu%18.10lf\n",icl,ipc,ip,elrpp[ip]);
      }  
      fclose(ELRDAT);
    }
  }

//======================================================================//

  void printEF(size_t n,int iout,double esr,double elr,vector<Vec3d>& force,VCStr& vcs)
  {
    FILE *FRCDAT;
    FRCDAT=fopen("TMPDAT/frc.dat","w");
//    FRCDAT=fopen("TMPDAT/frc.dat","a");
//    fprintf(FRCDAT,"Esr,Elr,Etot :%18.9lf%18.9lf%18.9lf\n",esr,elr,esr+elr);
//    fprintf(FRCDAT,"Esr,Elr,Etot :%18.10lf%18.10lf%18.10lf\n",esr,elr,esr+elr);
    fprintf(FRCDAT,"%20.10lf%20.10lf%20.10lf\n",esr,elr,esr+elr);

    if(iout==0)
    {
      for(size_t ip=0;ip<n;ip++)
      {
//        printf("ip,frcip :%6zu%16.8lf\n",ip,force[ip].x,force[ip].y,force[ip].z);
//        fprintf(FRCDAT,"%8zu%18.9lf%18.9lf%18.9lf\n",ip,force[ip].x,force[ip].y,force[ip].z);
        fprintf(FRCDAT,"%8zu%18.10lf%18.10lf%18.10lf\n",ip,force[ip].x,force[ip].y,force[ip].z);
      }  
    }  
    else if(iout==1)
    {
      for(size_t ip=0;ip<n;ip++)
      {
        size_t ivcl=vcs.ivclp[ip];
        int ipvc=ip-vcs.ip0vcl[ivcl];
        size_t icl=vcs.pmapinv[ivcl][ipvc].first;
        size_t ipc=vcs.pmapinv[ivcl][ipvc].second;
        fprintf(FRCDAT,"%8zu%8zu%20.10lf%20.10lf%20.10lf\n",icl,ipc,force[ip].x,force[ip].y,force[ip].z);
      }  
    }  
    fclose(FRCDAT);
  }

//======================================================================//

  void printVIR(size_t n,vector<Mat3d>& viriel)
  {

#ifdef VIRIELLOCAL
    FILE *VIRDAT;
    VIRDAT=fopen("TMPDAT/vir.dat","w");
    Mat3d virtot=viriel[0];
    for(size_t ip=1;ip<n;ip++){virtot=virtot+viriel[ip];}
    fprintf(VIRDAT,"%18.10lf%18.10lf%18.10lf\n",virtot.m11,virtot.m21,virtot.m31);
    fprintf(VIRDAT,"%18.10lf%18.10lf%18.10lf\n",virtot.m12,virtot.m22,virtot.m32);
    fprintf(VIRDAT,"%18.10lf%18.10lf%18.10lf\n",virtot.m13,virtot.m23,virtot.m33);
    fclose(VIRDAT);

    VIRDAT=fopen("TMPDAT/virloc.dat","w");
    for(size_t ip=0;ip<n;ip++)
    {
      fprintf(VIRDAT,"%8zu\n",ip);
      fprintf(VIRDAT,"%16.8lf%16.8lf%16.8lf\n",viriel[ip].m11,viriel[ip].m21,viriel[ip].m31);
      fprintf(VIRDAT,"%16.8lf%16.8lf%16.8lf\n",viriel[ip].m12,viriel[ip].m22,viriel[ip].m32);
      fprintf(VIRDAT,"%16.8lf%16.8lf%16.8lf\n",viriel[ip].m13,viriel[ip].m23,viriel[ip].m33);
    }  
    fclose(VIRDAT);
#else 
#ifdef VIRIEL
    FILE *VIRDAT;
    VIRDAT=fopen("TMPDAT/virglob.dat","w");
//    printf("Writing global virial\n");
//    cin.get();
#ifdef EnIneV
    double cfe=1.0/1.03641882007443324881e-04;
    viriel[0]=cfe*viriel[0];
#endif    
    fprintf(VIRDAT,"%18.10lf%18.10lf%18.10lf\n",viriel[0].m11,viriel[0].m21,viriel[0].m31);
    fprintf(VIRDAT,"%18.10lf%18.10lf%18.10lf\n",viriel[0].m12,viriel[0].m22,viriel[0].m32);
    fprintf(VIRDAT,"%18.10lf%18.10lf%18.10lf\n",viriel[0].m13,viriel[0].m23,viriel[0].m33);
    fclose(VIRDAT);
#endif
#endif
  }

//======================================================================//

  inline void printSRT(size_t ip1,size_t ip2,int iout,VCStr& vcs,NBStr& nbs)
  {
    FILE *SRDAT;
    SRDAT=fopen("TMPDAT/srt.dat","w");

    if(iout==0)
    {
      fprintf(SRDAT,"ip1,ip2 : %6zu%6zu\n",ip1,ip2);
//      fprintf(SRDAT,"ip1,ip2 : %6zu%6zu\n",nbs.ipmapinv[ip1],nbs.ipmapinv[ip2]);
      for(size_t ip=ip1;ip<ip2;ip++)
      {
        fprintf(SRDAT,"ip,nnbfl,nnbfr,nnbtot :%8zu%8zu%8zu%8zu\n",ip,nbs.nnbfull[ip]-nbs.nnb0[ip],
//        fprintf(SRDAT,"ip,nnbfl,nnbfr,nnbtot :%8zu%8zu%8zu%8zu\n",nbs.ipmapinv[ip],nbs.nnbfull[ip]-nbs.nnb0[ip],
        nbs.nnbfrac[ip]-nbs.nnbfull[ip],nbs.nnbfrac[ip]-nbs.nnb0[ip]);
        for(size_t inb=nbs.nnb0[ip];inb<nbs.nnbfrac[ip];inb++)
        {
          fprintf(SRDAT,"%8zu",nbs.nbtbl[inb]);
//          fprintf(SRDAT,"%8zu",nbs.ipmapinv[nbs.nbtbl[inb]]);
        }
//        for(size_t inb=nbs.nnb0[ip];inb<nbs.nnbfrac[ip];inb++)
//        {
//          fprintf(SRDAT,"%8zu",gp2cp[nbs.nbtbl[inb]]);
//        }
//        for(size_t inb=nbs.nnb0[ip];inb<nbs.nnbfrac[ip];inb++)
//        {
//          fprintf(SRDAT,"%8zu",  static_cast<size_t>( iflag[nbs.nbtbl[inb]] ) );
//        }
        fprintf(SRDAT,"\n\n");
      }
    }
    else if(iout==1)
    {
      fprintf(SRDAT,"ip1,ip2 : %6zu%6zu\n",ip1,ip2);
      for(size_t ip=ip1;ip<ip2;ip++)
      {
        size_t icl=iclipc(ip,vcs).first;
        size_t ipc=iclipc(ip,vcs).second;
//        if(icl==13&&ipc==165)
//        {
          fprintf(SRDAT,"ip,ispci,icl,ipc,nnbfl,nnbfr,nnbtot :%8zu%6d%8zu%8zu%8zu%8zu%8zu\n",
          ip,vcs.ispc[ip],icl,ipc,nbs.nnbfull[ip]-nbs.nnb0[ip],
          nbs.nnbfrac[ip]-nbs.nnbfull[ip],nbs.nnbfrac[ip]-nbs.nnb0[ip]);
//          fprintf(SRDAT,"statip,esrip : %8zu%18.10lf\n",vcs.istat[ip],nbs.esrpp[ip]);
          for(size_t inb=nbs.nnb0[ip];inb<nbs.nnbfrac[ip];inb++)
          {
            size_t jp=nbs.nbtbl[inb];
            size_t jcl=iclipc(jp,vcs).first;
            size_t jpc=iclipc(jp,vcs).second;
            double rij=nbs.rstor[inb];
            fprintf(SRDAT,"%8zu%8d%8zu%8zu%6d%16.10lf\n",jp,vcs.ispc[jp],jcl,jpc,vcs.istat[jp],rij);
          }
//          fprintf(SRDAT,"\n\n");
//        }
      }
    }
    fclose(SRDAT);
  }

//======================================================================//

  inline size_t ipvcl(size_t ivcl,int ipvc,VCStr& vcs)
  {
    size_t ip0=vcs.ip0vcl[ivcl];
    return ip0+ipvc;
  }

//======================================================================//

  inline pair<size_t,size_t> iclipc(size_t ip,VCStr& vcs)
  {
    size_t ivcl=vcs.ivclp[ip];
    int ipvc=ip-vcs.ip0vcl[ivcl];
    return vcs.pmapinv[ivcl][ipvc];
  }

//======================================================================//

  void printMRT0(size_t ip1,size_t ip2,NBStr& nbs)
  {
    printf("Pure MR bonds\n");
    for(size_t ip=ip1;ip<ip2;ip++)
    {
      printf("ip,nnbmr  :%8zu%8zud\n",ip,nbs.nnbmr[ip]);
      size_t imrij=nbs.indmr[ip];
      for(size_t inb0=0;inb0<nbs.nnbmr[ip];inb0++)
      {
        size_t jp=nbs.nbtblmr[imrij];
        size_t imrji=nbs.indmrji[imrij];
        printf("inb0,ip,jp,imrij,imrji :%6zud%6zud%6zud%6zud%6zud\n",inb0,ip,jp,imrij,imrji);
        imrij=nbs.nextind[imrij];
      }  
    }
  }

//======================================================================//

  void printMRT(size_t ip1,size_t ip2,NBStr& nbs,GhostTbls& gtbls,VCStr& vcs)
  {
    FILE *MRDAT;
//    MRDAT=fopen("TMPDAT/mrt.dat","a");
    MRDAT=fopen("TMPDAT/mrt.dat","w");

    vector<size_t>lstjp;
    vector<size_t>lstmr;
    vector<size_t>lstinbij;
    lstjp.resize(100);
    lstmr.resize(100);
    lstinbij.resize(100);

    for(size_t ip=ip1;ip<ip2;ip++)
    {
      size_t nnbfrmr=nbs.nnbfrac[ip]-nbs.nnbfull[ip];
      for(size_t inb=nbs.nnbfull[ip];inb<nbs.nnbfrac[ip];inb++)
      {
        if(nbs.isr2mr[inb]==0)nnbfrmr--;
      }
      size_t nnbpmr=nbs.nnbmr[ip];
      size_t nnbmr=nnbfrmr+nnbpmr;
//      printf(ip,nnbfrmr,nnbpmr,nnbmr :"%6zu%6zu%6zu%6zu\n",ip,nnbfrmr,nnbpmr,nnbmr);
//      fprintf(MRDAT,"%6zu%6zu%6zu%6zu\n",ip,nnbfrmr,nnbpmr,nnbmr);
      fprintf(MRDAT,"%6zu%6zu%6zu%6zu%6zu\n",iclipc(ip,vcs).first,iclipc(ip,vcs).second,nnbfrmr,nnbpmr,nnbmr);
    }
    fprintf(MRDAT,"\n");

    for(size_t ip=ip1;ip<ip2;ip++)
    {
      int inbnew=0;
      for(size_t inbij=nbs.nnbfull[ip];inbij<nbs.nnbfrac[ip];inbij++)
      {
        size_t jp=nbs.nbtbl[inbij];
        lstjp[inbnew++]=iclipc(jp,vcs).first*1000+iclipc(jp,vcs).second;
      }  
      orderlstjp(inbnew,lstjp,lstinbij);

      for(size_t inb0=nbs.nnbfull[ip];inb0<nbs.nnbfrac[ip];inb0++)
      {
        size_t inbij=nbs.nnbfull[ip]+lstinbij[inb0-nbs.nnbfull[ip]];
//        printf("inb0,inbij: %6zu%6zu\n",inb0,inbij);
//        printf("inb0-nbs.nnbfull[ip]: %6zu%6zu\n",inb0-nbs.nnbfull[ip]);
//        printf("lstinbij[inb0-nbs.nnbfull[ip]]: %6zu%6zu\n",lstinbij[inb0-nbs.nnbfull[ip]]);
        size_t inbji=nbs.indji[inbij];
        size_t jp=nbs.nbtbl[inbij];
//        printf("Frac SR: %6d%6d%6d%6d\n",ip,jp,nbs.isr2mr[inbij],nbs.isr2mr[inbji]);
//        if(jp>ip)continue;
//        if(ip!=171&&jp!=171)continue;
        size_t inblst[2]={inbij,inbji};
        for(size_t ilst=0;ilst<1;ilst++)
        {
          size_t inb=inblst[ilst];
          size_t imr=nbs.isr2mr[inb];
//          printf("Frac SR: %6zu%6zu%6zu\n",ip,jp,imr);
          if(imr>0)
          {
//            fprintf(MRDAT,"Frac MR: %6zu%6zu%6zu\n",ip,jp,gtbls.gp2cp[jp]);
            fprintf(MRDAT,"Frac MR             :%6zu%6zu%6zu%6zu\n",iclipc(ip,vcs).first,iclipc(ip,vcs).second,
            iclipc(jp,vcs).first,iclipc(jp,vcs).second);
            fprintf(MRDAT,"ispci,ispcj         :%6d%6d\n",vcs.ispc[ip],vcs.ispc[jp]);
            fprintf(MRDAT,"rij,smr,snsmr       :%18.10lf%18.10lf%18.10lf\n",nbs.rstor[inbij],nbs.smr[imr],nbs.snsmr[imr]);
            fprintf(MRDAT,"gij,rcmr,drcmrdndb  :%18.10lf%18.10lf%18.10lf\n",nbs.gij[imr],nbs.rcmr[imr],nbs.drcmrdndb[imr]);
            fprintf(MRDAT,"dsnsmrdr,dsnsmrdndb :%18.10lf%18.10lf\n",nbs.dsnsmrdr[imr],nbs.dsnsmrdndb[imr]);
            if(nbs.smr[imr]<1.0)
            {
              fprintf(MRDAT,"dsmrdgij,dgijdndb      :%18.10lf%18.10lf\n",nbs.dsmrdgij[imr],nbs.dgijdndb[imr]);
              fprintf(MRDAT,"dgijdsumk,dgijdnij     :%18.10lf%18.10lf\n",nbs.dgijdsumk[imr],nbs.dgijdnij[imr]);
              fprintf(MRDAT,"dsnsmrdnik,dsnsmrdsumk :%18.10lf%18.10lf\n",nbs.dsnsmrdnik[imr],nbs.dsnsmrdsumk[imr]);
            }
            fprintf(MRDAT,"\n");
//            fprintf(MRDAT,"Force contributions due to Ndb\n");
//            for(size_t ider=nbs.nderndb[imr-1];ider<nbs.nderndb[imr];ider++)
//            {
//              size_t kp=nbs.nderatndb[ider];
//              Vec3d dfr=nbs.frcndb[ider];
//              fprintf(MRDAT,"kp,frcndb :%6zu%16.10lf%16.10lf%16.10lf\n",kp,dfr.x,dfr.y,dfr.z);
//            } 
//            fprintf(MRDAT,"Force contributions due to sumk\n");
//            for(size_t ider=nbs.ndersumk[imr-1];ider<nbs.ndersumk[imr];ider++)
//            {
//              size_t kp=nbs.nderatsumk[ider];
//              Vec3d dfr=nbs.frcsumk[ider];
//              fprintf(MRDAT,"kp,frcsumk :%6zu%16.10lf%16.10lf%16.10lf\n",kp,dfr.x,dfr.y,dfr.z);
//            } 
          }
        }
      }

      inbnew=0;
      size_t imrij=nbs.indmr[ip];
      for(size_t inb0=0;inb0<nbs.nnbmr[ip];inb0++)
      {
        imrij=nbs.nextind[imrij];
        size_t jp=nbs.nbtblmr[imrij];
//        printf("Pure MR here: %6zu%6zu%6zu\n",imrij,ip,jp);
        lstjp[inbnew]=iclipc(jp,vcs).first*1000+iclipc(jp,vcs).second;
        lstmr[inbnew++]=imrij;
      }  
      orderlstjp(inbnew,lstjp,lstinbij);

//      size_t imrij=nbs.indmr[ip];
      for(size_t inb0=0;inb0<nbs.nnbmr[ip];inb0++)
      {
        size_t inbij=lstinbij[inb0];
        size_t imrij=lstmr[inbij];
//        imrij=nbs.nextind[imrij];
        size_t jp=nbs.nbtblmr[imrij];
        size_t imrji=nbs.indmrji[imrij];
//        printf("Pure MR here: %6d%6d%6d%6d%6d\n",inb0,ip,jp,imrij,imrji);
//        if(jp>ip)continue;
//        fprintf(MRDAT,"Pure MR here: %6d%6d%6d%6d%6d\n",inb0,ip,jp,imrij,imrji);
        size_t imrlst[2]={imrij,imrji};
        for(size_t ilst=0;ilst<1;ilst++)
        {
          size_t imr=imrlst[ilst];
          if(imr>0)
          {
//            fprintf(MRDAT,"Pure MR: %6zu%6zu%6zu\n",ip,jp,gtbls.gp2cp[jp]);
            fprintf(MRDAT,"Pure MR             :%6zu%6zu%6zu%6zu\n",iclipc(ip,vcs).first,iclipc(ip,vcs).second,
            iclipc(jp,vcs).first,iclipc(jp,vcs).second);
            fprintf(MRDAT,"ispci,ispcj         :%6d%6d\n",vcs.ispc[ip],vcs.ispc[jp]);
            fprintf(MRDAT,"rij,smr,snsmr       :%18.10lf%18.10lf%18.10lf\n",nbs.rstormr[imr],nbs.smr[imr],nbs.snsmr[imr]);
            fprintf(MRDAT,"gij,rcmr,drcmrdndb  :%18.10lf%18.10lf%18.10lf\n",nbs.gij[imr],nbs.rcmr[imr],nbs.drcmrdndb[imr]);
            fprintf(MRDAT,"dsnsmrdr,dsnsmrdndb :%18.10lf%18.10lf\n",nbs.dsnsmrdr[imr],nbs.dsnsmrdndb[imr]);
            if(nbs.smr[imr]<1.0)
            {
              fprintf(MRDAT,"dsmrdgij,dgijdndb      :%18.10lf%18.10lf\n",nbs.dsmrdgij[imr],nbs.dgijdndb[imr]);
              fprintf(MRDAT,"dgijdsumk,dgijdnij     :%18.10lf%18.10lf\n",nbs.dgijdsumk[imr],nbs.dgijdnij[imr]);
              fprintf(MRDAT,"dsnsmrdnik,dsnsmrdsumk :%18.10lf%18.10lf\n",nbs.dsnsmrdnik[imr],nbs.dsnsmrdsumk[imr]);
            }
            fprintf(MRDAT,"\n");
//            fprintf(MRDAT,"Force contributions due to Ndb\n");
//            for(size_t ider=nbs.nderndb[imr-1];ider<nbs.nderndb[imr];ider++)
//            {
//              size_t kp=nbs.nderatndb[ider];
//              Vec3d dfr=nbs.frcndb[ider];
//              fprintf(MRDAT,"kp,frcndb :%6zu%16.10lf%16.10lf%16.10lf\n",kp,dfr.x,dfr.y,dfr.z);
//            } 
//            fprintf(MRDAT,"Force contributions due to sumk\n");
//            for(size_t ider=nbs.ndersumk[imr-1];ider<nbs.ndersumk[imr];ider++)
//            {
//              size_t kp=nbs.nderatsumk[ider];
//              Vec3d dfr=nbs.frcsumk[ider];
//              fprintf(MRDAT,"kp,frcsumk :%6zu%16.10lf%16.10lf%16.10lf\n",kp,dfr.x,dfr.y,dfr.z);
//            } 
          }
        }
      }
    }
    fclose(MRDAT);
  }

//======================================================================//

    void orderlstjp(int inbnew,vector<size_t>& lstjp,vector<size_t>& lstinbij)
    {
      for(int inb=0;inb<inbnew;inb++){lstinbij[inb]=inb;}  

//      for(size_t inb=0;inb<inbnew;inb++)
//      {
//        size_t inbij=lstinbij[inb];
//        printf("inb,inbij,lstjp :%8zu%8zu%8zu\n",inb,inbij,lstjp[inbij]);
//      }  

      if(inbnew>0)
      {
        for(int i1=0;i1<inbnew-1;i1++)
        for(int i2=i1+1;i2<inbnew;i2++)
        {
          int inb1=lstinbij[i1];
          int inb2=lstinbij[i2];
          if(lstjp[inb2]<lstjp[inb1]){lstinbij[i1]=inb2;lstinbij[i2]=inb1;}
        }  
      }  

//      for(size_t inb=0;inb<inbnew;inb++)
//      {
//        size_t inbij=lstinbij[inb];
//        printf("inb,inbij,lstjp :%8zu%8zu%8zu\n",inb,inbij,lstjp[inbij]);
//      }  
//      abort();

    }

//======================================================================//

    template<typename GridT,typename CellT>
    void printCellDat(GridT& grid,CellT& cells,size_t nclcb,
    vector<size_t>& icllst,vector<size_t>& ip0cl)
    {
      FILE *CellCn;
      CellCn=fopen("TMPDAT/cell.dat","w");
      vector<double*> pri;
      pri.resize(3);  

      for(size_t ilst=0;ilst<nclcb;ilst++)
      {
        size_t icl=icllst[ilst];
        size_t ip0=ip0cl[icl];
        size_t npcli=cells[icl].size();
        if(npcli>0)
        {
//          printf("icl,npcli :%8d%8d\n",icl,npcli);
          fprintf(CellCn,"%8d%8d :",icl,npcli);
          pri={cells[icl][field::rx],cells[icl][field::ry],cells[icl][field::rz]};
          for(size_t ipc=0;ipc<npcli;ipc++) 
          {
            size_t ip=ip0+ipc;
            const Vec3d ri={pri[0][ipc],pri[1][ipc],pri[2][ipc]};
//            printf("%8d",ip);
            fprintf(CellCn,"%8d",ip);
          }
          fprintf(CellCn,"\n");
        }
      }  

//  Verification of atomic positions of ghost atoms
//  for cells or atoms outside central box
//      for(ssize_t igl=ngl-1;igl>=0;igl--)
//      {
//        size_t icnt=0;
//        ssize_t idcli=0;
//        for(ssize_t iclk=igl;iclk<dims.k-igl;iclk++)
//        for(ssize_t iclj=igl;iclj<dims.j-igl;iclj++)
//        {
//          if(iclk==igl||iclk==dims.k-igl-1){idcli=1;}
//          else if(iclj==igl||iclj==dims.j-igl-1){idcli=1;}
//          else{idcli=dims.i-2*igl-1;}
//          for(ssize_t icli=igl;icli<dims.i-igl;icli+=idcli)
//          {
//            IJK iclijk{icli,iclj,iclk};
////            cout << "iclijk =" << iclijk << endl;
//            icl=grid_ijk_to_index(dims,iclijk);
//            ip0=ip0cl[icl];
//            npcli=cells[icl].size();
//            if(npcli>0)
//            {
//              pri={cells[icl][field::rx],cells[icl][field::ry],cells[icl][field::rz]};
//              for(ipc=0;ipc<npcli;ipc++)
//              {
//                ip=ip0+ipc;
//                const Vec3d ri={pri[0][ipc],pri[1][ipc],pri[2][ipc]};
//                printf("ip,ri:%8d%14.8lf%14.8lf%14.8lf\n",ip+1,ri.x,ri.y,ri.z);
//              }
//            }  
//          }  
//        }  
//      }  
//      abort();
    fclose(CellCn);
  }

//======================================================================//

    void wrtatpos(VCStr& vcs)
    {
      FILE *ATPOS;
      ATPOS=fopen("TMPDAT/atpos.xyz","w");
//      ATPOS=fopen("TMPDAT/atpos.xyz","a");

//    char atspc[2][2]={"C ","H "}; 

    size_t ip=0;
    for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
    {
//      printf("wrtatpos; ilst :%8zu\n",ilst);
      size_t ivcl=vcs.icllst[ilst].first;
      int npcli=vcs.pmapinv[ivcl].size();
      for(int ipvc=0;ipvc<npcli;ipvc++,ip++)
      {  
//        printf("wrtatpos; ipvc,ip :%6d%8zu\n",ipvc,ip);
        Vec3d ri=vcs.atpos[ivcl][ipvc];
//        Vec3d ri=vcs.atpos[ivcl][ipvc]+vcs.xyz0;
        if(vcs.ispc[ip]==0)
        {
//          fprintf(ATPOS,"C %18.12lf%18.12lf%18.12lf\n",ri.x,ri.y,ri.z);
          fprintf(ATPOS,"C %6d%18.12lf%18.12lf%18.12lf\n",vcs.istat[ip],ri.x,ri.y,ri.z);
        }
        else
        {
//          fprintf(ATPOS,"H %18.12lf%18.12lf%18.12lf\n",ri.x,ri.y,ri.z);
          fprintf(ATPOS,"H %6d%18.12lf%18.12lf%18.12lf\n",vcs.istat[ip],ri.x,ri.y,ri.z);
        }
      }  
    }

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
      for(ssize_t iclk=iglsk;iclk<vcs.mvcl.k-iglsk;iclk+=idclk)
      {
        for(ssize_t iclj=iglsj;iclj<vcs.mvcl.j-iglsj;iclj+=idclj)
        {
          if(iclk==iglsk||iclk==vcs.mvcl.k-iglsk-1){idcli=1;}
          else if(iclj==iglsj||iclj==vcs.mvcl.j-iglsj-1){idcli=1;}
          else{idcli=vcs.mvcl.i-2*iglsi-1;}
          for(ssize_t icli=iglsi;icli<vcs.mvcl.i-iglsi;icli+=idcli)
          {
            IJK iclijk{icli,iclj,iclk};
//            IJK iclijk{static_cast<ssize_t>(icli),static_cast<ssize_t>(iclj),static_cast<ssize_t>(iclk)};
            size_t ivcl=grid_ijk_to_index(vcs.mvcl,iclijk);
            assert(ivcl>=0&&ivcl<vcs.nvcl);
            int npcli=vcs.pmapinv[ivcl].size();
            for(int ipvc=0;ipvc<npcli;ipvc++,ip++)
            {
              Vec3d ri=vcs.atpos[ivcl][ipvc];
//              Vec3d ri=vcs.atpos[ivcl][ipvc]+vcs.xyz0;
              if(vcs.ispc[ip]==0)
              {
//                fprintf(ATPOS,"C %18.12lf%18.12lf%18.12lf\n",ri.x,ri.y,ri.z);
                fprintf(ATPOS,"C %6d%18.12lf%18.12lf%18.12lf\n",vcs.istat[ip],ri.x,ri.y,ri.z);
              }
              else
              {
//                fprintf(ATPOS,"H %18.12lf%18.12lf%18.12lf\n",ri.x,ri.y,ri.z);
                fprintf(ATPOS,"H %6d%18.12lf%18.12lf%18.12lf\n",vcs.istat[ip],ri.x,ri.y,ri.z);
              }
            }
          }
        }  
      }  
    }
    fclose(ATPOS);
  }

//======================================================================//
//
//    template<typename GridT,typename CellT>
//    void wrtatpos(GridT& grid,CellT& cells,bool cbpos,size_t nclloc,
//    size_t nploc,vector<size_t>& icllst,vector<size_t>& ip0cl,
//    vector<uint8_t>& ispc,vector<Vec3d>& atpos)
//    {
//      FILE *ATPOS;
//      ATPOS=fopen("TMPDAT/atpos.xyz","w");
////      ATPOS=fopen("TMPDAT/atpos.xyz","a");
//      vector<double*>pri,pvi;
//      pri.resize(3);  
//      pvi.resize(3);  
//
//      size_t ngl=grid.ghost_layers();
//      IJK dims=grid.dimension();
//
//      double acl=grid.cell_size();
//      double xl=dims.i*acl;
//      Vec3d shift={0.0,0.0,0.0};
//      if(cbpos){xl-=2*ngl*acl;}
//      else{shift={4.0*acl,4.0*acl,4.0*acl};}
//
//      double box[3];
//      for(uint8_t i=0;i<3;i++){box[i]=xl;}
//
////      double box[3]={10.4660790,10.4660790,10.4660790};
//
//      fprintf(ATPOS,"%8zu\n",nploc);
//      fprintf(ATPOS,"%20.14lf%20.14lf%20.14lf\n",box[0],box[1],box[2]);
//      for(size_t ip=0;ip<nploc;ip++)
//      {
//        Vec3d ri=atpos[ip]+shift;
//        if(ispc[ip]==0)
//        {
//          fprintf(ATPOS,"C %18.12lf%18.12lf%18.12lf\n",ri.x,ri.y,ri.z);
//        }
//        else
//        {
//          fprintf(ATPOS,"H %18.12lf%18.12lf%18.12lf\n",ri.x,ri.y,ri.z);
//        }
//      }
//
////C Code for writing the velocities
////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
////      for(size_t ilst=0;ilst<216;ilst++)
////      {
////        size_t icl=icllst[ilst];
////        size_t npcli=cells[icl].size();
////        printf("icl,npcli :%8zu%8zu\n",icl,npcli);
////        if(npcli>0)
////        {
////          pvi={cells[icl][field::vx],cells[icl][field::vy],cells[icl][field::vz]};
////          for(size_t ipc=0;ipc<npcli;ipc++) 
////          {
////            Vec3d vi={pvi[0][ipc],pvi[1][ipc],pvi[2][ipc]};
////            fprintf(ATPOS,"C %18.12lf%18.12lf%18.12lf\n",vi.x,vi.y,vi.z);
////          }
////        }
////      }
////CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//
////      double box[3]={30.99428918,30.99428918,30.99428918};
////
////      fprintf(ATPOS,"%8zu\n",npcb);
////      fprintf(ATPOS,"%18.12lf%18.12lf%18.12lf\n",box[0],box[1],box[2]);
////      for(size_t ilst=0;ilst<nclcb;ilst++)
////      {
////        size_t icl=icllst[ilst];
//////        size_t ip0=ip0cl[icl];
////        size_t npcli=cells[icl].size();
//////        printf("icl,npcli :%8d%8d\n",icl,npcli);
//////        fprintf(ATPOS,"icl,npcli :%8d%8d\n",icl,npcli);
////        if(npcli>0)
////        {
//////          printf("icl,npcli :%8zu%8zu\n",icl,npcli);
//////          fprintf(ATPOS,"%8zu%8zu :",icl,npcli);
////          pri={cells[icl][field::rx],cells[icl][field::ry],cells[icl][field::rz]};
////          for(size_t ipc=0;ipc<npcli;ipc++) 
////          {
//////            size_t ip=ip0+ipc;
////            const Vec3d ri={pri[0][ipc],pri[1][ipc],pri[2][ipc]};
//////            printf("%8d",ip);
//////            fprintf(ATPOS,"%8d",ip);
////            fprintf(ATPOS,"C %18.12lf%18.12lf%18.12lf\n",ri.x,ri.y,ri.z);
//////            fprintf(ATPOS,"%8zu%16.10lf%16.10lf%16.10lf\n",ip,ri.x,ri.y,ri.z);
////          }
//////          printf("\n");
//////          fprintf(ATPOS,"\n");
////        }
////      }
////
//    fclose(ATPOS);
//  }
//
//======================================================================//
     
  template<typename GridT,typename CellT>
  inline void chkcells(GridT& grid,CellT& cells,size_t nclcb,size_t npcb,size_t np,
  vector<size_t>& icllst,vector<size_t>& ip0cl)
  {
    double acl=grid.cell_size();
    double rclmin=-(acl-2.2)/2.0;
    double rclmax=acl+(acl-2.2)/2.0;
    vector<double*>prgp,pri;
    pri.resize(3);  
    IJK dims=grid.dimension();

    for(size_t ilst=0;ilst<nclcb;ilst++)
    {
      size_t icl=icllst[ilst];
      size_t npcli=cells[icl].size();
      if(npcli==0)continue; 
      IJK ijkcl=grid_index_to_ijk(dims,static_cast<ssize_t>(icl));
      assert( grid.contains(ijkcl) );
      Vec3d rcl=grid.cell_position(ijkcl);
      size_t ip0=ip0cl[icl];
      pri={cells[icl][field::rx],cells[icl][field::ry],cells[icl][field::rz]};
    
      for(size_t ipc=0;ipc<npcli;ipc++) 
      {
        size_t ip=ip0+ipc;
	      assert( ipc>=0 && ipc<cells[icl].size() );
        double ri[3]={pri[0][ipc],pri[1][ipc],pri[2][ipc]};
        //assert( is_inside_threshold( grid.cell_bounds(ijkcl) , Vec3d{ri[0],ri[1],ri[2]} , grid.epsilon_cell_size2() ) );

//C Check whether ri is not too far outside cell
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
        double rimrcl[3]={ri[0]-rcl.x,ri[1]-rcl.y,ri[2]-rcl.z};  
        for(size_t i=0;i<3;i++)
        {
          if(rimrcl[i]<rclmin)
          {
            printf("Atom too far outside cell\n");
            printf("icl,i,ip,rclmin :%6zu%6zu%6zu%14.8lf\n",icl,i,ip,rclmin);
            printf("acl,rimrcl_i :%14.8lf%14.8lf\n",acl,rimrcl[i]);
            abort();
          }  
          else if(rimrcl[i]>rclmax)
          {
            printf("Atom too far outside cell\n");
            printf("icl,i,ip,rclmax :%6zu%6zu%6zu%14.8lf\n",icl,i,ip,rclmax);
            printf("acl,rimrcl_i :%14.8lf%14.8lf\n",acl,rimrcl[i]);
            abort();
          }  
        }
      }  
    }  
  }  

//======================================================================//
     
  template<typename GridT>
  inline void mklstgp2cp(GridT& grid,VCStr& vcs,GhostTbls& gtbls)
  {
    bool cpfound=false;
    Vec3d domsize=grid.dimension()*grid.cell_size();
    Vec3d boxxyz=domsize-2.0*grid.ghost_layers()*grid.cell_size();
    double box[3]={boxxyz.x,boxxyz.y,boxxyz.z};
//    Vec3d xyz0=grid.cell_position(vcs.ijkcb);
//    double acl=grid.cell_size();
////    double rclmin=0.0;
////    double rclmax=acl;
//    Vec3d rc{2.2/hmat.m11,2.2/hmat.m22,2.2/hmat.m33};
//    Vec3d drc{(acl-rc.x)/2.0,(acl-rc.y)/2.0,(acl-rc.z)/2.0};
//    double rclmin[3]={-drc.x,-drc.y,-drc.z};
//    double rclmax[3]={acl+drc.x,acl+drc.y,acl+drc.z};
//    double drclmin=-(acl-2.2)/2.0;
//    double drclmax=acl+(acl-2.2)/2.0;

//    printf("acl :%16.10lf\n",acl);
    printf("domsize :%16.10lf%16.10lf%16.10lf\n",domsize.x,domsize.y,domsize.z);
    printf("boxsize :%16.10lf%16.10lf%16.10lf\n",boxxyz.x,boxxyz.y,boxxyz.z);
//    printf("box     :%22.16lf%22.16lf%22.16lf\n",box[0],box[1],box[2]);
//    printf("rclmin  :%22.16lf%22.16lf%22.16lf\n",rclmin[0],rclmin[1],rclmin[2]);
//    printf("rclmax  :%22.16lf%22.16lf%22.16lf\n",rclmax[0],rclmax[1],rclmax[2]);
//    printf("drclmin  :%22.16lf\n",drclmin);
//    printf("drclmax  :%22.16lf\n",drclmax);

//    double boxsize=1.3283266791429;
//    double cboxsize=1.04660790;
//    box[0]=cboxsize;
//    box[1]=cboxsize;
//    box[2]=cboxsize;
//    printf("box     :%22.16lf%22.16lf%22.16lf\n",box[0],box[1],box[2]);
//    cin.get();

    for(size_t igcl=0;igcl<vcs.nvcl;igcl++)
    {
      int npgcli=vcs.pmapinv[igcl].size();
      if(npgcli==0)continue; 

      double drgp[3]={0.0,0.0,0.0};
      IJK ijkgcl=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(igcl));
      if(ijkgcl.i<vcs.mvclg.i){drgp[0]=box[0];}
      else if(ijkgcl.i>vcs.mvclg.i+vcs.mvclcb.i-1){drgp[0]=-box[0];} 
      if(ijkgcl.j<vcs.mvclg.j){drgp[1]=box[1];}
      else if(ijkgcl.j>vcs.mvclg.j+vcs.mvclcb.j-1){drgp[1]=-box[1];} 
      if(ijkgcl.k<vcs.mvclg.k){drgp[2]=box[2];}
      else if(ijkgcl.k>vcs.mvclg.k+vcs.mvclcb.k-1){drgp[2]=-box[2];} 

      size_t igp0=vcs.ip0vcl[igcl];
      for(int igpc=0;igpc<npgcli;igpc++) 
      {
        size_t igp=igp0+igpc;
        double rgp[3]={vcs.atpos[igcl][igpc].x,vcs.atpos[igcl][igpc].y,vcs.atpos[igcl][igpc].z};
        rgp[0]+=vcs.xyz0.x;rgp[1]+=vcs.xyz0.y;rgp[2]+=vcs.xyz0.z;
        for(size_t i=0;i<3;i++)rgp[i]+=drgp[i];

        cpfound=false;
        for(size_t ilst=0;ilst<vcs.nvclcb;ilst++)
        {
          size_t icl=vcs.icllst[ilst].first;
          int npcli=vcs.icllst[ilst].second;
          if(npcli==0)continue; 
//          IJK ijkcl=grid_index_to_ijk(vcs.mvcl,static_cast<ssize_t>(icl));
//          Vec3d rcl=grid.cell_position(ijkcl);
          size_t ip0=vcs.ip0vcl[icl];
       
          for(int ipc=0;ipc<npcli;ipc++) 
          {
            size_t ip=ip0+ipc;
            double ri[3]={vcs.atpos[icl][ipc].x,vcs.atpos[icl][ipc].y,vcs.atpos[icl][ipc].z};
            ri[0]+=vcs.xyz0.x;ri[1]+=vcs.xyz0.y;ri[2]+=vcs.xyz0.z;
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
//C Check whether ri is not too far outside cell
//            double rimrcl[3]={ri[0]-rcl.x,ri[1]-rcl.y,ri[2]-rcl.z};  
//            for(size_t i=0;i<3;i++)
//            {
//              if(rimrcl[i]<rclmin[i])
//              {
//                printf("Atom too far outside cell\n");
//                printf("icl,i,ip,rclmin :%6zu%6zu%6zu%14.8lf\n",icl,i,ip,rclmin[i]);
//                printf("acl,rimrcl_i :%14.8lf%14.8lf\n",acl,rimrcl[i]);
//                return;
//              }  
//              else if(rimrcl[i]>rclmax[i])
//              {
//                printf("Atom too far outside cell\n");
//                printf("icl,i,ip,rclmax :%6zu%6zu%6zu%14.8lf\n",icl,i,ip,rclmax[i]);
//                printf("acl,rimrcl_i :%14.8lf%14.8lf\n",acl,rimrcl[i]);
//                return;
//              }  
//            }
//CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC            
            double dr[3]={rgp[0]-ri[0],rgp[1]-ri[1],rgp[2]-ri[2]};
            if(abs(dr[0])<1e-9&&abs(dr[1])<1e-9&&abs(dr[2])<1e-9)
            {
              cpfound=true;
              gtbls.gp2cp[igp]=ip;
//              printf("Central image found; \n");
              break;
            }  
          }  
          if(cpfound)break;
        }  
        if(cpfound==false)
        {
          printf("A) Central image not found; igcl,igpc,igp,rgp :%6zu%6d%6zu%16.10lf%16.10lf%16.10lf\n",
          igcl,igpc,igp,rgp[0],rgp[1],rgp[2]);
          printf("box  :%16.10lf%16.10lf%16.10lf\n",box[0],box[1],box[2]);
          printf("cell :%6zd%6zd%6zd\n",ijkgcl.i,ijkgcl.j,ijkgcl.k);
          printf("drgp :%16.10lf%16.10lf%16.10lf\n",drgp[0],drgp[1],drgp[2]);
          return;
        }
      }
    }  
  }  

//======================================================================//
     
  inline void chknbgps(size_t npcb,NBStr& nbs,GhostTbls& gtbls)
  {
    for(size_t ilst1=0;ilst1<gtbls.nglst1;ilst1++)
    {
      size_t igp1=gtbls.glst1[ilst1];
      for(size_t ilst2=0;ilst2<gtbls.nglst2;ilst2++)
      {
        size_t igp2=gtbls.glst2[ilst2];
        if(igp2==igp1||igp1<npcb||igp2<npcb)
        {
          printf("Error; ilst1,ilst2,igp1,igp2 :%6zu%6zu%6zu%6zu\n",ilst1,ilst2,igp1,igp2);
          abort();
        }
      }
    }  

    for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
    {
      size_t igp=gtbls.glst1[ilst];
      chknbigp(igp,nbs,gtbls);
    }  
    printf("Neighbors ghost list 1 are OK\n");

    for(size_t ilst=0;ilst<gtbls.nglst2;ilst++)
    {
      size_t igp=gtbls.glst2[ilst];
      chknbigp(igp,nbs,gtbls);
    }  
    printf("Neighbors ghost list 2 are OK\n");
  }  

//======================================================================//
     
  inline void chknbigp(size_t igp,NBStr& nbs,GhostTbls& gtbls)
  {
    size_t ip=gtbls.gp2cp[igp];
    size_t nnbflgp=nbs.nnbfull[igp]-nbs.nnb0[igp];
    size_t nnbfrgp=nbs.nnbfrac[igp]-nbs.nnbfull[igp];
    //size_t nnbigp=nbs.nnbfrac[igp]-nbs.nnb0[igp];
    size_t nnbfl=nbs.nnbfull[ip]-nbs.nnb0[ip];
    size_t nnbfr=nbs.nnbfrac[ip]-nbs.nnbfull[ip];
    size_t nnbip=nbs.nnbfrac[ip]-nbs.nnb0[ip];
    uint8_t dnnbfl=nnbflgp-nnbfl;
    uint8_t dnnbfr=nnbfrgp-nnbfr;
    uint8_t iok=1;
    if(abs(dnnbfl)>0||abs(dnnbfr)>0)
    {
      printf("A) Ghost neighbor list 1 wrong \n");
      printf("igp,nnbflgp,nnbfrgp %6zu%6zu%6zu\n",igp,nnbflgp,nnbfrgp);
      printf("ip,nnbfl,nnbfr      %6zu%6zu%6zu\n",ip,nnbfl,nnbfr);
      iok=0;
//      abort();
    }
//    printf("igp,nnbflgp,nnbfrgp,nnbigp :%6zu%6zu%6zu%6zu\n",igp,nnbflgp,nnbfrgp,nnbigp);
//    printf("ip,nnbfl,nnbfr,nnbip       :%6zu%6zu%6zu%6zu\n",igp,nnbfl,nnbfr,nnbip);
    double tests1=nbs.crd[igp]-nbs.crd[ip];
    double tests2=nbs.crdspc[igp][0]-nbs.crdspc[ip][0];
    double tests3=nbs.crdspc[igp][1]-nbs.crdspc[ip][1];
    double sumabs=abs(tests1)+abs(tests2)+abs(tests3);
    if(sumabs>1e-10)
    {
      printf("B) Ghost neighbor list 1 wrong \n");
      printf("igp,ip  :%6zu%6zu\n",igp,ip);
      printf("test1,crd    :%16.10lf%16.10lf%16.10lf\n",tests1,nbs.crd[igp],nbs.crd[ip]);
      printf("test2,crdspc :%16.10lf%16.10lf%16.10lf\n",tests2,nbs.crdspc[igp][0],nbs.crdspc[ip][0]);
      printf("test3,crdspc :%16.10lf%16.10lf%16.10lf\n",tests3,nbs.crdspc[igp][1],nbs.crdspc[ip][1]);
      abort();
    }  
    iok=0;
    if(iok==0)
    {
      uint8_t iokall=0;
      for(size_t inb=nbs.nnb0[igp];inb<nbs.nnbfrac[igp];inb++)
      {
        size_t jgp=nbs.nbtbl[inb];
        size_t jp=gtbls.gp2cp[jgp];
//        printf("igp,jgp,jp  %6zu%6zu%6zu\n",igp,jgp,jp);
        uint8_t iok=0;
        for(size_t inbp=nbs.nnb0[ip];inbp<nbs.nnbfrac[ip];inbp++)
        {
          size_t jgp=nbs.nbtbl[inbp];
          size_t jpp=gtbls.gp2cp[jgp];
//          printf("ip,jgp,jpp   %6zu%6zu%6zu\n",ip,jgp,jpp);
          if(jpp==jp)
          {
            iok++;
            iokall++;
            double testr=nbs.rstor[inbp]-nbs.rstor[inb];
            Vec3d dsig=nbs.sigstor[inbp]-nbs.sigstor[inb];
            double testsig=sqrt(norm2(dsig));
            double tests0=nbs.snstor[inbp][0]-nbs.snstor[inb][0];
            double tests1=nbs.snstor[inbp][1]-nbs.snstor[inb][1];
            double tests2=nbs.snstor[inbp][2]-nbs.snstor[inb][2];
            sumabs=abs(testr)+abs(testsig)+abs(tests0)+abs(tests1)+abs(tests2);
            if(sumabs>1e-10)
            {
              printf("C) Ghost neighbor list 1 wrong \n");
              printf("igp,ip  :%6zu%6zu\n",igp,ip);
              printf("jgp,jp  :%6zu%6zu\n",jgp,jp);
              printf("testr   :%16.10lf\n",testr);
              printf("testsig :%16.10lf\n",testsig);
              printf("test0   :%16.10lf\n",tests0);
              printf("test1   :%16.10lf\n",tests1);
              printf("test2   :%16.10lf\n",tests2);
              abort();
            }  
          }  
        }
        if(iok==0)
        {
          printf("Error 1: neighbor not found \n");
          printf("igp,ip,jgp,jp  :%6zu%6zu%6zu%6zu\n",igp,ip,jgp,jp);
        }
      }
      if(iokall-nnbip!=0)
      {
        printf("Error 2: one or more neighbors not found \n");
        printf("iokall,nnbip : %6zu%6zu\n",static_cast<size_t>(iokall),static_cast<size_t>(nnbip));
        abort();
      }
    }
  }  

//======================================================================//
     
  inline void chkmrgps(VLCHBOP& par,vector<uint8_t>& ispc,
  vector<int8_t>& istat,NBStr& nbs,GhostTbls& gtbls)
  {
    for(size_t ilst=0;ilst<gtbls.nglst1;ilst++)
    {
      size_t igp=gtbls.glst1[ilst];
//      printf("ilst,igp,iflagigp :%6zu%6zu%6d\n",ilst,igp,iflag[igp]);
      if(istat[igp]<1)continue;
//      if(iflag[igp]<1)continue;
      chkmrigp(par,igp,ispc,nbs,gtbls);
    }  
    printf("MR tables ghost list 1 are OK\n");

    for(size_t ilst=0;ilst<gtbls.nglst2;ilst++)
    {
      size_t igp=gtbls.glst2[ilst];
      if(istat[igp]<1)continue;
//      if(iflag[igp]<1)continue;
      chkmrigp(par,igp,ispc,nbs,gtbls);
    }  
    printf("MR tables ghost list 2 are OK\n");
  }  

//======================================================================//
     
  inline void chkmrigp(VLCHBOP& par,size_t igp,vector<uint8_t>& ispc,
  NBStr& nbs,GhostTbls& gtbls)
  {
    size_t ip=gtbls.gp2cp[igp];
    size_t nnbmrigp=nbs.nnbmr[igp];
    size_t nnbmrip=nbs.nnbmr[ip];
    size_t dnnbmr=nnbmrigp-nnbmrip;
    double ndbigp=nbs.ndb[igp];
    double ndbip=nbs.ndb[ip];
    double testndb=ndbigp-ndbip;
    if(abs(testndb)>1e-10)
    {
      printf("Dangling bond number wrong\n");
      printf("igp,ip,ndbs :%6zu%6zu%16.10lf%16.10lf\n",igp,ip,ndbigp,ndbip);
      abort();
    }  
    if(labs(dnnbmr)>0)
    {
      printf("A) MR tables of ghost list 1 possible wrong \n");
      printf("igp,ip,nnbmrigp,nnbmrip : %6zu%6zu%6zu%6zu\n",igp,ip,nnbmrigp,nnbmrip);
//      size_t iok=0;

      int ispcc=0;
      int ispci=ispc[ip];
//      uint8_t ispcc=0;
//      uint8_t ispci=nbs.ispc[ip];
      double ndbij=nbs.ndb[ip];
      double rcmr=0.;
      double drcmrdndb=0.;
      size_t imrij=nbs.indmr[igp];
      for(size_t inb=0;inb<nbs.nnbmr[igp];inb++)
      {
//        printf("imrij,nextimrij: %6d%6d\n",imrij,nbs.nextind[imrij]);
        imrij=nbs.nextind[imrij];
        size_t jp=nbs.nbtblmr[imrij];
//        if(jp<216)
//        {
          printf("Neighbor ghost particle; jp = %6zu\n",jp);
//          continue;
//        }  
        double rij=nbs.rstormr[imrij];
        int ispcj=ispc[jp];
//        uint8_t ispcj=nbs.ispc[jp];
        par.pRC[ispci][ispcj].rcmrij(ispci,ispcj,ispcc,ndbij,rcmr,drcmrdndb);
        if(rij>=rcmr)nnbmrigp--;
      }  
      imrij=nbs.indmr[ip];
      for(size_t inb=0;inb<nbs.nnbmr[ip];inb++)
      {
//        printf("imrij,nextimrij: %6d%6d\n",imrij,nbs.nextind[imrij]);
        imrij=nbs.nextind[imrij];
        size_t jp=nbs.nbtblmr[imrij];
//        if(jp<216)
//        {
          printf("Neighbor central box particle; jp = %6zu\n",jp);
//          continue;
//        }  
        double rij=nbs.rstormr[imrij];
        int ispcj=ispc[jp];
//        uint8_t ispcj=nbs.ispc[jp];
        par.pRC[ispci][ispcj].rcmrij(ispci,ispcj,ispcc,ndbij,rcmr,drcmrdndb);
        if(rij>=rcmr)nnbmrip--;
      }  
      size_t dnnbmr=nnbmrigp-nnbmrip;
      if(labs(dnnbmr)>0)
      {
        printf("B) MR tables of ghost list 1 definitively wrong \n");
        printf("igp,ip,nnbmrigp,nnbmrip : %6zu%6zu%6zu%6zu\n",igp,ip,nnbmrigp,nnbmrip);
        abort();
      }  
    }
////    printf("igp,nnbflgp,nnbfrgp,nnbigp :%6zu%6zu%6zu%6zu\n",igp,nnbflgp,nnbfrgp,nnbigp);
////    printf("ip,nnbfl,nnbfr,nnbip       :%6zu%6zu%6zu%6zu\n",igp,nnbfl,nnbfr,nnbip);
//    iok=0;
//    if(iok==0)
//    {
//      uint8_t iokall=0;
//      for(size_t inb=nbs.nnb0[igp];inb<nbs.nnbfrac[igp];inb++)
//      {
//        size_t jgp=nbs.nbtbl[inb];
//        size_t jp=gtbls.gp2cp[jgp];
////        printf("igp,jgp,jp  %6zu%6zu%6zu\n",igp,jgp,jp);
//        uint8_t iok=0;
//        for(size_t inbp=nbs.nnb0[ip];inbp<nbs.nnbfrac[ip];inbp++)
//        {
//          size_t jgp=nbs.nbtbl[inbp];
//          size_t jpp=gtbls.gp2cp[jgp];
////          printf("ip,jgp,jpp   %6zu%6zu%6zu\n",ip,jgp,jpp);
//          if(jpp==jp)
//          {
//            iok++;
//            iokall++;
//            double testr=nbs.rstor[inbp]-nbs.rstor[inb];
//            Vec3d dsig=nbs.sigstor[inbp]-nbs.sigstor[inb];
//            double testsig=sqrt(norm2(dsig));
//            double tests0=nbs.snstor[inbp][0]-nbs.snstor[inb][0];
//            double tests1=nbs.snstor[inbp][1]-nbs.snstor[inb][1];
//            double tests2=nbs.snstor[inbp][2]-nbs.snstor[inb][2];
//            double tests3=nbs.crd[igp]-nbs.crd[ip];
//            double tests4=nbs.crdspc[igp][0]-nbs.crdspc[ip][0];
//            double tests5=nbs.crdspc[igp][1]-nbs.crdspc[ip][1];
//            double sumabs=abs(testr)+abs(testsig)+abs(tests0)+abs(tests1)+
//            abs(tests2)+abs(tests3)+abs(tests4)+abs(tests5);
//            if(sumabs>1e-10)
//            {
//              printf("B) Ghost neighbor list 1 wrong \n");
//              printf("testr   :%16.10lf\n",testr);
//              printf("testsig :%16.10lf\n",testsig);
//              printf("test0   :%16.10lf\n",tests0);
//              printf("test1   :%16.10lf\n",tests1);
//              printf("test2   :%16.10lf\n",tests2);
//              printf("test3   :%16.10lf\n",tests3);
//              printf("test4   :%16.10lf\n",tests4);
//              printf("test5   :%16.10lf\n",tests5);
//              abort();
//            }  
//          }  
//        }
//        if(iok==0)
//        {
//          printf("Error 1: neighbor not found \n");
//          printf("igp,ip,jgp,jp  :%6zu%6zu%6zu%6zu\n",igp,ip,jgp,jp);
//        }
//      }
//      if(iokall-nnbip!=0)
//      {
//        printf("Error 2: one or more neighbors not found \n");
//        printf("iokall,nnbip : %6zu%6zu\n",iokall,nnbip);
//        abort();
//      }
//    }
  }  

//======================================================================//
