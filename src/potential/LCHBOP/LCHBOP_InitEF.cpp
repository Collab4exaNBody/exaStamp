#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/algorithm.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/quantity.h>
#include <exanb/core/particle_id_codec.h>

#include <memory>
//#include <math.h>
#include <iostream>
#include <typeinfo>
#include <iomanip>
#include <type_traits>
//#include <mpi.h>
//#include <sys/time.h>

#include "lchbop_inits.h"
#include "lchbop_getEF.h"
#include "lchbop_utils.h"

using std::vector;
using std::cin;
//using std::clock;

template<typename GridT>
void InitEF(GridT&,vector<int>&,NbCells&,VCStr&);

inline void getstatcl(vector<int>&,VCStr&);

template<typename GridT>
inline void orderatpos(GridT&,VCStr& vcs);


 
namespace exaStamp
{
  template<
	typename GridT
        ,  class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz >
	>
  class LCHBOP_InitEF : public OperatorNode
  {
    using CellOwner  = std::vector<int>;

// ========= I/O slots =======================
//    ADD_SLOT( MPI_Comm  , mpi           , INPUT );
    ADD_SLOT( Domain    , domain     , INPUT , REQUIRED );
    ADD_SLOT( GridT     , grid       , INPUT );
    ADD_SLOT( CellOwner , cell_owner , INPUT , REQUIRED );
    ADD_SLOT( bool      , NPT        , INPUT );
    ADD_SLOT( double    , rcut_inc   , INPUT );
    ADD_SLOT( NbCells   , nbclist    , INPUT );
    ADD_SLOT( VCStr     , vcs        , INPUT_OUTPUT );

    public:
    
    void execute () override final
    {

//      MPI_Comm_rank(*mpi, &vcs.rank);
      GridT& grid = *(this->grid);
      CellOwner& cell_owner = *(this->cell_owner);
      NbCells& nbc = *(this->nbclist);
      VCStr& vcs = *(this->vcs);
      vcs.NPT = *(this->NPT);
      vcs.rcutinc0 = *(this->rcut_inc);
//      cout << "vcs.NPT =" << vcs.NPT << endl;

//C Initialisations for calculation of energy and forces
      InitEF(grid,cell_owner,nbc,vcs);
//      printf("Operator LCHBOP_InitEF done\n");
    }
  };

  template<class GridT> using LCHBOP_InitEF_Tmpl = LCHBOP_InitEF<GridT>;

// === register factories ===  
  __attribute__((constructor)) static void register_factories()
  {
    OperatorNodeFactory::instance()->register_factory( "LCHBOP_InitEF",
    make_grid_variant_operator< LCHBOP_InitEF_Tmpl > );
  }

} // end namespace exaStamp

//=========================================================================//

  template<typename GridT>
  void InitEF(GridT& grid,vector<int>& cell_owner,NbCells& nbc,VCStr& vcs)
  {
// Storing exaStamp cell data in vcs   
    vcs.np=grid.number_of_particles();
    vcs.dims=grid.dimension();
    size_t ngl=grid.ghost_layers();  
    vcs.ngls0=static_cast<ssize_t>(ngl);
    vcs.ncl=grid.number_of_cells();
    vcs.dimscb={vcs.dims.i-2*vcs.ngls0,vcs.dims.j-2*vcs.ngls0,vcs.dims.k-2*vcs.ngls0};
    vcs.ijkcb={vcs.ngls0,vcs.ngls0,vcs.ngls0};
    vcs.decl=grid.cell_size();

// Getting maximal number of atoms in an exaStamp cell
    auto cells=grid.cells();
    size_t npcli,npcm=0;
    for(size_t icl=0;icl<vcs.ncl;icl++)
    {
      npcli=cells[icl].size();
      if(npcli>0){if(npcli>npcm)npcm=npcli;}  
    }
    vcs.initvcs1(npcm);
    getstatcl(cell_owner,vcs);

// Verlet cell structure is still to be made
    vcs.mkvcs=true;
  }

//======================================================================//
     
  inline void getstatcl(vector<int>& cell_owner,VCStr& vcs)
  {
//C Declare central box block 
      GridBlock cbox;
      cbox.start={vcs.ngls0,vcs.ngls0,vcs.ngls0};
      cbox.end={vcs.dims.i-vcs.ngls0,vcs.dims.j-vcs.ngls0,vcs.dims.k-vcs.ngls0};
//      ijkcl={0,0,0};
//      ijkcl={ngls,ngls,ngls};
//      ijkcl={vcs.dims.i-ngls,vcs.dims.k-ngls,vcs.dims.k-ngls};
//      bool incbox=inside_block(cbox,ijkcl);
//      cin.get();

//C Get cell_owner of central box
      size_t icl=grid_ijk_to_index(vcs.dims,cbox.start);
      size_t iclowncb=cell_owner[icl];
//      printf("icl,iclowncb :%8zu%8zu\n",icl,iclowncb);
//      cin.get();
//      vcs.clown[icl]=iclowncb;                  // FOR TESTING ONLY

//C Get status of the exaStamp cells, i.e. central box, active ghost or non-active ghost.
    for(size_t icl=0;icl<vcs.ncl;icl++)
    {
      IJK ijkcl=grid_index_to_ijk(vcs.dims,static_cast<ssize_t>(icl));
      if(inside_block(cbox,ijkcl)){vcs.istatcl[icl]=4;continue;}
      size_t iclown=cell_owner[icl];
//      vcs.clown[icl]=iclown;    // FOR TESTING ONLY
      if(iclown<iclowncb){vcs.istatcl[icl]=-2;}
      else if(iclown>iclowncb){vcs.istatcl[icl]=-3;}
      else
      {
        if(ijkcl.k<vcs.ngls0){vcs.istatcl[icl]=-2;}
        else if(ijkcl.k<vcs.dims.k-vcs.ngls0)
        {
          if(ijkcl.j<vcs.ngls0){vcs.istatcl[icl]=-2;}
          else if(ijkcl.j<vcs.dims.j-vcs.ngls0)
          {
            if(ijkcl.i<vcs.ngls0){vcs.istatcl[icl]=-2;}
            else{vcs.istatcl[icl]=-3;}
          }
          else{vcs.istatcl[icl]=-3;}
        } 
        else{vcs.istatcl[icl]=-3;}
      }
    }

  } // end of routine getstatcl

//======================================================================//

  template<typename GridT>
  void orderatpos(GridT& grid,VCStr& vcs)
  {
    FILE *ATPOS0;
    ATPOS0=fopen("TMPDAT/atpos_0.xyz","r");
    int err;
    char texte[2];
    Vec3d ri;
  
    vector<pair<size_t,size_t>>pmapinvcl;
    pmapinvcl.resize(20);

    size_t ivcl0=vcs.ivclp[0];
    int ipvcnew=0;
    int istat=0;
    int npcli=vcs.pmapinv[ivcl0].size();
//    printf("ivcl,npcli,npvcm :%6zu%6d%6zu\n",ivcl0,npcli,vcs.npvcm);
    for(int ipvc=0;ipvc<npcli;ipvc++)
    {
      vcs.atposcl[ipvc]=vcs.atpos[ivcl0][ipvc];
      pmapinvcl[ipvc]=vcs.pmapinv[ivcl0][ipvc];
    }  

//    for(size_t ip=0;ip<4092;ip++)
    for(size_t ip=0;;ip++)
    {
//      err=fscanf(ATPOS0,"%s%lf%lf%lf\n",texte,&ri.x,&ri.y,&ri.z);	
      err=fscanf(ATPOS0,"%s%d%lf%lf%lf\n",texte,&istat,&ri.x,&ri.y,&ri.z);	
//      printf("ip,err :%6zu%6d\n",ip,err);
//      printf("%2s%6d%18.12lf%18.12lf%18.12lf\n",texte,istat,ri.x,ri.y,ri.z);
      if(err!=5)break;
//      cin.get();
//      if(ip==100){printf("ip,atspc_ip : %6zu%2s\n",ip,texte);cin.get();}

      size_t ivcl=vcs.ivclp[ip];
//      printf("ivcl :%6zu\n",ivcl0);
//      cin.get();
      if(ivcl!=ivcl0)
      {
        ivcl0=ivcl;
        ipvcnew=0;
        npcli=vcs.pmapinv[ivcl].size();
        for(int ipvc=0;ipvc<npcli;ipvc++)
        {
          vcs.atposcl[ipvc]=vcs.atpos[ivcl][ipvc];
          pmapinvcl[ipvc]=vcs.pmapinv[ivcl][ipvc];
        }  
      }
      int ipvc=0;
//      printf("ivcl,ipvcnew,npcli :%6zu%6d%6d\n",ivcl,ipvcnew,npcli);
      Vec3d dri;
      for(;ipvc<npcli;ipvc++)
      {
        Vec3d dri=ri-vcs.atposcl[ipvc];
//        if(ip==100)
//        {
//          printf("ivcl,ipvc,npcli :%6zu%6d%6d%18.12lf%18.12lf%18.12lf\n",
//          ipvc,ipvc,npcli,dri.x,dri.y,dri.z);
//        }  
        if(abs(dri.x)<1e-10&&abs(dri.y)<1e-10&&abs(dri.z)<1e-10)break;
      }
      if(abs(dri.x)>1e-10||abs(dri.y)>1e-10||abs(dri.z)>1e-10)
      {
//      printf("BUG; ipvc :%6d\n",ipvc);
        abort();
      }  
//      cin.get();
      vcs.atpos[ivcl][ipvcnew]=ri;
      size_t icl=pmapinvcl[ipvc].first;
      size_t ipc=pmapinvcl[ipvc].second;
      vcs.pmapinv[ivcl][ipvcnew]=pair<size_t,size_t>(icl,ipc);
      vcs.pmap[icl][ipc]=pair<size_t,int>(ivcl,ipvcnew);
//      if(icl==vcs.iclt&&ipc==vcs.ipct)
//      {
//        int ispci=grid.particle_type(icl,ipc);
//        printf("icl,ipc,ispci :%6zu%6zu%6d\n",icl,ipc,ispci);
//        size_t iclt=vcs.pmapinv[ivcl][ipvcnew].first;
//        size_t ipct=vcs.pmapinv[ivcl][ipvcnew].second;
//        int ispcit=grid.particle_type(iclt,ipct);
//        printf("iclt,ipct,ispcit :%6zu%6zu%6d\n",iclt,ipct,ispcit);
//        printf("ivcl,ipvcnew :%6zu%6d\n",ivcl,ipvcnew);
//        cin.get();
//      }
//      ipvcnew++;
    }  
//    abort();
  }

//======================================================================//
