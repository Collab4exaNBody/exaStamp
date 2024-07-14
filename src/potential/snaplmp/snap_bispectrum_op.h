#pragma once

#include "snap_compute_ui.h"
#include "snap_compute_zi.h"
#include "snap_compute_bi.h"

namespace exaStamp
{

  using namespace exanb;
//  using onika::memory::DEFAULT_ALIGNMENT;
//  using namespace SnapExt;

  // Force operator
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) BispectrumOp 
  {
    const LAMMPS_NS::ReadOnlySnapParameters snaconf;
    LAMMPS_NS::SnapScratchBuffers* m_scratch = nullptr;
    const size_t n_thread_ctx = 0;
    
    const size_t * const cell_particle_offset = nullptr;
    const double * const beta = nullptr;
    double * const bispectrum = nullptr;
    
    const double * const coeffelem = nullptr;
    const long ncoeff = 0;
    
    const double * const wjelem = nullptr; // data of m_factor in snap_ctx
    const double * const radelem = nullptr;
    const double * const sinnerelem = nullptr;
    const double * const dinnerelem = nullptr;
    const double rcutfac = 0.0;
    const double cutsq = rcutfac*rcutfac;
    const bool eflag = false;
    const bool quadraticflag = false;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () ( int jnum, ComputeBufferT& buf, int type, CellParticlesT cells ) const
    {
      assert( ONIKA_CU_BLOCK_IDX < n_thread_ctx );
      auto & snabuf = m_scratch[ ONIKA_CU_BLOCK_IDX ];

      // start of SNA
      // int type=0;
      const int itype = type; //type[i];
      const int ielem = itype; //map[itype];
      const double radi = radelem[ielem];
      // std::cout << "radius here for bispectrum = " << radi << std::endl;

      assert( jnum <= snabuf.nmax );
      // snabuf.grow_rij(jnum);
      
      int ninside = 0;
      for (int jj = 0; jj < jnum; jj++)
      {
        const int j = jj;// j = jlist[jj];
        // j &= NEIGHMASK;
        const double delx = buf.drx[jj]; //x[j][0] - xtmp;
        const double dely = buf.dry[jj]; //x[j][1] - ytmp;
        const double delz = buf.drz[jj]; //x[j][2] - ztmp;
        const double rsq = delx*delx + dely*dely + delz*delz;
	      const int jtype = buf.nbh_pt[j][field::type];
	      // const int jtype = 0;
        const int jelem = jtype; //map[jtype];

	      double cut_ij = (radelem[jtype]+radelem[itype])*rcutfac;
	      double cutsq_ij=cut_ij*cut_ij;
	      // cutsq[itype][jtype] = (radelem[itype]+radelem[jtype])*rcutfac;
	      // std::cout << "custq[" << itype << "][" << jtype << "] = " << cutsq_ij << std::endl;

        if (rsq < cutsq_ij/*[itype][jtype]*/&&rsq>1e-20) {
          //snabuf.rij[ninside][0] = delx;
          //snabuf.rij[ninside][1] = dely;
          //snabuf.rij[ninside][2] = delz;
          //snabuf.inside[ninside] = j;
          if( ninside != jj ) { buf.copy( jj , ninside ); }
          
          snabuf.wj[ninside] = wjelem[jelem];
          snabuf.rcutij[ninside] = (radi + radelem[jelem])*rcutfac;
          if (snaconf.switch_inner_flag) {
            snabuf.sinnerij[ninside] = 0.5*(sinnerelem[ielem]+sinnerelem[jelem]);
            snabuf.dinnerij[ninside] = 0.5*(dinnerelem[ielem]+dinnerelem[jelem]);
          }
          if (snaconf.chem_flag) snabuf.element[ninside] = jelem;
          ninside++;
        }
      }

#if 0
      /********** DEBUG ***************/
      // snabuf.reorder_rij(ninside);
      // const double chk_rij = snabuf.chk_rij(ninside);
      const bool chosen = false; //fabs(chk_rij-1.862026)<0.00001;
      // snabuf.set_dbg_chosen(chosen);
      if(chosen)
      {
        printf("SNAPDBG: compute_bispectrum: switchinnerflag=%d, cutsq=% .3e, rcutfac=% .3e, ninside=%d\n",switchinnerflag,cutsq,rcutfac,ninside);
        for(int jj=0;jj<ninside;jj++)
        {
          printf("SNAPDBG:\trij[%d]=(% .3e,% .3e,% .3e) wj=% .3e rcutij=% .3e sinnerij=% .3e dinnerij=% .3e\n",
                 jj,snabuf.rij[jj][0],snabuf.rij[jj][1],snabuf.rij[jj][2] , snabuf.wj[jj], snabuf.rcutij[jj], snabuf.sinnerij[jj], snabuf.dinnerij[jj] );          
        }
      }
      /*******************************/
#endif

      snap_compute_ui( snaconf.nelements, snaconf.twojmax, snaconf.idxu_max, snaconf.idxu_block
                     , snabuf.element, buf.drx,buf.dry,buf.drz , snabuf.rcutij, snaconf.rootpqarray, snabuf.sinnerij, snabuf.dinnerij, snabuf.wj
                     , snaconf.wselfall_flag, snaconf.switch_flag, snaconf.switch_inner_flag, snaconf.chem_flag
                     , snaconf.wself, snaconf.rmin0, snaconf.rfac0
                     , snabuf.ulist_r_ij, snabuf.ulist_i_ij, snabuf.ulisttot_r, snabuf.ulisttot_i
                     , ninside, snaconf.chem_flag ? ielem : 0);

      // snabuf.compute_zi();
      snap_compute_zi( snaconf.nelements, snaconf.idxz_max, snaconf.idxu_max, snaconf.idxu_block, snaconf.idxcg_block
                     , snaconf.idxz, snaconf.cglist, snabuf.ulisttot_r, snabuf.ulisttot_i
                     , snaconf.bnorm_flag, snabuf.zlist_r, snabuf.zlist_i );

//      if (chemflag) snabuf.compute_bi(ielem);
//      else          snabuf.compute_bi(0);
      snap_compute_bi( snaconf.nelements, snaconf.idxz_max, snaconf.idxb_max, snaconf.idxu_max
                     , snaconf.idxu_block, snaconf.idxz_block
                     , snaconf.idxz, snaconf.idxb
                     , snabuf.zlist_r, snabuf.zlist_i
                     , snabuf.ulisttot_r, snabuf.ulisttot_i
                     , snaconf.bzero , snaconf.bzero_flag, snaconf.wselfall_flag
                     , snabuf.blist
                     , snaconf.chem_flag ? ielem : 0 );

      const long bispectrum_ii_offset = ncoeff * ( cell_particle_offset[buf.cell] + buf.part );
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        bispectrum[ bispectrum_ii_offset + icoeff ] /*[ii][icoeff]*/ = snabuf.blist[icoeff];
      }

#if 0
      /********** DEBUG ***************/
      if( chosen )
      {
        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
        {
          if(icoeff%4 == 0) printf("\nSNAP//DBG: BI:\t");
          printf("% .6e ", snabuf.blist[icoeff] );
        }
        printf("\n");
      }
      /*******************************/
#endif

    }
    
  };

}


