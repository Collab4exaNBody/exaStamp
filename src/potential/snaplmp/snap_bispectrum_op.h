#pragma once

namespace exaStamp
{

  using namespace exanb;
//  using onika::memory::DEFAULT_ALIGNMENT;
//  using namespace SnapExt;

  // Force operator
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) BispectrumOp 
  {
    SnapLMPThreadContext* m_thread_ctx = nullptr;
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
    const bool switchinnerflag = false;
    const bool chemflag = false;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () ( int jnum, ComputeBufferT& buf, int type, CellParticlesT cells ) const
    {
      size_t tid = omp_get_thread_num();
      assert( tid < n_thread_ctx );
      SnapLMPThreadContext & snap_ctx = m_thread_ctx[tid];
      auto snaptr = snap_ctx.sna;

      // start of SNA
      // int type=0;
      const int itype = type; //type[i];
      const int ielem = itype; //map[itype];
      const double radi = radelem[ielem];
      // std::cout << "radius here for bispectrum = " << radi << std::endl;

      snaptr->grow_rij(jnum);
      
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
          snaptr->rij[ninside][0] = delx;
          snaptr->rij[ninside][1] = dely;
          snaptr->rij[ninside][2] = delz;
          snaptr->inside[ninside] = j;
          snaptr->wj[ninside] = wjelem[jelem];
          snaptr->rcutij[ninside] = (radi + radelem[jelem])*rcutfac;
          if (switchinnerflag) {
            snaptr->sinnerij[ninside] = 0.5*(sinnerelem[ielem]+sinnerelem[jelem]);
            snaptr->dinnerij[ninside] = 0.5*(dinnerelem[ielem]+dinnerelem[jelem]);
          }
          if (chemflag) snaptr->element[ninside] = jelem;
          ninside++;
        }
      }


      /********** DEBUG ***************/
      // snaptr->reorder_rij(ninside);
      // const double chk_rij = snaptr->chk_rij(ninside);
      const bool chosen = false; //fabs(chk_rij-1.862026)<0.00001;
      // snaptr->set_dbg_chosen(chosen);
      if(chosen)
      {
        printf("SNAPDBG: compute_bispectrum: switchinnerflag=%d, cutsq=% .3e, rcutfac=% .3e, ninside=%d\n",switchinnerflag,cutsq,rcutfac,ninside);
        for(int jj=0;jj<ninside;jj++)
        {
          printf("SNAPDBG:\trij[%d]=(% .3e,% .3e,% .3e) wj=% .3e rcutij=% .3e sinnerij=% .3e dinnerij=% .3e\n",
                 jj,snaptr->rij[jj][0],snaptr->rij[jj][1],snaptr->rij[jj][2] , snaptr->wj[jj], snaptr->rcutij[jj], snaptr->sinnerij[jj], snaptr->dinnerij[jj] );          
        }
      }
      /*******************************/


      if (chemflag) snaptr->compute_ui(ninside, ielem);
      else          snaptr->compute_ui(ninside, 0);

      snaptr->compute_zi();

      if (chemflag) snaptr->compute_bi(ielem);
      else          snaptr->compute_bi(0);

      const long bispectrum_ii_offset = ncoeff * ( cell_particle_offset[buf.cell] + buf.part );
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        bispectrum[ bispectrum_ii_offset + icoeff ] /*[ii][icoeff]*/ = snaptr->blist[icoeff];
      }


      /********** DEBUG ***************/
      if( chosen )
      {
        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
        {
          if(icoeff%4 == 0) printf("\nSNAP//DBG: BI:\t");
          printf("% .6e ", snaptr->blist[icoeff] );
        }
        printf("\n");
      }
      /*******************************/

    }
    
  };

}


