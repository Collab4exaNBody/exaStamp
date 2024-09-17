#pragma once

#include <exanb/core/physics_constants.h>

namespace exaStamp
{

  using namespace exanb;
//  using onika::memory::DEFAULT_ALIGNMENT;
//  using namespace SnapExt;

  // additional storage space added to compute buffer created by compute_cell_particle_pairs
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) SnapComputeBuffer
  {
    alignas(onika::memory::DEFAULT_ALIGNMENT) int type[exanb::MAX_PARTICLE_NEIGHBORS];
  };

  // functor that populate compute buffer's extended storage for particle charges
  struct CopyParticleType
  {
    template<class ComputeBufferT, class FieldArraysT>
    /* ONIKA_HOST_DEVICE_FUNC */ inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, FieldArraysT cells, size_t cell_b, size_t p_b, double weight) const noexcept
    {
      assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
      tab.ext.type[tab.count] = cells[cell_b][field::type][p_b];
      exanb::DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, weight );
    }
  };

  // Force operator
  struct alignas(onika::memory::DEFAULT_ALIGNMENT) ForceOp 
  {
    SnapLMPThreadContext* m_thread_ctx = nullptr;
    const size_t n_thread_ctx = 0;
    
    const size_t * const cell_particle_offset = nullptr;
    const double * const beta = nullptr;
    const double * const bispectrum = nullptr;
    
    const double * const coeffelem = nullptr;
    const long coeffelem_size = 0;
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
    
    // exaStamp conversion specific falgs
    const bool conv_energy_units = true;
    //    static inline constexpr double conv_energy_factor = 1e-4 * exanb::legacy_constant::elementaryCharge / exanb::legacy_constant::atomicMass;
    double conv_energy_factor = UnityConverterHelper::convert(1., "eV");

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      CellParticlesT cells
      ) const
    {
      FakeMat3d virial;
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      Mat3d& virial,
      CellParticlesT cells
      ) const
    {
      ComputePairOptionalLocks<false> locks;
      FakeParticleLock lock_a;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class GridCellLocksT, class ParticleLockT>
    inline void operator ()
      (
      size_t n,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      CellParticlesT cells,
      GridCellLocksT locks,
      ParticleLockT& lock_a
      ) const
    {
      FakeMat3d virial;
      this->operator () ( n,buf,en,fx,fy,fz,type,virial,cells,locks,lock_a);
    }

    template<class ComputeBufferT, class CellParticlesT, class Mat3dT,class GridCellLocksT, class ParticleLockT>
    inline void operator ()
      (
      int jnum ,
      ComputeBufferT& buf,
      double& en,
      double& fx,
      double& fy,
      double& fz,
      int type,
      Mat3dT& virial ,
      CellParticlesT cells,
      GridCellLocksT locks,
      ParticleLockT& lock_a
      ) const
    {
      static constexpr bool compute_virial = std::is_same_v< Mat3dT , Mat3d >;
      //static constexpr double scale[1][1] = { {1.0} };

      size_t tid = omp_get_thread_num();
      assert( tid < n_thread_ctx );
      SnapLMPThreadContext & snap_ctx = m_thread_ctx[tid];
      auto snaptr = snap_ctx.sna;

      // energy and force contributions to the particle
      Mat3dT _vir; // default constructor defines all elements to 0
      double _en = 0.;
      double _fx = 0.;
      double _fy = 0.;
      double _fz = 0.;

      // start of SNA

      //      const int itype = 0; //type[i];
      const int itype = type;
      const int ielem = itype; //map[itype];
      const double radi = radelem[ielem];

      //      snaptr->grow_rij(jnum);
      assert( jnum < 1024 );
      int ninside = 0;
      for (int jj = 0; jj < jnum; jj++)
      {
        const int j = jj;// j = jlist[jj];
        // j &= NEIGHMASK;
        const double delx = buf.drx[jj]; //x[j][0] - xtmp;
        const double dely = buf.dry[jj]; //x[j][1] - ytmp;
        const double delz = buf.drz[jj]; //x[j][2] - ztmp;
        const double rsq = delx*delx + dely*dely + delz*delz;
	const int jtype = buf.ext.type[jj]; //type[j];
        const int jelem = jtype; //map[jtype];

	double cut_ij = (radelem[jtype]+radelem[itype])*rcutfac;
	double cutsq_ij=cut_ij*cut_ij;

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
      
      //snaptr->reorder_rij(ninside);
      //const double chk_rij = snaptr->chk_rij(ninside);
      //const bool chosen = fabs(chk_rij-1.862026)<0.00001;
      //snaptr->set_dbg_chosen(chosen);

      // compute Ui, Yi for atom I

      if (chemflag)
        snaptr->compute_ui(ninside, ielem);
      else
        snaptr->compute_ui(ninside, 0);

      // for neighbors of I within cutoff:
      // compute Fij = dEi/dRj = -dEi/dRi
      // add to Fi, subtract from Fj
      // scaling is that for type I
/*
      std::vector<double> betaloc;
      betaloc.clear();
      betaloc.resize(ncoeff);
*/
      //      const double * const coeffi = coeffelem /*[ielem]*/;
/*
      for(int icoeff=0;icoeff<ncoeff;icoeff++)
      {
        long coeffelem_idx = itype * (ncoeff + 1 ) + icoeff + 1;
        assert( (coeffelem_idx >= 0) && (coeffelem_idx < coeffelem_size) );
	      betaloc[ icoeff ] = coeffelem[ coeffelem_idx ];
	    }
*/
      auto * betaloc = coeffelem + itype * (ncoeff + 1 ) + 1;
      //      std::cout << " My Type = " << itype << std::endl;

      // for(int icoeff=0;icoeff<ncoeff;icoeff++)
      // 	{
      // 	  std::cout << " Coeff " << icoeff << " = " << betaloc[icoeff] << std::endl;
      // 	}

      //      snaptr->compute_yi( /* beta[ii] ==> */ beta + ( ncoeff * ( cell_particle_offset[buf.cell] + buf.part ) ) );
      snaptr->compute_yi( betaloc /*.data()*/ );

      for (int jj = 0; jj < ninside; jj++)
      {
        int j = snaptr->inside[jj];
        snaptr->compute_duidrj(jj);

        double fij[3];
	fij[0]=0.;
	fij[1]=0.;
	fij[2]=0.;	
        snaptr->compute_deidrj(fij);
        
	//        if( conv_energy_units )
	//        {
	fij[0] *= conv_energy_factor;
	fij[1] *= conv_energy_factor;
	fij[2] *= conv_energy_factor;
	//        }

	//        [[maybe_unused]] const Mat3dT v_contrib = tensor( Vec3d{ fij[0] , fij[1] , fij[2] }, Vec3d{ snaptr->rij[jj][0] , snaptr->rij[jj][1] , snaptr->rij[jj][2] } );
	Mat3d v_contrib = tensor( Vec3d{ fij[0] , fij[1] , fij[2] }, Vec3d{ snaptr->rij[jj][0] , snaptr->rij[jj][1] , snaptr->rij[jj][2] } );
        if constexpr ( compute_virial ) { _vir += v_contrib * -1.0; }

        /*
        if(chosen)
        {
          printf("SNAPDBG: rij[%d]: x=% .3e y=% .3e z=% .3e , force: x=% .3e y=% .3e z=% .3e\n",jj,snaptr->rij[jj][0] , snaptr->rij[jj][1] , snaptr->rij[jj][2], fij[0] , fij[1] , fij[2] );
        }
        */

        _fx /*f[i][0]*/ += fij[0];//*scale[itype][itype];
        _fy /*f[i][1]*/ += fij[1];//*scale[itype][itype];
        _fz /*f[i][2]*/ += fij[2];//*scale[itype][itype];
        
        // f[j][0] -= fij[0]*scale[itype][itype];
        // f[j][1] -= fij[1]*scale[itype][itype];
        // f[j][2] -= fij[2]*scale[itype][itype];

        size_t cell_b=0, p_b=0;
        buf.nbh.get(j, cell_b, p_b);
        auto& lock_b = locks[cell_b][p_b];
        lock_b.lock();
        cells[cell_b][field::fx][p_b] -= fij[0];//*scale[itype][itype];
        cells[cell_b][field::fy][p_b] -= fij[1];//*scale[itype][itype];
        cells[cell_b][field::fz][p_b] -= fij[2];//*scale[itype][itype];
	// PL : for virial, no reciprocal contribution apparently.
	//        if constexpr ( compute_virial ) { cells[cell_b][field::virial][p_b] += v_contrib * 0.5; }
        lock_b.unlock();
      }

      if (eflag) {
	//        assert(ielem==0);
        const long bispectrum_ii_offset = ncoeff * ( cell_particle_offset[buf.cell] + buf.part );

        // evdwl = energy of atom I, sum over coeffs_k * Bi_k
        const double * const coeffi = coeffelem /*[ielem]*/;
	//	coeffelem[ itype * (ncoeff + 1 ) + icoeff + 1 ]
	double evdwl = coeffi[itype * (ncoeff + 1)];
	//	std::cout << "TYpe = " << itype << "e0 = " << evdwl << std::endl;

        // snaptr->copy_bi2bvec();

        // E = beta.B + 0.5*B^t.alpha.B

        // linear contributions

        for (int icoeff = 0; icoeff < ncoeff; icoeff++)
          evdwl += betaloc[icoeff] * bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;
	  //          evdwl += coeffi[itype * (ncoeff + 1) + icoeff+1] * bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;

        // quadratic contributions

        if (quadraticflag) {
          int k = ncoeff+1;
          for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
            double bveci = bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;
            evdwl += 0.5*coeffi[k++]*bveci*bveci;
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
              double bvecj = bispectrum[ bispectrum_ii_offset + jcoeff ] /*bispectrum[ii][jcoeff]*/ ;
              evdwl += coeffi[k++]*bveci*bvecj;
            }
          }
        }
        evdwl *= 1.;//scale[itype][itype];
        _en = evdwl; // ev_tally_full(i,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);
        
        if( conv_energy_units )
        {
          _en *= conv_energy_factor;
        }
      }

      // end of central atom calculation, atomically add contribution

      lock_a.lock();
      en += _en;
      fx += _fx;
      fy += _fy;
      fz += _fz;
      if constexpr ( compute_virial ) { virial += _vir; }
      lock_a.unlock();
      
    }
  };

}


