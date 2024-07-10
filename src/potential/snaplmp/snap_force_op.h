#pragma once

#include <exanb/core/physics_constants.h>
#include <cmath>

namespace exaStamp
{

  using namespace exanb;

  static inline void snap_zero_uarraytot( // READ ONLY
                                          int nelements
                                        , int twojmax
                                        , int const * idxu_block
                                        , int idxu_max
                                        , double wself
                                        , bool wselfall_flag
                                          // WRITE ONLY
                                        , double * ulisttot_r, double * ulisttot_i
                                          // ORIGINAL PARAMETERS
                                        , int ielem ) 
  {
    for (int jelem = 0; jelem < nelements; jelem++)
    for (int j = 0; j <= twojmax; j++) {
      int jju = idxu_block[j];
      for (int mb = 0; mb <= j; mb++) {
        for (int ma = 0; ma <= j; ma++) {
          ulisttot_r[jelem*idxu_max+jju] = 0.0;
          ulisttot_i[jelem*idxu_max+jju] = 0.0;

          // utot(j,ma,ma) = wself, sometimes
          if (jelem == ielem || wselfall_flag)
            if (ma==mb)
            ulisttot_r[jelem*idxu_max+jju] = wself; ///// double check this
          jju++;
        }
      }
    }
  }

  static inline void snap_compute_uarray( // READ ONLY
                                          int twojmax
                                        , int const * idxu_block
                                        , double const * const * rootpqarray
                                          // WRITE ONLY
                                        , double * const * ulist_r_ij, double * const * ulist_i_ij
                                          // ORIGINAL PARAMETERS
                                        , double x, double y, double z, double z0, double r, int jj )
  {
    double r0inv;
    double a_r, b_r, a_i, b_i;
    double rootpq;

    // compute Cayley-Klein parameters for unit quaternion

    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;

    // VMK Section 4.8.2
    double* ulist_r = ulist_r_ij[jj];
    double* ulist_i = ulist_i_ij[jj];

    ulist_r[0] = 1.0;
    ulist_i[0] = 0.0;

    for (int j = 1; j <= twojmax; j++) {
      int jju = idxu_block[j];
      int jjup = idxu_block[j-1];

      // fill in left side of matrix layer from previous layer

      for (int mb = 0; 2*mb <= j; mb++) {
        ulist_r[jju] = 0.0;
        ulist_i[jju] = 0.0;

        for (int ma = 0; ma < j; ma++) {
          rootpq = rootpqarray[j - ma][j - mb];
          ulist_r[jju] +=
            rootpq *
            (a_r * ulist_r[jjup] +
             a_i * ulist_i[jjup]);
          ulist_i[jju] +=
            rootpq *
            (a_r * ulist_i[jjup] -
             a_i * ulist_r[jjup]);

          rootpq = rootpqarray[ma + 1][j - mb];
          ulist_r[jju+1] =
            -rootpq *
            (b_r * ulist_r[jjup] +
             b_i * ulist_i[jjup]);
          ulist_i[jju+1] =
            -rootpq *
            (b_r * ulist_i[jjup] -
             b_i * ulist_r[jjup]);
          jju++;
          jjup++;
        }
        jju++;
      }

      // copy left side to right side with inversion symmetry VMK 4.4(2)
      // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

      jju = idxu_block[j];
      jjup = jju+(j+1)*(j+1)-1;
      int mbpar = 1;
      for (int mb = 0; 2*mb <= j; mb++) {
        int mapar = mbpar;
        for (int ma = 0; ma <= j; ma++) {
          if (mapar == 1) {
            ulist_r[jjup] = ulist_r[jju];
            ulist_i[jjup] = -ulist_i[jju];
          } else {
            ulist_r[jjup] = -ulist_r[jju];
            ulist_i[jjup] = ulist_i[jju];
          }
          mapar = -mapar;
          jju++;
          jjup--;
        }
        mbpar = -mbpar;
      }
    }
  }

  static inline double snap_compute_sfac( // READ ONLY
                                          double rmin0, bool switch_flag, bool switch_inner_flag
                                          // ORIGINAL PARAMTERS
                                        , double r, double rcut, double sinner, double dinner)
  {
    double sfac = 0.0;

    // calculate sfac = sfac_outer

    if (switch_flag == 0) sfac = 1.0;
    else if (r <= rmin0) sfac = 1.0;
    else if (r > rcut) sfac = 0.0;
    else {
      double rcutfac = M_PI / (rcut - rmin0);
      sfac = 0.5 * (cos((r - rmin0) * rcutfac) + 1.0);
    }

    // calculate sfac *= sfac_inner, rarely visited

    if (switch_inner_flag == 1 && r < sinner + dinner) {
      if (r > sinner - dinner) {
        double rcutfac = (M_PI/2) / dinner;
        sfac *= 0.5 * (1.0 - cos( (M_PI/2) + (r - sinner) * rcutfac));
      } else sfac = 0.0;
    }

    return sfac;
  }

  static inline void snap_add_uarraytot( // READ ONLY
                                         double rmin0, bool switch_flag, bool switch_inner_flag
                                       , int twojmax, bool chem_flag, int const * element
                                       , double const * rcutij, double const * sinnerij, double const * dinnerij, double const * wj
                                       , int const * idxu_block , int idxu_max
                                       , double const * const * ulist_r_ij, double const * const * ulist_i_ij
                                         // WRITE ONLY
                                       , double * ulisttot_r, double * ulisttot_i
                                         // ORIGINAL PARAMETERS
                                       , double r, int jj)
  {
    double sfac;
    int jelem;

    sfac = snap_compute_sfac( rmin0, switch_flag, switch_inner_flag, r, rcutij[jj], sinnerij[jj], dinnerij[jj]);
    sfac *= wj[jj];

    if (chem_flag) jelem = element[jj];
    else jelem = 0;

    double const * ulist_r = ulist_r_ij[jj];
    double const * ulist_i = ulist_i_ij[jj];

    for (int j = 0; j <= twojmax; j++) {
      int jju = idxu_block[j];
      for (int mb = 0; mb <= j; mb++)
        for (int ma = 0; ma <= j; ma++) {
          ulisttot_r[jelem*idxu_max+jju] +=
            sfac * ulist_r[jju];
          ulisttot_i[jelem*idxu_max+jju] +=
            sfac * ulist_i[jju];
          jju++;
        }
    }
  }


  static inline void snap_compute_ui( LAMMPS_NS::SNA * snaptr
                                    , double const * const * rij
                                    , double const * rcutij, double rmin0, double rfac0
                                    , /*parameters*/ int jnum, int ielem)
  {
    // utot(j,ma,mb) = 0 for all j,ma,ma
    // utot(j,ma,ma) = 1 for all j,ma
    // for j in neighbors of i:
    //   compute r0 = (x,y,z,z0)
    //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

    snap_zero_uarraytot( snaptr->nelements, snaptr->twojmax, snaptr->idxu_block, snaptr->idxu_max, snaptr->wself, snaptr->wselfall_flag
                         , snaptr->ulisttot_r, snaptr->ulisttot_i
                         , ielem );
    
    for (int j = 0; j < jnum; j++)
    {
      const double x = rij[j][0];
      const double y = rij[j][1];
      const double z = rij[j][2];
      const double rsq = x * x + y * y + z * z;
      const double r = sqrt(rsq);

      const double theta0 = (r - rmin0) * rfac0 * M_PI / (rcutij[j] - rmin0);
      const double z0 = r / tan(theta0);

      snap_compute_uarray( snaptr->twojmax, snaptr->idxu_block, snaptr->rootpqarray, snaptr->ulist_r_ij, snaptr->ulist_i_ij, x, y, z, z0, r, j);
      
      snap_add_uarraytot( snaptr->rmin0, snaptr->switch_flag, snaptr->switch_inner_flag, snaptr->twojmax, snaptr->chem_flag, snaptr->element
                        , snaptr->rcutij, snaptr->sinnerij, snaptr->dinnerij, snaptr->wj
                        , snaptr->idxu_block, snaptr->idxu_max, snaptr->ulist_r_ij, snaptr->ulist_i_ij, snaptr->ulisttot_r, snaptr->ulisttot_i
                        , r, j);
    }
  }


  // Force operator
  struct SnapLMPForceOp 
  {
    SnapLMPThreadContext* m_thread_ctx = nullptr;
    const size_t n_thread_ctx = 0;
    
    const size_t * const cell_particle_offset = nullptr;
    const double * const beta = nullptr;
    const double * const bispectrum = nullptr;
    
    const double * const coeffelem = nullptr;
    const unsigned int coeffelem_size = 0;
    const unsigned int ncoeff = 0;
    
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
        const int jtype = buf.nbh_pt[jj][field::type];; //buf.ext.type[jj]; 
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

      if (chemflag) snap_compute_ui(snaptr,snaptr->rij,snaptr->rcutij,snaptr->rmin0,snaptr->rfac0, ninside, ielem); // snaptr->compute_ui(ninside, ielem);
      else          snap_compute_ui(snaptr,snaptr->rij,snaptr->rcutij,snaptr->rmin0,snaptr->rfac0, ninside, 0);     // snaptr->compute_ui(ninside, 0);

      // for neighbors of I within cutoff:
      // compute Fij = dEi/dRj = -dEi/dRi
      // add to Fi, subtract from Fj
      // scaling is that for type I

      std::vector<double> betaloc;
      betaloc.clear();
      betaloc.resize(ncoeff);

      //      const double * const coeffi = coeffelem /*[ielem]*/;
      for(unsigned int icoeff=0;icoeff<ncoeff;icoeff++)
      {
        unsigned int coeffelem_idx = itype * (ncoeff + 1 ) + icoeff + 1;
        assert( (coeffelem_idx >= 0) && (coeffelem_idx < coeffelem_size) );
	      betaloc[ icoeff ] = coeffelem[ coeffelem_idx ];
	    }


      //      snaptr->compute_yi( /* beta[ii] ==> */ beta + ( ncoeff * ( cell_particle_offset[buf.cell] + buf.part ) ) );
      snaptr->compute_yi( betaloc.data() );

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

      if (eflag)
      {
	      // assert(ielem==0);
        const long bispectrum_ii_offset = ncoeff * ( cell_particle_offset[buf.cell] + buf.part );

        // evdwl = energy of atom I, sum over coeffs_k * Bi_k
        const double * const coeffi = coeffelem /*[ielem]*/;
	      //	coeffelem[ itype * (ncoeff + 1 ) + icoeff + 1 ]
	      double evdwl = coeffi[itype * (ncoeff + 1)];
	      //	std::cout << "TYpe = " << itype << "e0 = " << evdwl << std::endl;

        // snaptr->copy_bi2bvec();

        // E = beta.B + 0.5*B^t.alpha.B

        // linear contributions

        for (unsigned int icoeff = 0; icoeff < ncoeff; icoeff++)
        {
          evdwl += betaloc[icoeff] * bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;
	      // evdwl += coeffi[itype * (ncoeff + 1) + icoeff+1] * bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;
        }

        // quadratic contributions

        if (quadraticflag)
        {
          int k = ncoeff+1;
          for (unsigned int icoeff = 0; icoeff < ncoeff; icoeff++) 
          {
            double bveci = bispectrum[ bispectrum_ii_offset + icoeff ] /*bispectrum[ii][icoeff]*/ ;
            evdwl += 0.5*coeffi[k++]*bveci*bveci;
            for (unsigned int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) 
            {
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


