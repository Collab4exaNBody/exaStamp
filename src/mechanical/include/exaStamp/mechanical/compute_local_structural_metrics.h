#pragma once

#include <exanb/core/basic_types.h>
#include <exanb/compute/compute_pair_buffer.h>

#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>

#include <exaStamp/potential/snaplegacy/SnapLegacyCG.h>
#include <exaStamp/potential/snaplegacy/SnapLegacyBS.h>

#include <vector>
#include <algorithm>

#include <omp.h>

namespace exaStamp
{
  using namespace exanb;

    static inline bool compareNorm(Vec3d v1, Vec3d v2){return (norm(v1) < norm(v2));}

    using onika::memory::DEFAULT_ALIGNMENT;
    //using ComputeBuffer = ComputePairBuffer2<false,false>;
 
    struct PerThreadContext
    {
      std::shared_ptr<SnapLegacyBS> m_snapbs = nullptr;
    };

    struct alignas(DEFAULT_ALIGNMENT) BispectrumOp 
    {
      const SnapLegacyCG& m_cg;
      std::vector<PerThreadContext>& m_thread_ctx;
      const double m_rcut_pot;
      const double m_rcut_bs;
      const int m_nneigh_bs;
      const bool m_closest_bs;
      GridParticleLocalStructuralMetrics & m_local_structural_data;
      
      template<class CellParticlesT>
    	inline void operator ()
        (
    	 size_t n,
    	 ComputePairBuffer2<false,false>& buf,
    	 CellParticlesT*
    	 ) const	
      {
    	GridParticleLocalStructuralMetrics & local_structural_data = m_local_structural_data;

    	// get thread specific compute context
        size_t tid = omp_get_thread_num();
        assert( tid < m_thread_ctx.size() );
        SnapLegacyBS& snap_bs = * m_thread_ctx[tid].m_snapbs;

	vector<Vec3d> dist_neigh;
        dist_neigh.resize(n);
        for(size_t i=0;i<n;i++)
          {
            dist_neigh[i].x = buf.drx[i];
            dist_neigh[i].y = buf.dry[i];
            dist_neigh[i].z = buf.drz[i];
          }
        sort(dist_neigh.begin(), dist_neigh.end(), compareNorm);
        for(size_t i=0;i<n;i++)
          {
            buf.drx[i] = dist_neigh[i].x;
            buf.dry[i] = dist_neigh[i].y;
            buf.drz[i] = dist_neigh[i].z;
          }

	if (m_closest_bs)
	  {
	    double rcut = norm(dist_neigh[m_nneigh_bs-1])+0.01;
	    snap_bs.set_neighbours( buf.drx, buf.dry, buf.drz, nullptr, rcut , n);
	    snap_bs.compute_cmm(rcut);
	    snap_bs.compute_bs( 0, rcut, m_cg );
	  }
	else
	  {
	    double rcut = 0.;
	    if (m_rcut_bs > 0.1)
	      {
		rcut = std::min(m_rcut_bs,m_rcut_pot);
	      }
	    else if (m_nneigh_bs > 0)
	      {
		double nneigh_pot = n*1.;
		rcut = std::min(m_rcut_pot * pow(m_nneigh_bs/nneigh_pot,1./3.), m_rcut_pot);
	      }
	    
	    snap_bs.set_neighbours( buf.drx, buf.dry, buf.drz, nullptr, rcut , n);
	    snap_bs.compute_cmm(rcut);
	    snap_bs.compute_bs( 0, rcut, m_cg );
	  }

    	int nbs = SnapLegacyBS::n_idx_bs( m_cg.get_jmax()*2 );
        local_structural_data[ buf.cell ].bispectrum[ buf.part ].resize( nbs );
        for(int i=0;i<nbs;i++)
    	  {
    	    complex<double> bs = snap_bs.bs_val(i);
    	    local_structural_data[ buf.cell ].bispectrum[ buf.part ] [i]= bs.real();  
    	  }
	
      }
      
    };

  
    // std::vector<PerThreadContext> m_thread_ctx;
    // std::shared_ptr<SnapLegacyCG> m_cg = nullptr;
    // double m_rcut = 0.0;
    // int m_cg_nt = 2;    
}
