#pragma once
#include <onika/memory/allocator.h>
#include <onika/cuda/cuda.h>
#include <vector>
#include <utility>
#include <exanb/compute/compute_pair_buffer.h>

#include <exanb/fields.h>
#include <exanb/core/config.h>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  struct PhiRhoCutoff
  {
    double m_phi_cutoff;
    double m_rho_cutoff;
  };

  struct EamPotentialScratch
  {
    //onika::memory::CudaMMVector< size_t > m_offset; 
    onika::memory::CudaMMVector< double > m_emb;
    onika::memory::CudaMMVector< PhiRhoCutoff > m_phi_rho_cutoff;
    onika::memory::CudaMMVector< uint8_t > m_pair_enabled;
  };
  
/*
  // temporary storage for emb terms, used to pass emb across the 2 compute passes
  struct EamScratchBuffer
  {
    std::vector< std::vector<double> > m_emb;
    std::vector< std::pair<double,double> > m_phi_rho_cutoff;
    std::vector< bool > m_pair_enabled;
  };
*/

  // additional storage space added to compute buffer created by compute_cell_particle_pairs
  struct alignas(DEFAULT_ALIGNMENT) EamComputeBufferExt
  {
    alignas(DEFAULT_ALIGNMENT) double emb[exanb::MAX_PARTICLE_NEIGHBORS];
  };

  // additional storage space added to compute buffer created by compute_cell_particle_pairs
  struct alignas(DEFAULT_ALIGNMENT) EamComputeBufferEmbTypeExt
  {
    alignas(DEFAULT_ALIGNMENT) double emb[exanb::MAX_PARTICLE_NEIGHBORS];
    alignas(DEFAULT_ALIGNMENT) uint16_t type[exanb::MAX_PARTICLE_NEIGHBORS];
  };

  // additional storage space added to compute buffer created by compute_cell_particle_pairs
  struct alignas(DEFAULT_ALIGNMENT) EamComputeBufferTypeExt
  {
    alignas(DEFAULT_ALIGNMENT) uint16_t type[exanb::MAX_PARTICLE_NEIGHBORS];
  };
   
  // functor that populate compute buffer's extended storage for particle charges
  struct EamCopyParticleEmb
  {
    const size_t* m_cell_emb_offset = nullptr;
    const double* m_particle_emb = nullptr;

    template<typename ComputeBufferT, typename FieldArraysT, class NbhDataT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, const FieldArraysT * cells, size_t cell_b, size_t p_b, const NbhDataT& nbh_data) const noexcept
    {
      assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
      tab.ext.emb[tab.count] = m_particle_emb[ m_cell_emb_offset[cell_b] + p_b ];
      DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, nbh_data );
    }
  };

  // functor that populate compute buffer's extended storage for particle emb and atom type
  struct EamCopyParticleEmbType
  {
    const size_t* m_cell_emb_offset = nullptr;
    const double* m_particle_emb = nullptr;
    
    template<typename ComputeBufferT, typename FieldArraysT, class NbhDataT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, const FieldArraysT * cells, size_t cell_b, size_t p_b, const NbhDataT& nbh_data) const noexcept
    {
      assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
      tab.ext.emb[tab.count] = m_particle_emb[ m_cell_emb_offset[cell_b] + p_b ];
      tab.ext.type[tab.count] = cells[cell_b][field::type][p_b];
      DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, nbh_data );
    }
  };

  struct EamCopyParticleEmbInitFunc
  {
    const size_t* m_cell_emb_offset = nullptr;
    double* m_particle_emb = nullptr;
    template<bool W,bool N, class BufExt, class CopyEmbOp, size_t MaxN>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ComputePairBuffer2<W,N,BufExt,CopyEmbOp,MaxN>& cpbuf ) const
    {
      cpbuf.process_neighbor = { m_cell_emb_offset , m_particle_emb };
    }
  };

  // functor that populate compute buffer's extended storage for particle emb and atom type
  struct EamCopyParticleType
  {
    template<typename ComputeBufferT, typename FieldArraysT, class NbhDataT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, const FieldArraysT * cells, size_t cell_b, size_t p_b, const NbhDataT& nbh_data) const noexcept
    {
      assert( ssize_t(tab.count) < ssize_t(tab.MaxNeighbors) );
      tab.ext.type[tab.count] = cells[cell_b][field::type][p_b];
      DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, nbh_data );
    }
  };

  struct EAMSpecyPairInfo
  {
    double m_charge_a = 0.0;
    double m_charge_b = 0.0;
  };

  template<class EamParametersT>
  struct EamMultimatParameters
  {
    EamParametersT m_parameters;
    EAMSpecyPairInfo m_specy_pair;
    std::string m_type_a;
    std::string m_type_b;
  };

}

