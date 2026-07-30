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

#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/math/basic_types.h>
#include <onika/memory/allocator.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/compute/compute_cell_particle_pairs_chunk.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <ptm_constants.h>

#include <vector>

// Runs right after compute_ptm: recovers atoms whose own PTM match was rejected (ptm_type ==
// PTM_MATCH_NONE) by majority vote among their physical neighbors' CURRENT ptm_type -- if a clear
// majority already match some structure, this atom is very likely part of the same
// undistorted lattice, just individually noisier/more strained than PTM's own rmsd_cutoff
// tolerates, not a genuine defect. Repeated for n_iterations, letting newly-recovered atoms help
// recover their own neighbors in turn (matching a real, recovered atom's own neighbors are
// themselves probably fine too -- the region actually SHRINKS each pass, same idea, not the same
// mechanism, as the "majority-vote relaxation" step in reference DXA implementations checked
// against, see compute_dxa_tet_classification's own comment for that investigation).
//
// This exists because compute_dxa_tet_classification's edge-vector-closure test used to try to
// do this job implicitly and overclassified a dislocation's ordinary elastic strain field as
// defective; per-atom neighbor consensus, run once up front, is the simpler, better-targeted fix
// both reference DXA implementations actually use (see that operator's comment).
//
// Uses ComputePairBuffer2<false,true> (UseNeighbors=true) to get each neighbor's own (cell,part)
// identity directly from the buffer (DefaultComputePairBufferAppendFunc already stores this
// unconditionally, just a no-op when UseNeighbors=false -- see compute_pair_buffer.h) rather than
// compute_ptm's flat-buffer-index-lookup pattern, since here (unlike compute_ptm's own env
// construction) neighbor IDENTITY, not just relative position, is what's actually needed.
namespace exaStamp
{
  using namespace exanb;

  struct PTMShrinkDisorderFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const double * const __restrict__ m_struct_type_in = nullptr;  // this iteration's starting point (read-only)
    double * const __restrict__ m_struct_type_out = nullptr;       // next iteration's buffer, pre-seeded with a copy of m_struct_type_in
    double m_min_fraction = 0.5;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () ( int jnum, ComputeBufferT& buf, CellParticlesT /*cells*/ ) const
    {
      const size_t i = m_cell_particle_offset[buf.cell] + buf.part;
      if( m_struct_type_in[i] != static_cast<double>(PTM_MATCH_NONE) || jnum <= 0 ) { return; } // already matched, or isolated -- untouched

      // small fixed-size histogram: ENABLED_FLAGS in compute_ptm.cu only ever produces
      // PTM_MATCH_NONE(0)/FCC(1)/HCP(2)/BCC(3)/ICO(4)/SC(5)
      int counts[6] = {0,0,0,0,0,0};
      for(int j=0;j<jnum;j++)
      {
        size_t c=0, p=0;
        buf.nbh.get(j,c,p);
        const int t = static_cast<int>( m_struct_type_in[ m_cell_particle_offset[c] + p ] );
        if( t>=0 && t<6 ) { ++counts[t]; }
      }

      int best_t = -1, best_count = 0;
      for(int t=1;t<6;t++) { if( counts[t] > best_count ) { best_count = counts[t]; best_t = t; } } // skip t=0 (NONE)

      if( best_t >= 1 && static_cast<double>(best_count)/static_cast<double>(jnum) >= m_min_fraction )
      {
        m_struct_type_out[i] = static_cast<double>(best_t);
      }
    }
  };

  template<class GridT>
  class PTMShrinkDisorder : public OperatorNode
  {
    ADD_SLOT( GridT               , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain              , domain          , INPUT , REQUIRED );
    ADD_SLOT( double              , rcut            , INPUT , REQUIRED , DocString{"Neighbor search cutoff -- same convention as compute_ptm's own rcut, typically the exact same value"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double              , rcut_max        , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , ptm_struct_type , INPUT_OUTPUT , DocString{"compute_ptm's flat structure-type buffer, refined in place: an unmatched (PTM_MATCH_NONE) particle whose neighbor majority already matches some structure gets that structure's type"} );
    ADD_SLOT( double              , min_neighbor_fraction , INPUT , 0.5 , DocString{"Minimum fraction of an unmatched particle's neighbors that must agree on the same structure type for it to be recovered"} );
    ADD_SLOT( long                , n_iterations    , INPUT , 3 , DocString{"Number of majority-vote passes -- each pass only sees the PREVIOUS pass's result, so a recovered particle can help recover its own neighbors in a later pass, letting the recovered region grow outward from confidently-matched bulk toward a genuine defect core, not just one shell deep"} );
    ADD_SLOT( long                , n_recovered     , OUTPUT , DocString{"Number of particles recovered (were PTM_MATCH_NONE, now assigned a structure type) across all iterations"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      *rcut_max = std::max( *rcut , *rcut_max );
      if( grid->number_of_cells() == 0 ) return;

      if( chunk_neighbors->m_chunk_size != 1 )
      {
        fatal_error() << "ptm_shrink_disorder: requires chunk_neighbors built with chunk_size=1 (set random_access: true"
                       << " on the chunk_neighbors config), got chunk_size=" << chunk_neighbors->m_chunk_size << std::endl;
      }

      const size_t total_particles = grid->number_of_particles();
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();

      std::vector<double> buf_a( ptm_struct_type->begin(), ptm_struct_type->end() );
      std::vector<double> buf_b( buf_a );

      using ComputeBuffer = ComputePairBuffer2<false,true>; // UseNeighbors=true: buf.nbh.get(j,cell,part) per neighbor
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      ComputePairNullWeightIterator cp_weight{};
      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks );
      static constexpr onika::FlatTuple<> compute_field_set = {};
      static constexpr DefaultPositionFields posfields = {};
      static constexpr ComputeParticlePairOpts<false,true,false> cp_opts = {}; // Symmetric=false, PreferComputeBuffer=true

      const IJK dims = grid->dimension();
      const ssize_t gl = grid->ghost_layers();
      const auto cells = grid->cells();
      const double rcut2 = (*rcut) * (*rcut);
      const double min_fraction = *min_neighbor_fraction;
      const long n_iter = std::max( 1L, *n_iterations );

      for(long iter=0; iter<n_iter; iter++)
      {
        buf_b = buf_a; // unmatched stays unmatched this pass unless recovered below
        PTMShrinkDisorderFunctor compute_op = { cell_particle_offset, buf_a.data(), buf_b.data(), min_fraction };

#       pragma omp parallel for collapse(3) schedule(dynamic)
        for(ssize_t k=gl;k<dims.k-gl;k++)
        for(ssize_t j=gl;j<dims.j-gl;j++)
        for(ssize_t i=gl;i<dims.i-gl;i++)
        {
          const IJK cell_a_loc{i,j,k};
          const size_t cell_a = static_cast<size_t>( grid_ijk_to_index(dims,cell_a_loc) );
          compute_cell_particle_pairs_cell( cells, dims, cell_a_loc, cell_a, rcut2
                                           , compute_buf, optional, compute_op
                                           , compute_field_set, onika::UIntConst<1>{}, cp_opts
                                           , posfields, std::index_sequence<>{} );
        }

        buf_a.swap( buf_b );
      }

      size_t n_rec = 0;
      for(size_t p=0;p<total_particles;p++)
      {
        if( (*ptm_struct_type)[p] == static_cast<double>(PTM_MATCH_NONE) && buf_a[p] != static_cast<double>(PTM_MATCH_NONE) ) { ++n_rec; }
        (*ptm_struct_type)[p] = buf_a[p];
      }

      *n_recovered = static_cast<long>( n_rec );
      lout << "ptm_shrink_disorder: " << n_rec << " / " << total_particles << " particles recovered by neighbor consensus" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Recovers atoms whose own PTM match was rejected (compute_ptm's ptm_type == PTM_MATCH_NONE) by
majority vote among their physical neighbors' current structure type: if a clear majority
(min_neighbor_fraction) already match some structure, this atom is very likely part of the same
undistorted lattice -- just individually noisier/more strained than compute_ptm's own rmsd_cutoff
tolerates, not a genuine defect. Run right after compute_ptm, before ptm_fields/compute_delaunay.
Repeated for n_iterations so a recovered particle can help recover its own neighbors in a later
pass. A genuine dislocation core's own neighbors are themselves disordered too, so there's no
majority to recover it with -- this only shrinks spurious/noisy misclassification, not real
defects.

Usage example:

compute_ptm: { rcut: 5.0 ang, rmsd_cutoff: 0.2 }
ptm_shrink_disorder: { rcut: 5.0 ang, min_neighbor_fraction: 0.5, n_iterations: 3 }
ptm_fields: {}

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(ptm_shrink_disorder)
  {
    OperatorNodeFactory::instance()->register_factory( "ptm_shrink_disorder", make_grid_variant_operator< PTMShrinkDisorder > );
  }

}
