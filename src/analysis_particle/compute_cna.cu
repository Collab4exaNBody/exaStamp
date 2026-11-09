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
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exaStamp/analysis_particle/bounded_pair_buffer.h>

#include <omp.h>
#include <vector>
#include <cmath>
#include <algorithm>

// Adaptive-cutoff Common Neighbor Analysis: FCC, HCP, ICO (12 neighbors) and BCC (14 neighbors).
// Output uses the PTM_MATCH_* numbering of compute_ptm (0 none, 1 FCC, 2 HCP, 3 BCC, 4 ICO).
namespace exaStamp
{
  using namespace exanb;

  // same values as PTM_MATCH_*
  static constexpr double CNA_MATCH_NONE = 0.0;
  static constexpr double CNA_MATCH_FCC  = 1.0;
  static constexpr double CNA_MATCH_HCP  = 2.0;
  static constexpr double CNA_MATCH_BCC  = 3.0;
  static constexpr double CNA_MATCH_ICO  = 4.0;

  static constexpr int CNA_MAX_NBRS = 16; // BCC 14 + 1

  struct CNAFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    double * const __restrict__ m_struct_type_out = nullptr; // 1 per particle

    // nn nearest neighbors (sorted), adaptive cutoff = mean distance of the scaling_nbrs nearest
    // * scaling_rescale * (1+sqrt(2))/2
    template<class ComputeBufferT>
    static inline double classify_family( const ComputeBufferT& buf, int jnum, int nn, int scaling_nbrs, double scaling_rescale )
    {
      if( jnum < nn+1 ) { return CNA_MATCH_NONE; }

      double scaling = 0.0;
      for(int i=0;i<scaling_nbrs;i++) { scaling += std::sqrt( buf.d2[i] ); }
      scaling /= static_cast<double>(scaling_nbrs);
      const double cutoff = scaling * scaling_rescale * (1.0+std::sqrt(2.0)) * 0.5;
      const double cutoff2 = cutoff*cutoff;

      if( buf.d2[nn] <= cutoff2 ) { return CNA_MATCH_NONE; } // more than nn neighbors within cutoff

      bool bonded[CNA_MAX_NBRS][CNA_MAX_NBRS];
      for(int a=0;a<nn;a++)
      {
        bonded[a][a] = false;
        for(int b=a+1;b<nn;b++)
        {
          const double dx = buf.drx[a]-buf.drx[b], dy = buf.dry[a]-buf.dry[b], dz = buf.drz[a]-buf.drz[b];
          const bool is_bonded = (dx*dx+dy*dy+dz*dz) < cutoff2;
          bonded[a][b] = bonded[b][a] = is_bonded;
        }
      }

      int n421=0, n422=0, n555=0, n444=0, n666=0;
      for(int i=0;i<nn;i++)
      {
        int common[CNA_MAX_NBRS]; int nc=0;
        for(int j=0;j<nn;j++) { if( j!=i && bonded[i][j] ) { common[nc++] = j; } }

        // maxChainLength = largest number of bonds in a connected component of the common neighbors
        int numNeighborBonds = 0;
        int parent[CNA_MAX_NBRS];
        for(int k=0;k<nc;k++) { parent[k] = k; }
        struct { int* p; int operator()(int x){ while(p[x]!=x){p[x]=p[p[x]];x=p[x];} return x; } } find{parent};
        for(int a=0;a<nc;a++)
        {
          for(int b=a+1;b<nc;b++)
          {
            if( bonded[ common[a] ][ common[b] ] )
            {
              ++numNeighborBonds;
              const int ra = find(a), rb = find(b);
              if( ra != rb ) { parent[ra] = rb; }
            }
          }
        }
        int edgecount[CNA_MAX_NBRS] = {0};
        for(int a=0;a<nc;a++)
        {
          for(int b=a+1;b<nc;b++)
          {
            if( bonded[ common[a] ][ common[b] ] ) { ++edgecount[ find(a) ]; }
          }
        }
        int maxChainLength = 0;
        for(int k=0;k<nc;k++) { maxChainLength = std::max( maxChainLength, edgecount[k] ); }

        if     ( nc==4 && numNeighborBonds==2 && maxChainLength==1 ) { ++n421; }
        else if( nc==4 && numNeighborBonds==2 && maxChainLength==2 ) { ++n422; }
        else if( nc==5 && numNeighborBonds==5 && maxChainLength==5 ) { ++n555; }
        else if( nc==4 && numNeighborBonds==4 && maxChainLength==4 ) { ++n444; }
        else if( nc==6 && numNeighborBonds==6 && maxChainLength==6 ) { ++n666; }
      }

      if( nn==12 )
      {
        if( n421==12 ) { return CNA_MATCH_FCC; }
        if( n421==6 && n422==6 ) { return CNA_MATCH_HCP; }
        if( n555==12 ) { return CNA_MATCH_ICO; }
      }
      else // nn==14
      {
        if( n666==8 && n444==6 ) { return CNA_MATCH_BCC; }
      }
      return CNA_MATCH_NONE;
    }

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () ( int jnum, ComputeBufferT& buf, CellParticlesT /*cells*/ ) const
    {
      // buf holds the nearest neighbors, sorted by distance (NearestNeighborsComputeBuffer)
      const int n = std::min( jnum, CNA_MAX_NBRS );

      double struct_type = classify_family( buf, n, 12, 12, 1.0 );
      if( struct_type == CNA_MATCH_NONE )
      {
        struct_type = classify_family( buf, n, 14, 8, 2.0/std::sqrt(3.0) ); // BCC: first shell rescaled to the second
      }

      const size_t i = m_cell_particle_offset[buf.cell] + buf.part;
      m_struct_type_out[i] = struct_type;
    }
  };

  template<class GridT>
  class ComputeCNA : public OperatorNode
  {
    ADD_SLOT( GridT               , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain              , domain          , INPUT , REQUIRED );
    ADD_SLOT( double              , rcut            , INPUT , REQUIRED , DocString{"Neighbor search cutoff, must include 13 neighbors (FCC/HCP/ICO) or 15 (BCC)"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double              , rcut_max        , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , cna_struct_type , OUTPUT , DocString{"Per-particle structure type (0 none, 1 FCC, 2 HCP, 3 BCC, 4 ICO), indexed by cell_particle_offset[cell] + particle"} );
    ADD_SLOT( long , cna_n_matched , OUTPUT , DocString{"Number of matched particles on this rank"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      *rcut_max = std::max( *rcut , *rcut_max );
      if( grid->number_of_cells() == 0 ) return;

      if( chunk_neighbors->m_chunk_size != 1 )
      {
        fatal_error() << "compute_cna: requires chunk_neighbors built with chunk_size=1 (set random_access: true"
                       << " on the chunk_neighbors config), got chunk_size=" << chunk_neighbors->m_chunk_size << std::endl;
      }

      const size_t total_particles = grid->number_of_particles();
      cna_struct_type->resize( total_particles );

      using ComputeBuffer = NearestNeighborsComputeBuffer< CNA_MAX_NBRS >;
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      CNAFunctor compute_op = { grid->cell_particle_offset_data(), cna_struct_type->data() };

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

#     pragma omp parallel for collapse(3) schedule(dynamic)
      for(ssize_t k=gl;k<dims.k-gl;k++)
      for(ssize_t j=gl;j<dims.j-gl;j++)
      for(ssize_t i=gl;i<dims.i-gl;i++)
      {
        const IJK cell_a_loc{i,j,k};
        const size_t cell_a = static_cast<size_t>( grid_ijk_to_index(dims,cell_a_loc) );
        compute_cell_particle_pairs_cell( cells, dims, cell_a_loc, cell_a, rcut2
                                         , compute_buf, optional, compute_op
                                         , compute_field_set, cp_opts
                                         , posfields, std::index_sequence<>{} );
      }

      size_t n_matched = 0;
      for(size_t p=0;p<total_particles;p++) { if( (*cna_struct_type)[p] != CNA_MATCH_NONE ) { ++n_matched; } }

      *cna_n_matched = static_cast<long>(n_matched);
      lout << "compute_cna: " << n_matched << " / " << total_particles << " particles matched" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Adaptive-cutoff Common Neighbor Analysis: FCC, HCP, BCC and ICO. Writes a flat per-particle
buffer, cna_struct_type (0 none, 1 FCC, 2 HCP, 3 BCC, 4 ICO, same numbering as compute_ptm).
Use cna_fields to copy it into a grid field.

Usage example:

compute_cna: { rcut: 5.0 ang }

)EOF";
    }
  };

  // cna_struct_type -> grid field
  struct CNAFieldsFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    const double * const __restrict__ m_struct_type_in = nullptr;

    inline void operator () ( size_t cell, unsigned int part, double& struct_type_out ) const
    {
      struct_type_out = m_struct_type_in[ m_cell_particle_offset[cell] + part ];
    }
  };

  template<class GridT>
  class CNAFields : public OperatorNode
  {
    ADD_SLOT( GridT , grid , INPUT_OUTPUT );
    ADD_SLOT( onika::memory::CudaMMVector<double> , cna_struct_type , INPUT , REQUIRED );
    ADD_SLOT( std::string , struct_field , INPUT , std::string("cna_type") , DocString{"Name of the resulting per-particle structure-type grid field"} );

  public:
    inline void execute () override final
    {
      auto struct_acc = grid->field_accessor( field::mk_generic_real( *struct_field ) );
      CNAFieldsFunctor func = { grid->cell_particle_offset_data(), cna_struct_type->data() };
      compute_cell_particles( *grid, false, func, onika::make_flat_tuple( struct_acc ), parallel_execution_context() );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Copies compute_cna's cna_struct_type into a per-particle grid field. Owned particles only: run
ghost_update_opt on the field if ghosts need it.

Usage example:

compute_cna: { rcut: 5.0 ang }
cna_fields: { struct_field: cna_type }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(compute_cna)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_cna", make_grid_variant_operator< ComputeCNA > );
    OperatorNodeFactory::instance()->register_factory( "cna_fields", make_grid_variant_operator< CNAFields > );
  }

}
