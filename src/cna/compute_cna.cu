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
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <omp.h>
#include <vector>
#include <cmath>
#include <algorithm>

// Adaptive-cutoff Common Neighbor Analysis (CNA), reimplemented from scratch following OVITO's
// own structure-identification engine (src/ovito/crystalanalysis/modifier/structureanalysis/
// StructureAnalysis.cpp) -- the exact classifier DXA itself uses under the hood, per OVITO's own
// user manual ("This is done with the help of the Common Neighbor Analysis (CNA) method"). Added
// specifically because compute_ptm's PTM-based classification (continuous RMSD fit) turned out
// to have a genuinely different (wider) strain tolerance than CNA's discrete bond-topology test,
// which was inflating the DXA pipeline's "bad" (defective) region well beyond OVITO's own DXA
// output on an identical dislocation quadrupole test case -- see src/delaunay/README.md.
//
// Checks FCC, HCP, BCC, ICO (not SC or diamond -- classic CNA doesn't cover SC at all; diamond
// needs its own, differently-shaped 16-vector template and isn't needed for this pipeline's
// FCC/HCP/BCC focus, can be added later if needed). Output uses the exact same PTM_MATCH_*
// numbering convention as compute_ptm's own ptm_struct_type (0=none,1=FCC,2=HCP,3=BCC,4=ICO), so
// it's a drop-in replacement for compute_dxa_edge_vectors' struct_field -- orientation still
// needs to come from ptm_fields' ptm_orientation (CNA's own reference implementation only
// produces a per-CLUSTER, not per-atom, orientation -- a materially bigger feature, not
// implemented here; an atom's own PTM-fitted orientation remains a reasonable per-atom proxy
// even when its own PTM structure-type match was rejected, same reasoning as ptm_shrink_disorder).
//
// Algorithm (adaptive CNA, per-neighbor-family since FCC/HCP/ICO need 12 neighbors and BCC needs
// 14, with different cutoff formulas):
//   1. Sort the particle's neighbors by distance (need at least nn+1 to check the "next neighbor
//      must lie outside the adaptive cutoff sphere" condition that pins down exactly nn).
//   2. Compute the local adaptive cutoff from the nearest bond lengths (FCC/HCP/ICO: mean of the
//      12 nearest; BCC: mean of the nearest 8, rescaled by 2/sqrt(3) to the 2nd-shell length
//      scale), scaled by (1+sqrt(2))/2 -- the exact StructureAnalysis.cpp formula.
//   3. Build the pairwise "common neighbor" bond graph among the nn candidate neighbors (bonded
//      iff their mutual distance is under the adaptive cutoff).
//   4. For each of the nn neighbors, using ITS row of the bond graph: numCommonNeighbors = how
//      many of the OTHER nn-1 neighbors it's bonded to; numNeighborBonds = how many bonded PAIRS
//      exist among that common-neighbor set; maxChainLength = the largest connected component's
//      size within that set (union-find). Tally the classic 4-2-1/4-2-2/5-5-5 (12-neighbor) or
//      4-4-4/6-6-6 (14-neighbor) signature counts across all nn neighbors.
//   5. Accept FCC (12x 4-2-1), HCP (6x 4-2-1 + 6x 4-2-2), ICO (12x 5-5-5), or BCC (6x 4-4-4 + 8x
//      6-6-6) -- exactly OVITO's own acceptance table, no tolerance.
namespace exaStamp
{
  using namespace exanb;

  // same numbering as PTM_MATCH_* (ptm_constants.h) -- kept independent here (no PTM link at
  // all) so any downstream consumer keying off this convention (e.g. compute_dxa_edge_vectors'
  // target_structure name -> integer mapping) works transparently regardless of which operator
  // (compute_ptm or compute_cna) produced the structure-type field.
  static constexpr double CNA_MATCH_NONE = 0.0;
  static constexpr double CNA_MATCH_FCC  = 1.0;
  static constexpr double CNA_MATCH_HCP  = 2.0;
  static constexpr double CNA_MATCH_BCC  = 3.0;
  static constexpr double CNA_MATCH_ICO  = 4.0;

  static constexpr int CNA_MAX_NBRS = 16; // enough for BCC's 14 + 1 sanity-check neighbor

  struct CNAFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    double * const __restrict__ m_struct_type_out = nullptr; // 1 per particle

    // classify the nn nearest (already sorted) neighbors as one member of a CNA family (12: FCC/
    // HCP/ICO; 14: BCC). scaling_nbrs/scaling_rescale determine the adaptive cutoff (see file
    // comment): cutoff = mean(d2[0..scaling_nbrs)) * scaling_rescale * (1+sqrt(2))/2.
    template<class ComputeBufferT>
    static inline double classify_family( const ComputeBufferT& buf, int jnum, int nn, int scaling_nbrs, double scaling_rescale )
    {
      if( jnum < nn+1 ) { return CNA_MATCH_NONE; } // need one extra neighbor for the sanity check below

      double scaling = 0.0;
      for(int i=0;i<scaling_nbrs;i++) { scaling += std::sqrt( buf.d2[i] ); }
      scaling /= static_cast<double>(scaling_nbrs);
      const double cutoff = scaling * scaling_rescale * (1.0+std::sqrt(2.0)) * 0.5;
      const double cutoff2 = cutoff*cutoff;

      if( buf.d2[nn] <= cutoff2 ) { return CNA_MATCH_NONE; } // (nn+1)-th neighbor too close: not exactly nn neighbors at this cutoff

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

        // maxChainLength is the largest EDGE count among connected components of the common-
        // neighbor bond graph -- not the node count. For closed rings (444/555/666) edges==nodes
        // so the two coincide, but for the open-chain FCC/HCP signatures (421: two disjoint edges,
        // 2 nodes/1 edge each; 422: one 2-edge/3-node path + 1 isolated node) they don't -- must
        // count edges per component, not nodes.
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
        // second pass, after all unions have settled: tally edges per final component root.
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
      // partial selection sort: true globally-nearest neighbors first (inner loop scans the full
      // jnum, not just the first n) -- same technique as compute_ptm.cu, same reasoning: CNA's
      // own "which neighbors are the nn nearest" question needs the true nearest set, not
      // whichever happen to land first in the buffer's own traversal order.
      const int n = std::min( jnum, CNA_MAX_NBRS );
      for(int i=0;i<n;i++)
      {
        int m = i;
        for(int j=i+1;j<jnum;j++) { if( buf.d2[j] < buf.d2[m] ) { m = j; } }
        if( m != i )
        {
          std::swap( buf.drx[i], buf.drx[m] );
          std::swap( buf.dry[i], buf.dry[m] );
          std::swap( buf.drz[i], buf.drz[m] );
          std::swap( buf.d2[i], buf.d2[m] );
        }
      }

      double struct_type = classify_family( buf, n, 12, 12, 1.0 ); // FCC/HCP/ICO: mean of 12 nearest, no rescale
      if( struct_type == CNA_MATCH_NONE )
      {
        struct_type = classify_family( buf, n, 14, 8, 2.0/std::sqrt(3.0) ); // BCC: mean of nearest 8 (1st shell), rescaled to 2nd-shell length scale
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
    ADD_SLOT( double              , rcut            , INPUT , REQUIRED , DocString{"Neighbor search cutoff -- must be generous enough for at least BCC's 14 first+second-shell neighbors plus one more (the largest requirement among the checked structures: FCC/HCP/ICO need 12+1, BCC needs 14+1)"} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double              , rcut_max        , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );
    ADD_SLOT( onika::memory::CudaMMVector<double> , cna_struct_type , OUTPUT , DocString{"Flat per-particle structure type buffer (0=none, 1=FCC, 2=HCP, 3=BCC, 4=ICO -- same numbering as compute_ptm's ptm_struct_type): cna_struct_type[ cell_particle_offset[cell] + particle ], see grid->cell_particle_offset_data()"} );
    ADD_SLOT( long , cna_n_matched , OUTPUT , DocString{"Number of particles CNA found a match for (out of grid->number_of_particles())"} );

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

      using ComputeBuffer = ComputePairBuffer2<false,false>;
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
                                         , compute_field_set, onika::UIntConst<1>{}, cp_opts
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

Per-particle structure identification via adaptive-cutoff Common Neighbor Analysis (CNA),
reimplemented from scratch following OVITO's own structure-identification engine (the exact
classifier its DXA modifier uses under the hood, per OVITO's own user manual). Checks FCC, HCP,
BCC and ICO (not SC or diamond -- see compute_cna.cu). Writes a single flat per-particle buffer,
indexed via grid->cell_particle_offset_data() (same convention as compute_ptm's own
ptm_struct_type, including its 0=none/1=FCC/2=HCP/3=BCC/4=ICO numbering) -- can be used as a
drop-in replacement for compute_dxa_edge_vectors' struct_field, mixed with compute_ptm/ptm_fields'
own ptm_orientation for orientation (CNA itself, per this reference implementation, only produces
a per-cluster orientation, not per-atom).

Usage example:

compute_cna: { rcut: 5.0 ang }

)EOF";
    }
  };

  // Materializes compute_cna's flat CudaMMVector buffer into a named per-particle grid field
  // (field::mk_generic_real), so any consumer expecting a real grid field (compute_dxa_edge_
  // vectors' struct_field, write_xyz, write_delaunay_vtk's color_fields, ...) can read it
  // directly. Purely pointwise, no neighbor search -- same exanb::compute_cell_particles pattern
  // as ptm_fields.
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

Copies compute_cna's flat per-particle output buffer (cna_struct_type) into a named per-particle
grid field, so any consumer expecting a real grid field (compute_dxa_edge_vectors' struct_field,
write_xyz, write_delaunay_vtk's color_fields, ...) can read it directly. Purely pointwise (no
neighbor search) -- run any time after compute_cna. Only touches owned particles -- ghost copies
of this field are not synchronized across MPI ranks.

Usage example:

compute_cna: { rcut: 5.0 ang }
cna_fields: { struct_field: cna_type }

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_cna)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_cna", make_grid_variant_operator< ComputeCNA > );
    OperatorNodeFactory::instance()->register_factory( "cna_fields", make_grid_variant_operator< CNAFields > );
  }

}
