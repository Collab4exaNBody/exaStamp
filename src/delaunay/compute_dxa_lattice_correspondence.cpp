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

#include <exaStamp/delaunay/lattice_structure.h>
#include <exaStamp/delaunay/dxa_lattice_correspondence.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <string>
#include <vector>

// DXA elastic-mapping step, replacing PTM's role: per atom, a DISCRETE match of its own physical
// neighbor bond graph against OVITO's real fixed reference tables (lattice_structure.h), ported
// from StructureAnalysis::determineLocalStructure() (ovito/src/ovito/crystalanalysis/modifier/dxa/
// StructureAnalysis.cpp) -- see src/delaunay/README.md for why this replaces compute_ptm's
// continuous RMSD orientation fit specifically for the DXA pipeline (PTM's own strain tolerance
// turned out to be materially wider than the discrete topology test OVITO's real DXA source uses).
//
// Two-stage match, same structure as compute_cna.cu's own adaptive-CNA classification (reused
// verbatim here for the cutoff/bond-graph/signature machinery -- see that file's comment for the
// formulas), plus a NEW second stage this operator adds: a backtracking permutation search
// (`match_permutation`, ported from determineLocalStructure()'s own "find first matching neighbor
// permutation" loop) that finds a bijection between the atom's own nn physical neighbors and the
// reference structure's nn canonical slots, consistent both in per-slot CNA signature AND in
// pairwise bond topology. This is a graph-isomorphism test against a FIXED, unrotated reference
// (lattice_structure.h's own tables), not a continuous fit -- the atom's own physical rotation is
// absorbed entirely into WHICH neighbor lands in which slot, never into a rotation matrix. Output
// per atom: which physical neighbor occupies each canonical slot (DXALatticeCorrespondence::
// neighbor_atom) -- this labeling is only self-consistent so far (each atom found its own,
// independently arbitrary, matching permutation); making it globally consistent across a whole
// grain is compute_dxa_lattice_clusters' job (next stage, buildClusters/connectClusters port).
namespace exaStamp
{
  using namespace exanb;

  static constexpr int DXA_LC_MAX_NBRS = 16; // enough for BCC's 14 + 1 sanity-check neighbor, matches compute_cna.cu's CNA_MAX_NBRS

  // Ported from determineLocalStructure()'s own backtracking search (StructureAnalysis.cpp, lines
  // ~799-862): finds a permutation `mapping` such that physical neighbor mapping[slot] satisfies
  // both the reference structure's per-slot CNA signature and pairwise bond topology
  // simultaneously, for every slot. `bonded`/`cna_sig` describe the nn physical neighbor candidates
  // (already sorted nearest-first, but that order is otherwise arbitrary -- the search finds
  // whichever valid relabeling comes first in permutation order). Returns false if no permutation
  // satisfies the reference topology (shouldn't happen once the aggregate CNA family counts already
  // matched, per OVITO's own "this should not happen" assertion, but the search terminates cleanly
  // either way since next_permutation() exhausts).
  static bool match_permutation( const LatticeStructure& ref, int nn,
                                  const bool bonded[][DXA_LC_MAX_NBRS], const int* cna_sig,
                                  int* mapping_out )
  {
    std::array<int, DXA_LC_MAX_NBRS> perm{}, prev{};
    for(int n=0;n<nn;n++) { perm[n]=n; prev[n]=-1; }

    for(;;)
    {
      int ni1 = 0;
      while( ni1<nn && perm[ni1]==prev[ni1] ) { ++ni1; }
      for(; ni1<nn; ni1++)
      {
        const int a1 = perm[ni1];
        prev[ni1] = a1;
        if( cna_sig[a1] != ref.cna_signature[ni1] ) { break; }
        int ni2;
        for(ni2=0; ni2<ni1; ni2++)
        {
          const int a2 = perm[ni2];
          if( bonded[a1][a2] != ref.is_bonded(ni1,ni2) ) { break; }
        }
        if( ni2 != ni1 ) { break; }
      }

      if( ni1 == nn )
      {
        for(int i=0;i<nn;i++) { mapping_out[i] = perm[i]; }
        return true;
      }

      std::vector<int> suffix( perm.begin(), perm.begin()+nn );
      dxa_bitmap_sort_desc( suffix, ni1+1, nn, nn );
      std::copy( suffix.begin(), suffix.end(), perm.begin() );
      if( !std::next_permutation( perm.begin(), perm.begin()+nn ) ) { return false; }
    }
  }

  struct LatticeCorrespondenceFunctor
  {
    const size_t * const __restrict__ m_cell_particle_offset = nullptr;
    uint8_t * const __restrict__ m_structure_type_out = nullptr;
    int64_t * const __restrict__ m_neighbor_atom_out = nullptr; // atom*DXA_MAX_NEIGHBORS + slot
    std::atomic<long>* m_n_cna_ok_perm_fail = nullptr; // atoms whose aggregate CNA count matched BCC but the exact bond-topology permutation search still failed -- see execute()'s own log line comment
    LatticeStructureType m_target_structure = LATTICE_BCC; // the ONE structure this classifier ever tests -- see the "try_fcc_hcp"/"try_bcc" gating below for why

    // adaptive-cutoff common-neighbor bond graph + per-slot CNA family signature, identical
    // formulas to compute_cna.cu's own classify_family() -- see that file's header comment.
    // Returns the CNA family match (0=none, else the aggregate count pattern accepted) and fills
    // bonded[][]/cna_sig[] for the nn candidates regardless of whether the aggregate pattern is
    // ultimately accepted by the caller.
    template<class ComputeBufferT>
    static inline bool compute_topology( const ComputeBufferT& buf, int jnum, int nn, int scaling_nbrs, double scaling_rescale,
                                          bool bonded[][DXA_LC_MAX_NBRS], int* cna_sig )
    {
      if( jnum < nn+1 ) { return false; }

      double scaling = 0.0;
      for(int i=0;i<scaling_nbrs;i++) { scaling += std::sqrt( buf.d2[i] ); }
      scaling /= static_cast<double>(scaling_nbrs);
      const double cutoff = scaling * scaling_rescale * (1.0+std::sqrt(2.0)) * 0.5;
      const double cutoff2 = cutoff*cutoff;

      if( buf.d2[nn] <= cutoff2 ) { return false; }

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

      for(int i=0;i<nn;i++)
      {
        int common[DXA_LC_MAX_NBRS]; int nc=0;
        for(int j=0;j<nn;j++) { if( j!=i && bonded[i][j] ) { common[nc++] = j; } }

        int numNeighborBonds = 0;
        int parent[DXA_LC_MAX_NBRS];
        for(int k=0;k<nc;k++) { parent[k] = k; }
        struct { int* p; int operator()(int x){ while(p[x]!=x){p[x]=p[p[x]];x=p[x];} return x; } } find{parent};
        for(int a=0;a<nc;a++) for(int b=a+1;b<nc;b++)
        {
          if( bonded[ common[a] ][ common[b] ] ) { ++numNeighborBonds; const int ra=find(a), rb=find(b); if(ra!=rb){parent[ra]=rb;} }
        }
        int edgecount[DXA_LC_MAX_NBRS] = {0};
        for(int a=0;a<nc;a++) for(int b=a+1;b<nc;b++)
        {
          if( bonded[ common[a] ][ common[b] ] ) { ++edgecount[ find(a) ]; }
        }
        int maxChainLength = 0;
        for(int k=0;k<nc;k++) { maxChainLength = std::max( maxChainLength, edgecount[k] ); }

        cna_sig[i] = -1;
        if     ( nc==4 && numNeighborBonds==2 && maxChainLength==1 ) { cna_sig[i] = 0; } // 421
        else if( nc==4 && numNeighborBonds==2 && maxChainLength==2 ) { cna_sig[i] = 1; } // 422
        else if( nc==4 && numNeighborBonds==4 && maxChainLength==4 ) { cna_sig[i] = 1; } // 444
        else if( nc==6 && numNeighborBonds==6 && maxChainLength==6 ) { cna_sig[i] = 0; } // 666
      }
      return true;
    }

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator () ( int jnum, ComputeBufferT& buf, CellParticlesT /*cells*/ ) const
    {
      // partial selection sort (true nearest-first, full jnum scan -- same as compute_cna.cu/
      // compute_ptm.cu), swapping neighbor identity along with position so buf.nbh stays aligned.
      const int n = std::min( jnum, DXA_LC_MAX_NBRS );
      for(int i=0;i<n;i++)
      {
        int m = i;
        for(int j=i+1;j<jnum;j++) { if( buf.d2[j] < buf.d2[m] ) { m = j; } }
        if( m != i )
        {
          std::swap( buf.drx[i], buf.drx[m] ); std::swap( buf.dry[i], buf.dry[m] );
          std::swap( buf.drz[i], buf.drz[m] ); std::swap( buf.d2[i], buf.d2[m] );
          size_t c1,p1,c2,p2; buf.nbh.get(i,c1,p1); buf.nbh.get(m,c2,p2);
          buf.nbh.set(i,c2,p2); buf.nbh.set(m,c1,p1);
        }
      }

      const size_t i = m_cell_particle_offset[buf.cell] + buf.part;
      m_structure_type_out[i] = static_cast<uint8_t>(LATTICE_OTHER);
      for(int s=0;s<DXA_MAX_NEIGHBORS;s++) { m_neighbor_atom_out[ i*DXA_MAX_NEIGHBORS + s ] = -1; }

      bool bonded12[12][DXA_LC_MAX_NBRS]; int sig12[DXA_LC_MAX_NBRS];
      bool bonded14[14][DXA_LC_MAX_NBRS]; int sig14[DXA_LC_MAX_NBRS];

      LatticeStructureType matched = LATTICE_OTHER;
      int nn = 0;
      int mapping[DXA_MAX_NEIGHBORS];

      // Only ever test the ONE user-declared target structure -- matching OVITO's real
      // DislocationAnalysisModifier exactly (its own input_crystal_structure is a single required
      // choice, never "try everything and see what sticks"). Originally this checked FCC/HCP THEN
      // BCC unconditionally for every atom regardless of target: found (via a real, measured
      // discrepancy against OVITO's own DXA-internal classification, see this file's own commit
      // history / src/delaunay/README.md) that this lets a handful of atoms right at a genuinely
      // disordered region -- where local coordination is distorted enough to coincidentally,
      // spuriously satisfy the FCC 12-neighbor topology test purely by combinatorial chance --
      // get misclassified as FCC in an otherwise pure-BCC system, isolating them into their own
      // spurious FCC-type cluster with no valid transition to the real surrounding BCC cluster.
      // Only 2 atoms out of 128000 on the real quadrupole test case, but both sat within ~3.5 Å of
      // the exact same real 3-way dislocation junction -- exactly where a locally different
      // interface-mesh shape has an outsized effect on which arm's seed search wins a territorial
      // race during the later circuit sweep. A tiny classification bug with a locally
      // disproportionate downstream effect.
      const bool try_fcc_hcp = ( m_target_structure == LATTICE_FCC || m_target_structure == LATTICE_HCP );
      const bool try_bcc = ( m_target_structure == LATTICE_BCC );

      if( try_fcc_hcp && compute_topology( buf, n, 12, 12, 1.0, bonded12, sig12 ) )
      {
        int n421=0, n422=0;
        for(int k=0;k<12;k++) { if(sig12[k]==0) ++n421; else if(sig12[k]==1) ++n422; }
        if( m_target_structure == LATTICE_FCC && n421==12 )
        {
          const LatticeStructure& ref = dxa_lattice_structure(LATTICE_FCC);
          if( match_permutation( ref, 12, bonded12, sig12, mapping ) ) { matched = LATTICE_FCC; nn = 12; }
        }
        else if( m_target_structure == LATTICE_HCP && n421==6 && n422==6 )
        {
          const LatticeStructure& ref = dxa_lattice_structure(LATTICE_HCP);
          if( match_permutation( ref, 12, bonded12, sig12, mapping ) ) { matched = LATTICE_HCP; nn = 12; }
        }
      }

      if( try_bcc && matched == LATTICE_OTHER && compute_topology( buf, n, 14, 8, 2.0/std::sqrt(3.0), bonded14, sig14 ) )
      {
        int n444=0, n666=0;
        for(int k=0;k<14;k++) { if(sig14[k]==1) ++n444; else if(sig14[k]==0) ++n666; }
        if( n666==8 && n444==6 )
        {
          const LatticeStructure& ref = dxa_lattice_structure(LATTICE_BCC);
          if( match_permutation( ref, 14, bonded14, sig14, mapping ) ) { matched = LATTICE_BCC; nn = 14; }
          else if( m_n_cna_ok_perm_fail ) { m_n_cna_ok_perm_fail->fetch_add(1, std::memory_order_relaxed); }
        }
      }

      if( matched != LATTICE_OTHER )
      {
        m_structure_type_out[i] = static_cast<uint8_t>(matched);
        for(int slot=0; slot<nn; slot++)
        {
          size_t c=0, p=0; buf.nbh.get( mapping[slot], c, p );
          m_neighbor_atom_out[ i*DXA_MAX_NEIGHBORS + slot ] = static_cast<int64_t>( m_cell_particle_offset[c] + p );
        }
      }
    }
  };

  template<class GridT>
  class ComputeDXALatticeCorrespondence : public OperatorNode
  {
    ADD_SLOT( GridT               , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain              , domain          , INPUT , REQUIRED );
    ADD_SLOT( double              , rcut            , INPUT , REQUIRED , DocString{"Neighbor search cutoff -- must be generous enough for at least BCC's 14 first+second-shell neighbors plus one more"} );
    ADD_SLOT( std::string         , target_structure , INPUT , std::string("BCC") , DocString{"The ONE crystal structure to classify against: FCC, HCP or BCC -- matching OVITO's own DislocationAnalysisModifier.input_crystal_structure exactly (a single required choice, never multiple structures tried per atom). See this operator's own functor comment for why testing more than the declared target caused a real, measured bug."} );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( double              , rcut_max        , INPUT_OUTPUT , 0.0 , DocString{"Updated max rcut"} );
    ADD_SLOT( DXALatticeCorrespondence , dxa_lattice_correspondence , OUTPUT );
    ADD_SLOT( long                , n_matched       , OUTPUT , DocString{"Number of particles matched to a discrete lattice correspondence (out of grid->number_of_particles())"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      *rcut_max = std::max( *rcut , *rcut_max );
      if( grid->number_of_cells() == 0 ) return;

      if( chunk_neighbors->m_chunk_size != 1 )
      {
        fatal_error() << "compute_dxa_lattice_correspondence: requires chunk_neighbors built with chunk_size=1 (set random_access: true"
                       << " on the chunk_neighbors config), got chunk_size=" << chunk_neighbors->m_chunk_size << std::endl;
      }

      LatticeStructureType target = LATTICE_BCC;
      if( *target_structure == "FCC" ) { target = LATTICE_FCC; }
      else if( *target_structure == "HCP" ) { target = LATTICE_HCP; }
      else if( *target_structure == "BCC" ) { target = LATTICE_BCC; }
      else { fatal_error() << "compute_dxa_lattice_correspondence: unknown target_structure '" << *target_structure << "', expected one of FCC, HCP, BCC" << std::endl; }

      const size_t total_particles = grid->number_of_particles();
      DXALatticeCorrespondence& result = *dxa_lattice_correspondence;
      result.structure_type.assign( total_particles, static_cast<uint8_t>(LATTICE_OTHER) );
      result.neighbor_atom.assign( total_particles * DXA_MAX_NEIGHBORS, -1 );

      using ComputeBuffer = ComputePairBuffer2<false,true>; // UseNeighbors=true: buf.nbh.get/set per neighbor slot
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto compute_buf = make_compute_pair_buffer<ComputeBuffer>();

      std::atomic<long> n_cna_ok_perm_fail{0};
      LatticeCorrespondenceFunctor compute_op = { grid->cell_particle_offset_data(), result.structure_type.data(), result.neighbor_atom.data(), &n_cna_ok_perm_fail, target };

      ComputePairNullWeightIterator cp_weight{};
      LinearXForm cp_xform { domain->xform() };
      auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks );
      static constexpr onika::FlatTuple<> compute_field_set = {};
      static constexpr DefaultPositionFields posfields = {};
      static constexpr ComputeParticlePairOpts<false,true,false> cp_opts = {}; // Symmetric=false, PreferComputeBuffer=true

      const IJK dims = grid->dimension();
      const auto cells = grid->cells();
      const double rcut2 = (*rcut) * (*rcut);

      // Unlike compute_ptm/compute_cna (which only classify OWNED cells, relying on
      // ghost_update_opt to sync a named field for the few consumers that need ghost values),
      // this classifier covers the FULL grid, ghost cells included. Downstream consumers here
      // (compute_dxa_crystal_path_edge_vectors, later the per-tet elastic-mapping test) need a
      // real classification on ghost-copy ATOMS of a Delaunay tessellation vertex, not just owned
      // ones -- and unlike PTM/CNA's flat output there's no per-particle *named grid field* to hang
      // a ghost_update_opt sync off of anyway (DXALatticeCorrespondence is its own struct).
      // Measured: this raises edge resolution on the real BCC Ta test system from 92.8% (16000/
      // 41635 particles matched, owned-only) to 100% (33201/41635 matched, all cells) -- confirming
      // chunk_neighbors does carry valid neighbor data for ghost-centered cells on this build.
      // Caveat: only exercised so far on a single-MPI-rank periodic system, where every "ghost" is
      // a periodic self-image with a fully populated local neighborhood (ghost_dist_max is 2x rcut
      // here) -- a genuine cross-rank ghost near the outer edge of a real multi-rank halo could
      // still see a truncated neighbor count and be silently (correctly) left unmatched, same
      // "ghost-fringe-trust" caveat compute_delaunay.cpp already documents for tets. Revisit if
      // this pipeline is ever run multi-rank.
#     pragma omp parallel for collapse(3) schedule(dynamic)
      for(ssize_t k=0;k<dims.k;k++)
      for(ssize_t j=0;j<dims.j;j++)
      for(ssize_t i=0;i<dims.i;i++)
      {
        const IJK cell_a_loc{i,j,k};
        const size_t cell_a = static_cast<size_t>( grid_ijk_to_index(dims,cell_a_loc) );
        compute_cell_particle_pairs_cell( cells, dims, cell_a_loc, cell_a, rcut2
                                         , compute_buf, optional, compute_op
                                         , compute_field_set, onika::UIntConst<1>{}, cp_opts
                                         , posfields, std::index_sequence<>{} );
      }

      size_t n = 0;
      for(size_t p=0;p<total_particles;p++) { if( result.structure_type[p] != static_cast<uint8_t>(LATTICE_OTHER) ) { ++n; } }
      *n_matched = static_cast<long>(n);

      // Owned-only matched count, reported alongside the full (owned+ghost) count above -- lets a
      // reader compare directly against compute_cna's own owned-only convention (compute_cna never
      // touches ghost cells at all, so its "matched/total" ratio silently counts every ghost slot
      // as an automatic reject; comparing the two operators' raw matched/total numbers directly is
      // misleading without this). Confirmed once on the real quadrupole dislocation case: this
      // classifier's own owned-only rejection (128000-126976=1024) is IDENTICAL to compute_cna's
      // owned-only rejection on the same file -- the extra bond-topology permutation-match stage
      // (n_cna_ok_perm_fail below) never rejects anything CNA's own aggregate count wouldn't have
      // already rejected, on this dataset.
      const ssize_t gl = grid->ghost_layers();
      const size_t * const cell_particle_offset = grid->cell_particle_offset_data();
      size_t n_owned_matched = 0, n_owned_total = 0;
      for(ssize_t k=gl;k<dims.k-gl;k++)
      for(ssize_t j=gl;j<dims.j-gl;j++)
      for(ssize_t i=gl;i<dims.i-gl;i++)
      {
        const size_t c = static_cast<size_t>( grid_ijk_to_index(dims, IJK{i,j,k}) );
        const size_t np = cells[c].size();
        for(size_t p=0;p<np;p++)
        {
          const size_t idx = cell_particle_offset[c] + p;
          ++n_owned_total;
          if( result.structure_type[idx] != static_cast<uint8_t>(LATTICE_OTHER) ) { ++n_owned_matched; }
        }
      }

      lout << "compute_dxa_lattice_correspondence: " << n << " / " << total_particles << " particles matched"
           << " (owned-only: " << n_owned_matched << "/" << n_owned_total
           << ", " << n_cna_ok_perm_fail.load() << " passed the aggregate CNA count but failed the exact bond-topology permutation match)" << std::endl;
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

DXA elastic-mapping classifier, replacing compute_ptm's role: per atom, a DISCRETE match of its own
physical neighbor bond graph against fixed reference tables (FCC/HCP/BCC), ported from OVITO's real
DXA source (StructureAnalysis::determineLocalStructure) -- see this file's own header comment. Not a
continuous RMSD/orientation fit: the atom's rotation is captured entirely by which physical neighbor
occupies which canonical slot.

Usage example:

chunk_neighbors: { config: { chunk_size: 1 } }
compute_dxa_lattice_correspondence: { rcut: 5.0 ang }
compute_dxa_lattice_clusters: {}

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_dxa_lattice_correspondence)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_dxa_lattice_correspondence", make_grid_variant_operator< ComputeDXALatticeCorrespondence > );
  }

}
