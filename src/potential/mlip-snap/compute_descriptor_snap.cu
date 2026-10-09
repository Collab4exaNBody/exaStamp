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

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/core/concurent_add_contributions.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/file_utils.h>
#include <onika/cuda/cuda_math.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <md/snap/snap_config.h>
#include <md/snap/snap_context.h>
#include <md/snap/snap_compute_buffer.h>
#include <md/snap/snap_bispectrum_op.h>
#include <md/snap/snap_compute_dbidrj.h>
#include <md/snap/snap_check_bispectrum.h>
#include <md/snap/sna.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <mpi.h>

// GPU-compatible bispectrum descriptor pass: same SNA setup and BispectrumOpRealT compute-pair
// call as md::SnapForceRealT (snap_force.h), without the force/energy pass. BispectrumOpRealT only
// writes to its central particle's slot, so no particle locks are needed, and only radelem/wjelem
// and the SNA configuration are used (the SNAP coefficients are not). With compute_derivative, a
// CSR-shaped output holds the per-neighbor-pair bispectrum Jacobian, computed by
// md::snap_compute_neighbor_dbidrj (mono-element) or md::snap_compute_neighbor_dbidrj_multi
// (chem_flag, see snap_compute_dbidrj.h).
namespace exaStamp
{
  using namespace exanb;
  using namespace md; // snap_compute_ui/zi/bi, snap_compute_sfac, SnapBSExtStorage, ...

  struct ResetBispectrumCPBuf
  {
    template<class CPBufT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () (CPBufT & buf) const { buf.ext.reset(); }
  };

  // Bispectrum with an adaptive, per-particle cutoff chosen so each particle sees
  // (approximately) a constant number of neighbors, instead of a fixed physical radius --
  // mirrors compute_local_structural_metrics.h's "closest_bispectrum" / density-scaled
  // constant-neighbor-count modes. md::BispectrumOpRealT's own neighbor cutoff is a single
  // (radelem[i]+radelem[j])*rcutfac scalar shared by the whole grid pass, with no per-particle
  // hook, so this is a copy of its operator() (snap_bispectrum_op.h) with the neighbor-cutoff
  // step replaced by either an exact Nth-nearest-neighbor-distance cutoff (sort-based, "closest"
  // mode) or a density-scaled estimate (no sort needed, assumes uniform local density), reusing
  // exaNBody's underlying GPU-tagged SNA math (snap_compute_ui/zi/bi) unchanged.
  template<class RealT, class RijRealT, class SnapConfParamT>
  struct ConstNeighBispectrumOpRealT
  {
    const SnapConfParamT snaconf;
    const size_t * const __restrict__ cell_particle_offset = nullptr;
    RealT * const __restrict__ bispectrum = nullptr;
    const RealT * const __restrict__ wjelem = nullptr; // per-species weight (m_factor)
    const long ncoeff = 0;
    const int nneigh = 0;                                 // target (constant) neighbor count
    const RijRealT margin = static_cast<RijRealT>(0.01);  // "closest" mode: added past the Nth-neighbor distance
    const RijRealT outer_rcut = 0.0;                      // physical search radius neighbors were gathered within
    const bool closest = true;                            // true: exact Nth-nearest (sort); false: density-scaled estimate

    // descriptor-derivative pass -- see BispectrumOpRealT (snap_bispectrum_op.h) for the
    // two-dispatch (count then fill) CSR scheme this mirrors.
    long   * const __restrict__ deriv_row_offset = nullptr;
    double * const __restrict__ deriv_buffer      = nullptr;
    uint64_t * const __restrict__ deriv_nbh_id    = nullptr;
    const bool count_only_pass    = false;
    const bool compute_derivative = false;
    double * const * const __restrict__ deriv_agg_ptrs = nullptr; // optional: see BispectrumOpRealT (snap_bispectrum_op.h)

    template<class ComputeBufferT>
    ONIKA_HOST_DEVICE_FUNC
    static inline void swap_neighbor( ComputeBufferT& buf, int i, int j )
    {
      double t;
      t=buf.drx[i]; buf.drx[i]=buf.drx[j]; buf.drx[j]=t;
      t=buf.dry[i]; buf.dry[i]=buf.dry[j]; buf.dry[j]=t;
      t=buf.drz[i]; buf.drz[i]=buf.drz[j]; buf.drz[j]=t;
      t=buf.d2[i];  buf.d2[i] =buf.d2[j];  buf.d2[j] =t;
      auto tpt = buf.nbh_pt[i]; buf.nbh_pt[i]=buf.nbh_pt[j]; buf.nbh_pt[j]=tpt;
    }

    template<class ComputeBufferT, class CellParticlesT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( int jnum, ComputeBufferT& buf, int itype, CellParticlesT cells) const
    {
      assert( ncoeff == static_cast<unsigned int>(snaconf.ncoeff) );
      buf.ext.init( snaconf );

      int ninside = 0;
      RijRealT rcutij = margin;

      if( closest )
      {
        // selection sort by squared distance -- jnum is small (bounded by MAX_PARTICLE_NEIGHBORS)
        for(int i=0;i<jnum-1;i++)
        {
          int m = i;
          for(int j=i+1;j<jnum;j++) if( buf.d2[j] < buf.d2[m] ) m = j;
          if( m != i ) { swap_neighbor(buf,i,m); }
        }
        ninside = (jnum < nneigh) ? jnum : nneigh;
        rcutij = ( ninside > 0 ) ? ( sqrt(buf.d2[ninside-1]) + margin ) : margin;
      }
      else
      {
        // density-scaled estimate: assumes neighbor count scales as r^3 with uniform density,
        // so no sort is needed -- just shrink outer_rcut by (nneigh/jnum)^(1/3), capped at outer_rcut.
        const RijRealT n_ratio = ( jnum > 0 ) ? ( static_cast<RijRealT>(nneigh) / static_cast<RijRealT>(jnum) ) : static_cast<RijRealT>(0);
        rcutij = ( jnum > 0 ) ? onika::cuda::min( outer_rcut , outer_rcut * cbrt(n_ratio) ) : static_cast<RijRealT>(0);
        const RijRealT rcutij_sq = rcutij * rcutij;
        for(int jj=0;jj<jnum;jj++)
        {
          if( buf.d2[jj] < rcutij_sq ) { if( ninside != jj ) { buf.copy(jj,ninside); } ninside++; }
        }
      }

      if( count_only_pass )
      {
        deriv_row_offset[ cell_particle_offset[buf.cell] + buf.part ] = ninside + 1;
        return;
      }

      /************ begin of UiTot computation ******************/
      snap_uarraytot_zero( snaconf.nelements, snaconf.idxu_max, buf.ext.m_UTot_array.r(), buf.ext.m_UTot_array.i() );
      snap_uarraytot_init_wself( snaconf.nelements, snaconf.twojmax, snaconf.idxu_max, snaconf.wself, snaconf.wselfall_flag, buf.ext.m_UTot_array.r(), buf.ext.m_UTot_array.i(), snaconf.chem_flag ? itype : 0 );
      for (int jj = 0; jj < ninside; jj++)
      {
        const int jtype = buf.nbh_pt[jj][field::type];
        const int jelem = snaconf.chem_flag ? jtype : 0 ;

        const RijRealT x = buf.drx[jj];
        const RijRealT y = buf.dry[jj];
        const RijRealT z = buf.drz[jj];
        const RijRealT rsq = buf.d2[jj];
        const RijRealT r = sqrt(rsq);
        const RijRealT theta0 = (r - snaconf.rmin0) * snaconf.rfac0 * M_PI / (rcutij - snaconf.rmin0);
        const RijRealT z0 = r / tan(theta0);

        const RijRealT wj_jj = wjelem[jtype];
        const RijRealT sfac_jj = snap_compute_sfac( static_cast<RijRealT>(snaconf.rmin0), snaconf.switch_flag, false, r, rcutij, static_cast<RijRealT>(0), static_cast<RijRealT>(0) );

        snap_add_nbh_contrib_to_uarraytot( snaconf.twojmax, sfac_jj*wj_jj, x,y,z,z0,r, snaconf.rootpqarray, buf.ext.m_UTot_array.r() + snaconf.idxu_max * jelem, buf.ext.m_UTot_array.i() + snaconf.idxu_max * jelem, buf.ext );
      }
      /****************** end of UiTot computation **********************/

      snap_compute_zi( snaconf.nelements, snaconf.idxz_max, snaconf.idxu_max, snaconf.twojmax
                     , snaconf.idxcg_block
                     , snaconf.idxz, snaconf.cglist, buf.ext.m_UTot_array.r(), buf.ext.m_UTot_array.i()
                     , snaconf.bnorm_flag, buf.ext.m_Z_array.r(), buf.ext.m_Z_array.i() );

      const long bispectrum_ii_offset = snaconf.ncoeff * ( cell_particle_offset[buf.cell] + buf.part );
      snap_compute_bi( snaconf.nelements, snaconf.idxz_max, snaconf.idxb_max, snaconf.idxu_max, snaconf.twojmax
                     , snaconf.idxz_block
                     , snaconf.idxz, snaconf.idxb
                     , buf.ext.m_Z_array.r(), buf.ext.m_Z_array.i()
                     , buf.ext.m_UTot_array.r(), buf.ext.m_UTot_array.i()
                     , snaconf.bzero , snaconf.bzero_flag, snaconf.wselfall_flag
                     , bispectrum + bispectrum_ii_offset
                     , snaconf.chem_flag ? itype : 0 );

      if( compute_derivative )
      {
        const long p = cell_particle_offset[buf.cell] + buf.part;
        const long row0 = deriv_row_offset[p];
        const int ncoeff3 = static_cast<int>(snaconf.ncoeff) * 3;
        // multi-type widening: this pair's contributions (both sides) land in the CENTRAL atom's
        // (itype's) own ncoeff3-wide slot -- see snap_bispectrum_op.h's deriv_agg_ptrs comment.
        const int typeoffset = ncoeff3 * itype;

        double * const __restrict__ self_row = deriv_buffer + row0 * ncoeff3;
        for( int i=0; i<ncoeff3; i++ ) self_row[i] = 0.0;
        deriv_nbh_id[row0] = cells[buf.cell][field::id][buf.part];

        for (int jj = 0; jj < ninside; jj++)
        {
          const int jtype = buf.nbh_pt[jj][field::type];

          const RijRealT x = buf.drx[jj];
          const RijRealT y = buf.dry[jj];
          const RijRealT z = buf.drz[jj];
          const RijRealT rsq = buf.d2[jj];
          const RijRealT r = sqrt(rsq);
          const RijRealT theta0 = (r - snaconf.rmin0) * snaconf.rfac0 * M_PI / (rcutij - snaconf.rmin0);
          const RijRealT z0 = r / tan(theta0);
          const RijRealT wj_jj = wjelem[jtype];

          double * const __restrict__ row = deriv_buffer + (row0 + 1 + jj) * ncoeff3;
          if( snaconf.chem_flag )
          {
            const int jelem = jtype;
            snap_compute_neighbor_dbidrj_multi( snaconf.twojmax, snaconf.idxu_max, snaconf.idxb_max, snaconf.idxz_max, snaconf.nelements, jelem,
                                                 wj_jj, rcutij, static_cast<RijRealT>(0), static_cast<RijRealT>(0), x, y, z, z0, r,
                                                 snaconf.rootpqarray, snaconf.idxz_block, snaconf.idxb,
                                                 buf.ext.m_Z_array.r(), buf.ext.m_Z_array.i(),
                                                 snaconf.rmin0, snaconf.rfac0, snaconf.switch_flag, false, snaconf.bnorm_flag,
                                                 row, buf.ext );
          }
          else
          {
            snap_compute_neighbor_dbidrj( snaconf.twojmax, snaconf.idxu_max, snaconf.idxb_max,
                                           wj_jj, rcutij, static_cast<RijRealT>(0), static_cast<RijRealT>(0), x, y, z, z0, r,
                                           snaconf.rootpqarray, snaconf.idxz_block, snaconf.idxb,
                                           buf.ext.m_Z_array.r(), buf.ext.m_Z_array.i(),
                                           snaconf.rmin0, snaconf.rfac0, snaconf.switch_flag, false, snaconf.bnorm_flag,
                                           row, buf.ext );
          }
          size_t nbh_cell=0, nbh_part=0;
          buf.nbh.get( jj, nbh_cell, nbh_part );
          deriv_nbh_id[row0 + 1 + jj] = cells[nbh_cell][field::id][nbh_part];
          for( int i=0; i<ncoeff3; i++ ) self_row[i] -= row[i];

          if( deriv_agg_ptrs != nullptr )
          {
            const size_t nbh_p = cell_particle_offset[nbh_cell] + nbh_part;
            for( int i=0; i<ncoeff3; i++ ) atomic_add_contribution( deriv_agg_ptrs[typeoffset+i][nbh_p], -row[i] );
          }
        }

        if( deriv_agg_ptrs != nullptr )
        {
          for( int i=0; i<ncoeff3; i++ ) atomic_add_contribution( deriv_agg_ptrs[typeoffset+i][p], -self_row[i] );
        }
      }
    }
  };
}

namespace exanb
{
  template<class RealT, class RijRealT, class SnapConfParamT>
  struct ComputePairTraits< exaStamp::ConstNeighBispectrumOpRealT<RealT,RijRealT,SnapConfParamT> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool ComputeBufferCompatible      = true;
    static inline constexpr bool BufferLessCompatible         = false;
    static inline constexpr bool CudaCompatible               = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDescriptorSnap : public OperatorNode
  {
    using RealT = double;
    using SnapContext = md::SnapXSContextRealT<RealT>;
    template<int I> using ICST = onika::IntConst<I>;
    template<int jm> using ROParamsMonoElem = SnapInternal::ReadOnlySnapParametersRealT<RealT,ICST<jm>,ICST<1>,false>;

    ADD_SLOT( MPI_Comm                 , mpi               , INPUT , REQUIRED );
    ADD_SLOT( double                   , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors  , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                     , ghost             , INPUT , false );
    ADD_SLOT( GridT                    , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                   , domain            , INPUT , REQUIRED );
    ADD_SLOT( long                     , timestep          , INPUT , REQUIRED );
    ADD_SLOT( std::string              , bispectrumchkfile , INPUT , OPTIONAL , DocString{"file with reference values to check bispectrum correctness"} );
    ADD_SLOT( double                   , check_bs_max_error, INPUT , 1.e-12 );
    ADD_SLOT( long                     , nneigh_bispectrum , INPUT , 0 , DocString{"If >0, switches to constant-neighbor-count mode: for each particle, the cutoff is adapted so (approximately) this many neighbors contribute, instead of the fixed physical rcut from the SNAP param file. 0 (default) keeps the fixed-rcut behavior. The param file's rcutfac must still be generous enough to gather at least this many neighbors."} );
    ADD_SLOT( double                   , neigh_margin      , INPUT , 0.01 , DocString{"closest_bispectrum=true only: small margin added past the Nth-nearest-neighbor distance, in the same length unit as positions, so that neighbor sits inside rather than exactly on the smooth cutoff's edge."} );
    ADD_SLOT( bool                     , closest_bispectrum, INPUT , true , DocString{"When nneigh_bispectrum>0: true (default) picks the exact Nth-nearest-neighbor distance as cutoff (sort-based, exact but O(n^2) per particle); false uses a density-scaled estimate instead -- assumes neighbor count scales as r^3 with uniform local density, so rcut = min(rcutfac, rcutfac*(nneigh_bispectrum/n_within_rcutfac)^(1/3)) with no sorting needed, cheaper but approximate."} );

    ADD_SLOT( bool                     , compute_derivative, INPUT , false , DocString{"If true, also computes the per-neighbor-pair bispectrum derivative Jacobian (mono-element and chem_flag=true multi-element SNAP configs both supported)."} );
    ADD_SLOT( std::string              , deriv_agg_field_prefix, INPUT , std::string("sda_") , DocString{"compute_derivative only: name prefix for the ncoeff*3*ntypes dynamically-named generic-real grid fields ('<prefix>0'..'<prefix>{ncoeff*3*ntypes-1}', widened by ntypes -- one ncoeff*3-wide slot per possible CENTRAL atom type, mono-type is the trivial ntypes==1 case) holding the force-signed per-atom derivative aggregate. Stored as grid fields so update_opt_from_ghost can reduce them ghost->owner: add 'update_opt_from_ghost: { opt_fields: [\"<prefix>.*\"] }' after this operator on multi-rank runs (after any *_global operator). Keep it short: field names are limited to 15 characters, this operator aborts if prefix+index is longer."} );

    ADD_SLOT( SnapContext              , snap_ctx          , INPUT , REQUIRED , DocString{"SNAP context built by snap_init"} );
    ADD_SLOT( onika::memory::CudaMMVector<RealT> , bispectrum , OUTPUT , DocString{"Flat per-particle bispectrum buffer: bispectrum[ ncoeff*(cell_particle_offset[cell]+particle) + component ], see grid->cell_particle_offset_data()"} );
    ADD_SLOT( long                     , ncoeff            , OUTPUT , DocString{"Number of bispectrum coefficients per particle (stride of the bispectrum buffer)"} );
    ADD_SLOT( onika::memory::CudaMMVector<long>    , bispectrum_deriv_offset , OUTPUT , DocString{"compute_derivative only: CSR row offset per particle (size total_particles+1); particle p's rows span [offset[p],offset[p+1]), row 0 of each particle's block is its own self/negative-sum term, rows 1..ninside are its neighbors in compacted order."} );
    ADD_SLOT( onika::memory::CudaMMVector<double>  , bispectrum_deriv        , OUTPUT , DocString{"compute_derivative only: flat per-neighbor-pair bispectrum Jacobian, bispectrum_deriv[ (bispectrum_deriv_offset[p]+row)*ncoeff*3 + k*3 + xyz ]."} );
    ADD_SLOT( onika::memory::CudaMMVector<uint64_t>, bispectrum_deriv_nbh_id , OUTPUT , DocString{"compute_derivative only: field::id of the atom each derivative row belongs to (row 0 = particle p's own id)."} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      // snap_ctx (param/coef files, per-material tables, SNA setup) is built by snap_init
      *rcut_max = std::max( double(*rcut_max) , double(snap_ctx->m_rcut) );

      if( grid->number_of_cells() == 0 ) { *ncoeff = 0; return; }

      // ncoeff (bispectrum stride) comes from twojmax/nelements (sna.h's compute_coeff_count):
      // only the per-material header lines (name, radelem, weight) of the coefficient file matter
      const bool quadraticflag = snap_ctx->m_config.quadraticflag();
      const int ncoeff_local = snap_ctx->sna->ncoeff;
      *ncoeff = ncoeff_local;

      const size_t total_particles = grid->number_of_particles();
      bispectrum->clear();
      bispectrum->resize( total_particles * ncoeff_local );
      if( ! *compute_derivative )
      {
        bispectrum_deriv_offset->clear();
        bispectrum_deriv->clear();
        bispectrum_deriv_nbh_id->clear();
      }

      ComputePairNullWeightIterator cp_weight{};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      LinearXForm cp_xform { domain->xform() };
      ComputePairOptionalLocks<false> cp_locks{};
      static constexpr FieldSet<field::_type> compute_bispectrum_field_set{};

      auto snap_compute_bispectrum_specialized = [&]( const auto & snapconf )
      {
        using SnapConfParamsT = std::remove_cv_t< std::remove_reference_t< decltype(snapconf) > >;
        using ComputeBufferBS = ComputePairBuffer2< false, true
                                        , md::SnapBSExtStorage<SnapConfParamsT,RealT>, DefaultComputePairBufferAppendFunc
                                        , exanb::MAX_PARTICLE_NEIGHBORS, ComputePairBuffer2Weights
                                        , FieldSet<field::_type> >;

        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks
                        , ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{}, grid->field_accessors_from_field_set(FieldSet<field::_type>{}) );
        auto cp_fields = grid->field_accessors_from_field_set( compute_bispectrum_field_set );

        // one dispatch of the compute-pair driver; called once (no derivative) or twice
        // (derivative: a cheap count_only pass to size the CSR buffers, then the real fill pass)
        auto dispatch = [&]( bool count_only, bool do_deriv, long* deriv_off, double* deriv_buf, uint64_t* deriv_id, double* const * deriv_agg_ptrs )
        {
          auto bs_buf = make_compute_pair_buffer< ComputeBufferBS , ResetBispectrumCPBuf >();
          if( *nneigh_bispectrum > 0 )
          {
            if( snap_ctx->m_config.switchinnerflag() )
            {
              fatal_error() << "compute_descriptor_snap: nneigh_bispectrum (constant-neighbor-count) mode does not support switchinnerflag" << std::endl;
            }
            ConstNeighBispectrumOpRealT<RealT,RealT,SnapConfParamsT> bispectrum_op {
                                 snapconf, grid->cell_particle_offset_data(), bispectrum->data(),
                                 snap_ctx->m_factor.data(), ncoeff_local, static_cast<int>(*nneigh_bispectrum), *neigh_margin,
                                 snap_ctx->m_rcut, *closest_bispectrum,
                                 deriv_off, deriv_buf, deriv_id, count_only, do_deriv, deriv_agg_ptrs };
            compute_cell_particle_pairs2( *grid, snap_ctx->m_rcut, *ghost, optional, bs_buf, bispectrum_op, cp_fields
                                         , DefaultPositionFields{}, parallel_execution_context() );
          }
          else
          {
            md::BispectrumOpRealT<RealT,RealT,SnapConfParamsT> bispectrum_op {
                                 snapconf, grid->cell_particle_offset_data(), nullptr, bispectrum->data(),
                                 nullptr, ncoeff_local, snap_ctx->m_factor.data(), snap_ctx->m_radelem.data(),
                                 // rcutfac is the per-pair scale factor (cut_ij = (radelem[i]+radelem[j])*rcutfac),
                                 // not the neighbor-list cutoff snap_ctx->m_rcut
                                 nullptr, nullptr, snap_ctx->m_config.rcutfac(), false, quadraticflag,
                                 deriv_off, deriv_buf, deriv_id, count_only, do_deriv, deriv_agg_ptrs };
            compute_cell_particle_pairs2( *grid, snap_ctx->m_rcut, *ghost, optional, bs_buf, bispectrum_op, cp_fields
                                         , DefaultPositionFields{}, parallel_execution_context() );
          }
        };

        if( *compute_derivative )
        {
          bispectrum_deriv_offset->clear();
          bispectrum_deriv_offset->resize( total_particles + 1, 0 );
          dispatch( true, false, bispectrum_deriv_offset->data(), nullptr, nullptr, nullptr );

          long total_rows = 0;
          for( size_t i=0; i<total_particles; i++ )
          {
            const long cnt = (*bispectrum_deriv_offset)[i];
            (*bispectrum_deriv_offset)[i] = total_rows;
            total_rows += cnt;
          }
          (*bispectrum_deriv_offset)[total_particles] = total_rows;

          bispectrum_deriv->clear();
          bispectrum_deriv->resize( total_rows * ncoeff_local * 3 );
          bispectrum_deriv_nbh_id->clear();
          bispectrum_deriv_nbh_id->resize( total_rows );

          // aggregate: one dynamically-named generic-real grid field per component (not a
          // private buffer) so the generic update_opt_from_ghost operator can reduce ghost
          // contributions back into their real owner across MPI ranks -- see
          // deriv_agg_field_prefix's DocString. Scattered into directly by the fill pass below;
          // this operator does NOT itself do any ghost/cross-rank reduction.
          // one ncoeff*3-wide slot per possible CENTRAL atom type (see snap_bispectrum_op.h)
          const size_t ntypes = snap_ctx->m_config.materials().size();
          const size_t nc3 = static_cast<size_t>(ncoeff_local) * 3 * ntypes;
          // onika::soatl::FieldId's dynamic-field name storage is a fixed char[16] (incl. null
          // terminator), silently strncpy-truncated -- a too-long prefix+index would alias
          // multiple components onto the same field with no error, so check instead of guessing.
          static constexpr size_t FIELD_NAME_MAX_LEN = 16;
          const size_t max_index_digits = std::to_string(nc3-1).size();
          if( deriv_agg_field_prefix->size() + max_index_digits + 1 > FIELD_NAME_MAX_LEN )
          {
            fatal_error() << "compute_descriptor_snap: deriv_agg_field_prefix '"<<*deriv_agg_field_prefix<<"' is too long -- "
                          <<"prefix ("<<deriv_agg_field_prefix->size()<<" chars) + largest index ("<<max_index_digits<<" digits) "
                          <<"+ null terminator must fit in "<<FIELD_NAME_MAX_LEN<<" characters (onika::soatl::FieldId's fixed name buffer)" << std::endl;
          }
          onika::memory::CudaMMVector<double*> deriv_agg_ptrs;
          deriv_agg_ptrs.resize( nc3 );
          for( size_t k=0; k<nc3; k++ )
          {
            double * const ptr = grid->flat_array_data( field::mk_generic_real( *deriv_agg_field_prefix + std::to_string(k) ) );
            std::fill_n( ptr, total_particles, 0.0 );
            deriv_agg_ptrs[k] = ptr;
          }

          dispatch( false, true, bispectrum_deriv_offset->data(), bispectrum_deriv->data(), bispectrum_deriv_nbh_id->data(), deriv_agg_ptrs.data() );
        }
        else
        {
          dispatch( false, false, nullptr, nullptr, nullptr, nullptr );
        }

        if( bispectrumchkfile.has_value() )
        {
          std::ostringstream oss; oss << *bispectrumchkfile << "." << *timestep;
          std::string file_name = onika::data_file_path( oss.str() );
          md::snap_check_bispectrum( *mpi, *grid, file_name, ncoeff_local, bispectrum->data(), *check_bs_max_error );
        }
      };

      bool fallback_to_generic = false;
      const int JMax = snap_ctx->sna->twojmax / 2;
      if( snap_ctx->sna->nelements == 1 )
      {
             if( JMax == 3 ) snap_compute_bispectrum_specialized( ROParamsMonoElem<3>(snap_ctx->sna) );
        else if( JMax == 4 ) snap_compute_bispectrum_specialized( ROParamsMonoElem<4>(snap_ctx->sna) );
        else fallback_to_generic = true;
      }
      else fallback_to_generic = true;

      if( fallback_to_generic )
      {
        snap_compute_bispectrum_specialized( SnapInternal::ReadOnlySnapParametersRealT<RealT,int,int,false>( snap_ctx->sna ) );
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Per-atom SNAP bispectrum descriptors (GPU-compatible). Uses the SNAP context built by snap_init,
without the force/energy pass. Output: a flat per-particle buffer and its stride,

  bispectrum[ ncoeff * ( cell_particle_offset[cell] + particle ) + component ]

compute_derivative: true (mono-element and chem_flag multi-element configurations) also computes:

- the per-neighbor-pair Jacobian dB_k/dr_j, in a CSR layout:
    bispectrum_deriv_offset[p] .. bispectrum_deriv_offset[p+1]-1  -- rows of particle p
    row 0 of that range                                           -- self term
    rows 1..ninside                                               -- one row per neighbor
    bispectrum_deriv[ (bispectrum_deriv_offset[p]+row)*ncoeff*3 + k*3 + xyz ]
    bispectrum_deriv_nbh_id[ bispectrum_deriv_offset[p]+row ]      -- particle id of that row

- the per-atom aggregate: for atom m, minus the sum of dB_i/dr_m over every atom i that has m as a
  neighbor (or i = m), i.e. force-signed (F = +coeff . aggregate). It is stored as
  ncoeff*3*ntypes grid fields '<prefix>0'..'<prefix>{ncoeff*3*ntypes-1}' (default prefix "sda_"),
  one ncoeff*3 slot per CENTRAL atom type, read by compute_descriptor_snap_global. To use the
  per-atom aggregate itself on more than one MPI rank, reduce the ghost contributions with

    update_opt_from_ghost: { opt_fields: [ "sda_.*" ] }

  after compute_descriptor_snap (and after compute_descriptor_snap_global, which needs the
  un-reduced fields).

The coefficient values of the coefficient file are not used: ncoeff only depends on twojmax and
the number of elements, so the file can be reduced to its per-material header lines, e.g.

  1 0
  Ta 0.5 1

By default the neighbor cutoff is the fixed rcutfac of the parameter file. nneigh_bispectrum > 0
switches to a constant-neighbor-count mode:

  - closest_bispectrum: true (default) -- exact: keeps the nneigh_bispectrum nearest neighbors
    within rcutfac, the cutoff is set to the farthest of them plus neigh_margin.
  - closest_bispectrum: false -- density-scaled estimate: cutoff = rcutfac*(nneigh_bispectrum/n)^(1/3)
    for n neighbors within rcutfac, capped at rcutfac. Cheaper, approximate.

rcutfac must contain at least nneigh_bispectrum neighbors everywhere; this mode does not support
switchinnerflag.

Usage example:

init_parameters:
  - species
  - snap_init:
      parameters: { param: "W.snapparam", coef: "W.snapcoeff" }

compute_descriptor_snap:
  nneigh_bispectrum: 48       # optional, constant-neighbor-count mode
  closest_bispectrum: false   # optional, density-scaled estimate instead of exact sort
  compute_derivative: false   # optional, descriptor derivatives

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_descriptor_snap)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_descriptor_snap", make_grid_variant_operator< ComputeDescriptorSnap > );
  }

}
