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
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/file_utils.h>
#include <onika/cuda/cuda_math.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <md/snap/snap_params.h>
#include <md/snap/snap_read_lammps.h>
#include <md/snap/snap_config.h>
#include <md/snap/snap_context.h>
#include <md/snap/snap_compute_buffer.h>
#include <md/snap/snap_bispectrum_op.h>
#include <md/snap/snap_check_bispectrum.h>
#include <md/snap/sna.h>

#include <cmath>
#include <sstream>
#include <mpi.h>

// GPU-compatible bispectrum-only pass, lifted out of md::SnapForceRealT (snap_force.h):
// same LAMMPS param/coef loading + SNA setup + BispectrumOpRealT compute-pair call, but
// with the force/energy pass dropped entirely. BispectrumOpRealT only ever writes to its
// own central particle's slot (see snap_bispectrum_op.h:134, "bispectrum + bispectrum_ii_offset"),
// so unlike snap_force this needs no particle locks. coeffelem/beta/eflag/quadraticflag are
// accepted by BispectrumOpRealT for interface parity with the force op but are never read
// inside its operator() -- only radelem/wjelem (species radius/weight) and snaconf feed the
// actual bispectrum values -- so the LAMMPS coefficient array itself is never populated here.
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
  class ComputeBispectrum : public OperatorNode
  {
    using RealT = double;
    using SnapContext = md::SnapXSContextRealT<RealT>;
    template<int I> using ICST = onika::IntConst<I>;
    template<int jm> using ROParamsMonoElem = SnapInternal::ReadOnlySnapParametersRealT<RealT,ICST<jm>,ICST<1>,false>;

    ADD_SLOT( MPI_Comm                 , mpi               , INPUT , REQUIRED );
    ADD_SLOT( md::SnapParms            , parameters        , INPUT , REQUIRED , DocString{"LAMMPS-format SNAP parameter/coefficient files (param/coef), see snap_force"} );
    ADD_SLOT( double                   , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors  , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                     , ghost             , INPUT , false );
    ADD_SLOT( bool                     , conv_coef_units   , INPUT , false );
    ADD_SLOT( GridT                    , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                   , domain            , INPUT , REQUIRED );
    ADD_SLOT( long                     , timestep          , INPUT , REQUIRED );
    ADD_SLOT( std::string              , bispectrumchkfile , INPUT , OPTIONAL , DocString{"file with reference values to check bispectrum correctness"} );
    ADD_SLOT( double                   , check_bs_max_error, INPUT , 1.e-12 );
    ADD_SLOT( long                     , nneigh_bispectrum , INPUT , 0 , DocString{"If >0, switches to constant-neighbor-count mode: for each particle, the cutoff is adapted so (approximately) this many neighbors contribute, instead of the fixed physical rcut from the SNAP param file. 0 (default) keeps the fixed-rcut behavior. The param file's rcutfac must still be generous enough to gather at least this many neighbors."} );
    ADD_SLOT( double                   , neigh_margin      , INPUT , 0.01 , DocString{"closest_bispectrum=true only: small margin added past the Nth-nearest-neighbor distance, in the same length unit as positions, so that neighbor sits inside rather than exactly on the smooth cutoff's edge."} );
    ADD_SLOT( bool                     , closest_bispectrum, INPUT , true , DocString{"When nneigh_bispectrum>0: true (default) picks the exact Nth-nearest-neighbor distance as cutoff (sort-based, exact but O(n^2) per particle); false uses a density-scaled estimate instead -- assumes neighbor count scales as r^3 with uniform local density, so rcut = min(rcutfac, rcutfac*(nneigh_bispectrum/n_within_rcutfac)^(1/3)) with no sorting needed, cheaper but approximate."} );

    ADD_SLOT( SnapContext              , snap_ctx          , PRIVATE );
    ADD_SLOT( onika::memory::CudaMMVector<RealT> , bispectrum , OUTPUT , DocString{"Flat per-particle bispectrum buffer: bispectrum[ ncoeff*(cell_particle_offset[cell]+particle) + component ], see grid->cell_particle_offset_data()"} );
    ADD_SLOT( long                     , ncoeff            , OUTPUT , DocString{"Number of bispectrum coefficients per particle (stride of the bispectrum buffer)"} );

  public:
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      if( snap_ctx->m_rcut == 0.0 )
      {
        std::string lammps_param = onika::data_file_path( parameters->lammps_param );
        std::string lammps_coef = onika::data_file_path( parameters->lammps_coef );
        ldbg << "compute_bispectrum: read lammps files "<<lammps_param<<" and "<<lammps_coef<<std::endl;
        SnapExt::snap_read_lammps(lammps_param, lammps_coef, snap_ctx->m_config, *conv_coef_units );
        snap_ctx->m_rcut = snap_ctx->m_config.rcutfac();
      }
      *rcut_max = std::max( double(*rcut_max) , double(snap_ctx->m_rcut) );

      if( grid->number_of_cells() == 0 ) { *ncoeff = 0; return; }

      if( snap_ctx->m_factor.empty() )
      {
        int nmat = snap_ctx->m_config.materials().size();
        snap_ctx->m_factor.assign( nmat, 1.0 );
        snap_ctx->m_radelem.assign( nmat, 0.0 );
        int cnt=0;
        for ( const auto& mat : snap_ctx->m_config.materials() )
        {
          snap_ctx->m_factor[cnt] = mat.weight();
          snap_ctx->m_radelem[cnt] = mat.radelem();
          cnt+=1;
        }
      }

      if( snap_ctx->sna == nullptr )
      {
        snap_ctx->sna = new SnapInternal::SNARealT<RealT>( new SnapInternal::Memory()
                                          , snap_ctx->m_config.rfac0(), snap_ctx->m_config.twojmax(), snap_ctx->m_config.rmin0()
                                          , snap_ctx->m_config.switchflag(), snap_ctx->m_config.bzeroflag(), snap_ctx->m_config.chemflag()
                                          , snap_ctx->m_config.bnormflag(), snap_ctx->m_config.wselfallflag(), snap_ctx->m_config.nelements()
                                          , snap_ctx->m_config.switchinnerflag() );
        snap_ctx->sna->init();
      }

      // ncoeff (bispectrum stride) comes straight from twojmax/nelements via the SNA config
      // itself (see sna.h's compute_coeff_count) -- it never depends on the coefficient file,
      // so the coef file's actual coefficient values (and even their count) are irrelevant here;
      // only its per-material header lines (name, radelem, weight) matter for this operator.
      const bool quadraticflag = snap_ctx->m_config.quadraticflag();
      const int ncoeff_local = snap_ctx->sna->ncoeff;
      *ncoeff = ncoeff_local;

      const size_t total_particles = grid->number_of_particles();
      bispectrum->clear();
      bispectrum->resize( total_particles * ncoeff_local );

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
        auto bs_buf = make_compute_pair_buffer< ComputeBufferBS , ResetBispectrumCPBuf >();
        auto cp_fields = grid->field_accessors_from_field_set( compute_bispectrum_field_set );

        if( *nneigh_bispectrum > 0 )
        {
          if( snap_ctx->m_config.switchinnerflag() )
          {
            fatal_error() << "compute_bispectrum: nneigh_bispectrum (constant-neighbor-count) mode does not support switchinnerflag" << std::endl;
          }
          ConstNeighBispectrumOpRealT<RealT,RealT,SnapConfParamsT> bispectrum_op {
                               snapconf, grid->cell_particle_offset_data(), bispectrum->data(),
                               snap_ctx->m_factor.data(), ncoeff_local, static_cast<int>(*nneigh_bispectrum), *neigh_margin,
                               snap_ctx->m_rcut, *closest_bispectrum };
          compute_cell_particle_pairs2( *grid, snap_ctx->m_rcut, *ghost, optional, bs_buf, bispectrum_op, cp_fields
                                       , DefaultPositionFields{}, parallel_execution_context() );
        }
        else
        {
          md::BispectrumOpRealT<RealT,RealT,SnapConfParamsT> bispectrum_op {
                               snapconf, grid->cell_particle_offset_data(), nullptr, bispectrum->data(),
                               nullptr, ncoeff_local, snap_ctx->m_factor.data(), snap_ctx->m_radelem.data(),
                               nullptr, nullptr, snap_ctx->m_rcut, false, quadraticflag };
          compute_cell_particle_pairs2( *grid, snap_ctx->m_rcut, *ghost, optional, bs_buf, bispectrum_op, cp_fields
                                       , DefaultPositionFields{}, parallel_execution_context() );
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

GPU-compatible, bispectrum-only SNAP pass. Reuses the same LAMMPS param/coef loading and
SNA setup as snap_force, but only runs BispectrumOpRealT -- no force/energy pass, no
particle locks needed. Output is a flat per-particle buffer (bispectrum) plus its
per-particle stride (ncoeff), indexed the same way snap_force's internal bispectrum
buffer is (see md/snap/snap_check_bispectrum.h):

  bispectrum[ ncoeff * ( cell_particle_offset[cell] + particle ) + component ]

where cell_particle_offset comes from grid->cell_particle_offset_data().

The coefficient file's actual coefficient values (and their count) are never used --
ncoeff is derived purely from twojmax/nelements (see sna.h's compute_ncoeff), not from
the coefficient file. Only each material's header line (name, radelem, weight) matters,
so the coef file can be reduced to just that, e.g. for a single-element potential:

  1 0
  Ta 0.5 1

(nmat=1, 0 coefficients per material, followed by one "name radelem weight" line per
material and no coefficient lines at all).

By default the neighbor cutoff is the fixed physical rcutfac from the param file. Setting
nneigh_bispectrum>0 switches to a constant-neighbor-count mode instead (mirrors
compute_local_structural_metrics.h's constant-neighbor-count modes), with two variants
selected by closest_bispectrum:

  - closest_bispectrum: true (default) -- exact mode. Sorts all neighbors within rcutfac
    by distance and keeps the nearest nneigh_bispectrum of them; the smooth cutoff radius
    is set to their farthest distance plus neigh_margin. O(n^2) sort per particle, but exact.

  - closest_bispectrum: false -- density-scaled estimate. No sorting: given n neighbors
    found within rcutfac, assumes neighbor count scales as r^3 with uniform local density
    and sets the cutoff to rcutfac * (nneigh_bispectrum/n)^(1/3), capped at rcutfac. Cheaper,
    but only approximate, and can end up with somewhat more or fewer than nneigh_bispectrum
    neighbors in non-uniform/anisotropic local environments.

Either way, rcutfac must be generous enough to contain at least nneigh_bispectrum neighbors
everywhere in the system; this mode does not support switchinnerflag.

Usage example:

compute_bispectrum:
  parameters: { param: "W.snapparam", coef: "W.snapcoeff" }
  nneigh_bispectrum: 48       # optional, constant-neighbor-count mode
  closest_bispectrum: false   # optional, density-scaled estimate instead of exact sort

)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_bispectrum)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_bispectrum", make_grid_variant_operator< ComputeBispectrum > );
  }

}
