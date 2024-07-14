// #pragma xstamp_cuda_enable    // DO NOT NOT REMOVE THIS LINE
// #pragma xstamp_grid_variant   // DO NOT NOT REMOVE THIS LINE

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/log.h>
#include <exanb/core/cpp_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/file_utils.h>

#include <exaStamp/potential/snap/snap_params.h>
#include <exaStamp/potential/snap/snap_read_lammps.h>
#include <exaStamp/potential/snap/snap_config.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <vector>
#include <memory>
#include <iostream>

//#include "sna.h"
//#include "memory.h"

#include "snap_context.h"
#include "snap_force_op.h"
#include "snap_bispectrum_op.h"

namespace exaStamp
{

  using namespace exanb;
  using onika::memory::DEFAULT_ALIGNMENT;
//  using namespace SnapExt;


  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class SnapNewForce : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( SnapParms             , parameters        , INPUT , REQUIRED );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT , false );
    ADD_SLOT( bool                  , conv_coef_units   , INPUT , false );
    ADD_SLOT( bool                  , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
    ADD_SLOT( GridParticleLocks     , particle_locks    , INPUT , OPTIONAL , DocString{"particle spin locks"} );

    ADD_SLOT( SnapLMPContext        , snap_ctx          , PRIVATE );

    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    static constexpr bool UseWeights = false;
    static constexpr bool UseNeighbors = true;

//    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors,SnapComputeBuffer,CopyParticleType>; // original
    using ComputeBuffer = ComputePairBuffer2< UseWeights, UseNeighbors
                                            , NoExtraStorage, DefaultComputePairBufferAppendFunc
                                            , exanb::MAX_PARTICLE_NEIGHBORS, ComputePairBuffer2Weights
                                            , FieldSet<field::_type> >;

    using CellParticles = typename GridT::CellParticles;

    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    // using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz >;
    // using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_virial>;
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_type ,field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};
    static constexpr FieldSet< field::_type> compute_bispectrum_field_set{};
        
  public:
    
    // Operator execution
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      //ldbg << "rcut="<<snap_ctx->m_rcut <<std::endl << std::flush;
      if( snap_ctx->m_rcut == 0.0 )
      {
        std::string lammps_param = data_file_path( parameters->lammps_param );
        std::string lammps_coef = data_file_path( parameters->lammps_coef ); 
        ldbg << "Snap: read lammps files "<<lammps_param<<" and "<<lammps_coef<<std::endl << std::flush;
        SnapExt::snap_read_lammps(lammps_param, lammps_coef, snap_ctx->m_config, *conv_coef_units );
        ldbg <<"rfac0="<<snap_ctx->m_config.rfac0() <<", rmin0="<<snap_ctx->m_config.rmin0() <<", rcutfac="<<snap_ctx->m_config.rcutfac() 
             <<", twojmax="<<snap_ctx->m_config.twojmax()<<", bzeroflag="<<snap_ctx->m_config.bzeroflag()<<", nmat="<<snap_ctx->m_config.materials().size() <<std::endl;
        snap_ctx->m_rcut = snap_ctx->m_config.rcutfac(); // because LAMMPS uses angstrom while exastamp uses nm
      }

      *rcut_max = std::max( *rcut_max , snap_ctx->m_rcut );
      
      size_t n_cells = grid->number_of_cells();
      if( n_cells==0 )
      {
        return ;
      }

      if( ! particle_locks.has_value() )
      {
        fatal_error() << "No particle locks" << std::endl;
      }

      if( snap_ctx->m_coefs.empty() )
      {
        /*
        for( const auto& mat : snap_ctx->m_config.materials() )
        {
          ldbg << '\t' << mat.name() << ": radelem="<<mat.radelem()<<", weight="<<mat.weight()<<", ncoefs="<<mat.number_of_coefficients()<<std::endl;
          for(size_t i=0;i<mat.number_of_coefficients();i++)
          {
            ldbg << "\t\t" << mat.coefficient(i) << std::endl;
          }
        }
        */
        
        int nmat = snap_ctx->m_config.materials().size();
	      
        // temporay, enable mutiple species if they all have weight=1. modifications needed for true multimaterial
        snap_ctx->m_factor.assign( nmat, 1.0 );
        snap_ctx->m_radelem.assign( nmat, 0.0 );

        int cnt=0;
        for ( const auto& mat : snap_ctx->m_config.materials() )
        {
          snap_ctx->m_factor[cnt] = mat.weight();	    
          snap_ctx->m_radelem[cnt] = mat.radelem();
          cnt+=1;
        }

        size_t ncoefs_per_specy = snap_ctx->m_config.materials()[0].number_of_coefficients();
        snap_ctx->m_coefs.resize( nmat * ncoefs_per_specy );
        for(int j=0;j<nmat;j++)
        {
          const auto& mat = snap_ctx->m_config.materials()[j];
          for(size_t i=0;i<ncoefs_per_specy;i++)
          {
            snap_ctx->m_coefs[ j * ncoefs_per_specy + i ] = mat.coefficient(i);
          }
        }

      }

      if( snap_ctx->sna == nullptr )
      {
        snap_ctx->sna = new LAMMPS_NS::SNA( new LAMMPS_NS::Memory()
                                          , snap_ctx->m_config.rfac0() 
                                          , snap_ctx->m_config.twojmax() 
                                          , snap_ctx->m_config.rmin0()
                                          , snap_ctx->m_config.switchflag()
                                          , snap_ctx->m_config.bzeroflag()
                                          , snap_ctx->m_config.chemflag()
                                          , snap_ctx->m_config.bnormflag()
                                          , snap_ctx->m_config.wselfallflag()
                                          , snap_ctx->m_config.nelements()
                                          , snap_ctx->m_config.switchinnerflag()
                                          );
        snap_ctx->sna->init();
      }
      ldbg << "Constant config allocation" << std::endl;
      snap_ctx->sna->memory->print( ldbg );

      size_t gpu_block_size = 1;
      size_t gpu_blocks = 1;
      size_t nt = omp_get_max_threads();
      ldbg << "cuda_ctx = " << global_cuda_ctx() << std::endl;
      if( global_cuda_ctx() != nullptr ) ldbg << "cuda has devices = " << global_cuda_ctx()->has_devices() << std::endl;
      if( global_cuda_ctx() != nullptr && global_cuda_ctx()->has_devices() )
      {
        gpu_block_size = onika::parallel::ParallelExecutionContext::gpu_block_size();
        gpu_blocks = global_cuda_ctx()->m_devices[0].m_deviceProp.multiProcessorCount
                                * onika::parallel::ParallelExecutionContext::gpu_sm_mult()
                                + onika::parallel::ParallelExecutionContext::gpu_sm_add();
      }
      if( gpu_blocks > nt ) nt = gpu_blocks;
      ldbg << "Execution buffer scratch: thread blocks = "<<nt<<" , block size = "<<gpu_block_size<<std::endl;
      
      if( nt > snap_ctx->m_scratch.size() )
      {
        size_t cursize = snap_ctx->m_scratch.size();
        ldbg << "grow number of scratch buffers : "<<cursize<<" -> "<<nt<<std::endl;
        snap_ctx->m_scratch.resize( nt );
        for(;cursize<nt;cursize++) snap_ctx->m_scratch[cursize].initialize( snap_ctx->sna , gpu_block_size );
      }

      ldbg << "Max number of neighbors = "<< chunk_neighbors->m_max_neighbors << std::endl;
      for(size_t i=0;i<nt;i++)
      {
        snap_ctx->m_scratch[i].grow_rij( chunk_neighbors->m_max_neighbors );
      }
      ldbg << "Scratch buffers allocation" << std::endl;
      if( nt >= 1 ) snap_ctx->m_scratch[0].memory->print( ldbg );

      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
      {
        ldbg << "trigger_thermo_state = " << *trigger_thermo_state << std::endl;
        log_energy = *trigger_thermo_state ;
      }
      else
      {
        ldbg << "trigger_thermo_state missing " << std::endl;
      }

      const double cutsq = snap_ctx->m_rcut * snap_ctx->m_rcut;
      const bool eflag = log_energy;
      const bool quadraticflag = snap_ctx->m_config.quadraticflag();

      assert( snap_ctx->m_config.switchinnerflag() == snap_ctx->sna->switch_inner_flag );
      assert( snap_ctx->m_config.chemflag() == snap_ctx->sna->chem_flag );

      // exanb objects to perform computations with neighbors      
      ComputePairNullWeightIterator cp_weight{};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>(); // make_compute_pair_buffer< SimpleNbhComputeBuffer< FieldSet<field::_type> > >();
      LinearXForm cp_xform { domain->xform() };

      // constants to resize bispectrum and beta intermediate terms
      const size_t total_particles = grid->number_of_particles();
      size_t ncoefs_per_specy = snap_ctx->m_config.materials()[0].number_of_coefficients();
      int ncoeffall = ncoefs_per_specy; //_per_specysnap_ctx->m_coefs.size() ;
      int ncoeff = -1;
      
      if (!quadraticflag)
        ncoeff = ncoeffall - 1;
      else {
        ncoeff = sqrt(2*ncoeffall)-1;
        int ncoeffq = (ncoeff*(ncoeff+1))/2;
        int ntmp = 1+ncoeff+ncoeffq;
        if (ntmp != ncoeffall) {
          lerr << "Incorrect SNAP coeff file" << std::endl;
          std::abort();
        }
      }

      ldbg << "snap: quadratic="<<quadraticflag<<", eflag="<<eflag<<", ncoeff="<<ncoeff<<", ncoeffall="<<ncoeffall<<std::endl;

      if (quadraticflag || eflag)
      {
        // ************ compute_bispectrum(); ****************
        snap_ctx->m_bispectrum.clear();
        snap_ctx->m_bispectrum.resize( total_particles * ncoeff );

        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, ComputePairOptionalLocks<false>{} );
        BispectrumOp bispectrum_op {
                           LAMMPS_NS::ReadOnlySnapParameters(snap_ctx->sna),
                           snap_ctx->m_scratch.data(), snap_ctx->m_scratch.size(),
                           grid->cell_particle_offset_data(), snap_ctx->m_beta.data(), snap_ctx->m_bispectrum.data(),
                           snap_ctx->m_coefs.data(), ncoeff,
                           snap_ctx->m_factor.data(), snap_ctx->m_radelem.data(),
                           nullptr, nullptr,
                           snap_ctx->m_rcut, cutsq,
                           eflag, quadraticflag };

        auto cp_fields = grid->field_accessors_from_field_set( compute_bispectrum_field_set );
        onika::parallel::BlockParallelForOptions bpfor_opts = {};
        bpfor_opts.max_block_size = gpu_block_size;
        bpfor_opts.fixed_gpu_grid_size = true;
        compute_cell_particle_pairs2( *grid, snap_ctx->m_rcut, *ghost, optional, force_buf, bispectrum_op , cp_fields
                                    , DefaultPositionFields{}, parallel_execution_context(), std::false_type{}, bpfor_opts );

        // *********************************************
      }

      auto compute_opt_locks = [&](auto cp_locks)
      {
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks
                      , ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{}, grid->field_accessors_from_field_set(FieldSet<field::_type>{}) );
        SnapLMPForceOp force_op { LAMMPS_NS::ReadOnlySnapParameters(snap_ctx->sna), snap_ctx->m_scratch.data(), snap_ctx->m_scratch.size(),
                           grid->cell_particle_offset_data(), snap_ctx->m_beta.data(), snap_ctx->m_bispectrum.data(),
                           snap_ctx->m_coefs.data(), static_cast<unsigned int>(snap_ctx->m_coefs.size()), static_cast<unsigned int>(ncoeff),
                           snap_ctx->m_factor.data(), snap_ctx->m_radelem.data(),
                           nullptr, nullptr,
                           snap_ctx->m_rcut, cutsq,
                           eflag, quadraticflag,
                           ! (*conv_coef_units) // if coefficients were not converted, then output energy/force must be converted
                           };
        auto cp_fields = grid->field_accessors_from_field_set( compute_force_field_set );
        onika::parallel::BlockParallelForOptions bpfor_opts = {};
        bpfor_opts.max_block_size = gpu_block_size;
        bpfor_opts.fixed_gpu_grid_size = true;
        compute_cell_particle_pairs2( *grid, snap_ctx->m_rcut, *ghost, optional, force_buf, force_op , cp_fields
                                    , DefaultPositionFields{}, parallel_execution_context(), std::false_type{}, bpfor_opts );
      };
      ldbg << "snap: nthreads="<< omp_get_max_threads() <<std::endl;
      
      bool need_locks = ( omp_get_max_threads() > 1 ) && ( gpu_block_size == 1 );
      if( need_locks ) compute_opt_locks( ComputePairOptionalLocks<true>{ particle_locks->data() } );
      else             compute_opt_locks( ComputePairOptionalLocks<false>{} );

      ldbg << "snap: done"<<std::endl;
    }

  };

  template<class GridT> using SnapNewForceTmpl = SnapNewForce<GridT>;

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "snap_force" ,make_grid_variant_operator< SnapNewForceTmpl > );
  }

}


