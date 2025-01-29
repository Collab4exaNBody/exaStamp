//#/* */pragma xstamp



#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <exanb/core/cpp_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/file_utils.h>

#include <exaStamp/potential/snap/snap_params.h>
#include <exaStamp/potential/snap/snap_read_lammps.h>
#include <exaStamp/potential/snap/snap_config.h>
#include <exaStamp/potential/snap/snap_check_bispectrum.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <vector>
#include <memory>
#include <iostream>

#include <mpi.h>

#include "sna.h"
#include "memory.h"

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
  class SnapLMPForce : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( MPI_Comm              , mpi               , INPUT , REQUIRED);
    ADD_SLOT( SnapParms             , parameters        , INPUT , REQUIRED );
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT , false );
    ADD_SLOT( bool                  , conv_coef_units   , INPUT , false );
    ADD_SLOT( bool                  , trigger_thermo_state, INPUT , OPTIONAL );
    ADD_SLOT( GridT                 , grid              , INPUT_OUTPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
    ADD_SLOT( GridParticleLocks     , particle_locks    , INPUT , OPTIONAL , DocString{"particle spin locks"} );

    ADD_SLOT( long                  , timestep          , INPUT , REQUIRED , DocString{"Iteration number"} );
    ADD_SLOT( std::string           , bispectrumchkfile , INPUT , OPTIONAL , DocString{"file with reference values to check bispectrum correctness"} );

    ADD_SLOT( SnapLMPContext        , snap_ctx          , PRIVATE );

    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs
    static constexpr bool UseWeights = false;
    static constexpr bool UseNeighbors = true;
    //static constexpr bool UseLocks = true;
    //    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors>;
    using ComputeBuffer = ComputePairBuffer2<UseWeights,UseNeighbors,SnapComputeBuffer,CopyParticleType>;

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
        std::string lammps_param = onika::data_file_path( parameters->lammps_param );
        std::string lammps_coef = onika::data_file_path( parameters->lammps_coef ); 
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
        //snap_ctx->m_cg_nt = parameters->nt;
        
        for( const auto& mat : snap_ctx->m_config.materials() )
        {
          ldbg << '\t' << mat.name() << ": radelem="<<mat.radelem()<<", weight="<<mat.weight()<<", ncoefs="<<mat.number_of_coefficients()<<std::endl;
          for(size_t i=0;i<mat.number_of_coefficients();i++)
          {
            ldbg << "\t\t" << mat.coefficient(i) << std::endl;
          }
        }
        
        //double jmax = snap_ctx->m_config.twojmax()*0.5;
        int nmat = snap_ctx->m_config.materials().size();
	
        // if( nmat != 1 )
        // {
        //   lerr << "Snap: ERROR: only 1 material is allowed" << std::endl;
	//   //          std::abort();
        // }
      
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

        // snap_ctx->m_factor.assign( MAX_PARTICLE_SPECIES, 1.0 );
        // snap_ctx->m_factor[0] = mat.weight();
        // snap_ctx->m_radelem.assign( MAX_PARTICLE_SPECIES, 0.0 );
        // snap_ctx->m_radelem[0] = mat.radelem();
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
      
      if( snap_ctx->ptr == nullptr )
      {
        snap_ctx->ptr = new LAMMPS_NS::LAMMPS;
        snap_ctx->ptr->error = new LAMMPS_NS::ErrorLogWrapper;
        snap_ctx->ptr->comm = new LAMMPS_NS::CommunicatorInfo;
        snap_ctx->ptr->memory = new LAMMPS_NS::Memory(snap_ctx->ptr);
      }
          
      size_t nt = omp_get_max_threads();
      if( nt > snap_ctx->m_thread_ctx.size() )
      {
        size_t old_nt = snap_ctx->m_thread_ctx.size();
        snap_ctx->m_thread_ctx.resize( nt );
        for(size_t i=old_nt;i<nt;i++)
        {
          assert( snap_ctx->m_thread_ctx[i].sna == nullptr );
          snap_ctx->m_thread_ctx[i].sna =
            new LAMMPS_NS::SNA( snap_ctx->ptr
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
          snap_ctx->m_thread_ctx[i].sna->init();
          snap_ctx->m_thread_ctx[i].sna->grow_rij(1024);
        }
      }

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
      const bool eflag = log_energy || bispectrumchkfile.has_value();
      const bool quadraticflag = snap_ctx->m_config.quadraticflag();
      const bool switchinnerflag = snap_ctx->m_config.switchinnerflag();
      const bool chemflag = snap_ctx->m_config.chemflag();

      // exanb objects to perform computations with neighbors      
      ComputePairNullWeightIterator cp_weight{};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();      
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

      ldbg << "snaplmp: quadratic="<<quadraticflag<<", eflag="<<eflag<<", ncoeff="<<ncoeff<<", ncoeffall="<<ncoeffall<<std::endl;

      if (quadraticflag || eflag)
      {
        // ************ compute_bispectrum(); ****************
        snap_ctx->m_bispectrum.clear();
        snap_ctx->m_bispectrum.resize( total_particles * ncoeff );

        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, ComputePairOptionalLocks<false>{} );
        BispectrumOp bispectrum_op {
                           snap_ctx->m_thread_ctx.data(), snap_ctx->m_thread_ctx.size(),
                           grid->cell_particle_offset_data(), snap_ctx->m_beta.data(), snap_ctx->m_bispectrum.data(),
                           snap_ctx->m_coefs.data(), ncoeff,
                           snap_ctx->m_factor.data(), snap_ctx->m_radelem.data(),
                           nullptr, nullptr,
                           snap_ctx->m_rcut, cutsq,
                           eflag, quadraticflag,
                           switchinnerflag, chemflag
                           };
        compute_cell_particle_pairs( *grid, snap_ctx->m_rcut, *ghost, optional, force_buf, bispectrum_op , compute_bispectrum_field_set , parallel_execution_context() );
        // *********************************************
        if( bispectrumchkfile.has_value() )
        {
          std::ostringstream oss; oss << *bispectrumchkfile << "." << *timestep;
          std::string file_name = onika::data_file_path( oss.str() );
          ldbg << "bispectrumchkfile is set, check bispectrum from file "<< file_name << std::endl;
          snap_check_bispectrum(*mpi, *grid, file_name, ncoeff, snap_ctx->m_bispectrum.data() );
        }
      }

      // // ************ compute_beta(); ****************
      // {        
      //   snap_ctx->m_beta.clear();
      //   snap_ctx->m_beta.resize( total_particles * ncoeff );
      //   for(size_t ii=0;ii<total_particles;ii++)
      //   {
      //     for(int icoeff=0;icoeff<ncoeff;icoeff++)
      //     {
      // 	    snap_ctx->m_beta[ ii * ncoeff + icoeff ] = snap_ctx->m_coefs[icoeff+1];

      // 	    // Here we need to know the particle type to get the proper coefficients (ex: ntypes 0 or 1
      // 	    // We might need to shift the compute_beta into the snap force operator since it is thread dependent.
      // 	    // const int iitype = type[ii];
      // 	    // snap_ctx->m_beta[ ii * ncoeff + icoeff ] = snap_ctx->m_coefs[typeii * (ncoeff+1) + icoeff+1];

      //     }
      //     if (quadraticflag)
      //     {
      //       const double * const coeffi = snap_ctx->m_coefs.data();

      //       int k = ncoeff+1;
      //       for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
      //         double bveci = snap_ctx->m_bispectrum[ ii * ncoeff + icoeff ]; // bispectrum[ii][icoeff];
      //         snap_ctx->m_beta[ ii * ncoeff + icoeff ] /*beta[ii][icoeff]*/ += coeffi[k]*bveci;
      //         k++;
      //         for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
      //           double bvecj = snap_ctx->m_bispectrum[ ii * ncoeff + jcoeff ]; //bispectrum[ii][jcoeff];
      //           snap_ctx->m_beta[ ii * ncoeff + icoeff ] /*beta[ii][icoeff]*/ += coeffi[k]*bvecj;
      //           snap_ctx->m_beta[ ii * ncoeff + jcoeff ] /*beta[ii][jcoeff]*/ += coeffi[k]*bveci;
      //           k++;
      //         }
      //       }
      //     }
      //     //printf("SNAPDBG: beta[%d]:",int(ii));
      //     //for(int icoeff=0;icoeff<ncoeff;icoeff++) printf(" % .3e",snap_ctx->m_beta[ ii * ncoeff + icoeff ]);
      //     //printf("\n");
      //   }
      // }
      // // *********************************************

      auto compute_opt_locks = [&](auto cp_locks)
      {
        auto optional = make_compute_pair_optional_args( nbh_it, cp_weight , cp_xform, cp_locks );
        ForceOp force_op { snap_ctx->m_thread_ctx.data(), snap_ctx->m_thread_ctx.size(),
                           grid->cell_particle_offset_data(), snap_ctx->m_beta.data(), snap_ctx->m_bispectrum.data(),
                           snap_ctx->m_coefs.data(), static_cast<long>(snap_ctx->m_coefs.size()), ncoeff,
                           snap_ctx->m_factor.data(), snap_ctx->m_radelem.data(),
                           nullptr, nullptr,
                           snap_ctx->m_rcut, cutsq,
                           eflag, quadraticflag,
                           switchinnerflag, chemflag,
                           ! (*conv_coef_units) // if coefficients were not converted, then output energy/force must be converted
                           };      
        compute_cell_particle_pairs( *grid, snap_ctx->m_rcut, *ghost, optional, force_buf, force_op , compute_force_field_set , parallel_execution_context() );
      };
      ldbg << "snaplmp: nthreads="<< omp_get_max_threads() <<std::endl;
      
      if( omp_get_max_threads() > 1 ) compute_opt_locks( ComputePairOptionalLocks<true>{ particle_locks->data() } );
      else                            compute_opt_locks( ComputePairOptionalLocks<false>{} );

      ldbg << "snaplmp: done"<<std::endl;
    }

  };

  template<class GridT> using SnapLMPForceTmpl = SnapLMPForce<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(snaplmp)
  {
    OperatorNodeFactory::instance()->register_factory( "snaplmp_force" ,make_grid_variant_operator< SnapLMPForceTmpl > );
  }

}


