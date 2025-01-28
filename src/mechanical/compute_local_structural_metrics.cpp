#pragma xstamp_grid_variant

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <exanb/core/cpp_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/file_utils.h>

#include <exaStamp/potential/snaplegacy/SnapLegacyCG.h>
#include <exaStamp/potential/snaplegacy/SnapLegacyBS.h>
#include <exaStamp/potential/snaplegacy/snap_legacy_read_lammps.h>
#include <exaStamp/potential/snaplegacy/snap_legacy_config.h>

#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>
#include <exaStamp/mechanical/compute_local_structural_metrics.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include <vector>
#include <memory>
#include <iostream>

// this allows for parallel compilation of templated operator for each available field set

//

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;
  
  template< class GridT,
	    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
	    >
  class ComputeLocalStructuralMetricsOperator : public OperatorNode
  {    

    using VariablesVec = std::vector<std::string>;
    using DegreesVec = std::vector<int>;
    //DegreesVec default_degree = {5, 4, 6, 8, 10, 12};
    
    // ========= I/O slots =======================
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT , false );
    ADD_SLOT( GridT                 , grid              , INPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
    
    ADD_SLOT( GridParticleLocalStructuralMetrics, local_structural_data , OUTPUT );
    ADD_SLOT( VariablesVec                      , per_atom_data         , INPUT );
    
    ADD_SLOT( double                            , rcut_bispectrum       , INPUT , 0.0);    
    ADD_SLOT( uint64_t                          , nnn_bispectrum        , INPUT , 48);    
    ADD_SLOT( double                            , jmax_bispectrum       , INPUT , 3);
    ADD_SLOT( bool                              , closest_bispectrum    , INPUT , true);
    
    ADD_SLOT( double                            , rcut_steinhardt       , INPUT , 0.0);    
    ADD_SLOT( uint64_t                          , nnn_steinhardt        , INPUT , 12);    
    ADD_SLOT( DegreesVec                        , degree_steinhardt     , INPUT , DegreesVec({5, 4, 6, 8, 10, 12}) );

    // shortcut to the Compute buffer used (and passed to functor) by compute_cell_particle_pairs    
    using ComputeBuffer = ComputePairBuffer2<false,false>;
    //using CellParticles = typename GridT::CellParticles;
    //using ParticleLock = decltype( ComputePairOptionalLocks<false>{}[0][0] );
    // compile time constant indicating if grid has virial field
    static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

    // attributes processed during computation
    //using ComputeFieldsWithoutVirial = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz >;
    //using ComputeFieldsWithVirial    = FieldSet< field::_ep ,field::_fx ,field::_fy ,field::_fz ,field::_virial >;
    using ComputeFields = FieldSet< > ; //std::conditional_t< has_virial_field , ComputeFieldsWithVirial , ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};
    
  public:
  
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      bool has_chunk_neighbors = chunk_neighbors.has_value();
      if( !has_chunk_neighbors )
      {
        lerr << "No neighbors input data available" << std::endl;
        std::abort();
      }

      //      GridT& grid = *(this->grid);      
      size_t n_cells = grid->number_of_cells();
      auto cells = grid->cells();
      
      GridParticleLocalStructuralMetrics& local_structural_data = *(this->local_structural_data);
      local_structural_data.resize( n_cells );
      
      VariablesVec per_atom_data = *(this->per_atom_data);
      std::sort(per_atom_data.begin(), per_atom_data.end());

      // Computation of per-atom bispectrum only if asked     
      if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "bispectrum"))
	{
	  lout << "\t- Gathering per-atom bispectrum" << std::endl;
	  for(size_t c=0; c<n_cells;++c)
	    {
	      local_structural_data[c].bispectrum.resize( cells[c].size() );
	    }

	  // We declare and resize the structure that contains the per-atom bispectrum
	  double rrBispectrum = *(this->rcut_bispectrum);
	  double rrPotential = *(this->rcut_max);
	  int nnBispectrum = *(this->nnn_bispectrum);
	  double jmax = (*this->jmax_bispectrum);
	  bool closest = (*this->closest_bispectrum);
      
	  if( m_cg == nullptr )
	    {
	      m_cg_nt = 2;
	      m_cg = std::make_shared<SnapLegacyCG>( jmax , m_cg_nt );	
	      m_cg->compute();        
	    }
      
	  size_t nt = omp_get_max_threads();
	  if( nt > m_thread_ctx.size() )
	    {
	      size_t old_nt = m_thread_ctx.size();
	      m_thread_ctx.resize( nt );
	      for(size_t i=old_nt;i<nt;i++)
		{
		  assert( m_thread_ctx[i].m_snapbs == nullptr );
		  m_thread_ctx[i].m_snapbs = std::make_shared<SnapLegacyBS>( m_cg->get_jmax(), nullptr , nullptr );
		}
	    }

	  BispectrumOp force_op { *m_cg, m_thread_ctx, rrPotential, rrBispectrum, nnBispectrum, closest, local_structural_data};
	  //BispectrumOpBis<GridT> force_op { *m_cg, m_thread_ctx, rrPotential, rrBispectrum, nnBispectrum, local_structural_data};	  	  

	  ComputePairNullWeightIterator cp_weight{};
	  ComputePairOptionalLocks<false> cp_locks {};
	  exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
	  auto force_buf = make_compute_pair_buffer<ComputeBuffer>();
	  ComputePairTrivialCellFiltering cpu_cell_filter = {}; 	  
	  
	  if( domain->xform_is_identity() )
	    {
	      NullXForm cp_xform;
	      auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter );
	      compute_cell_particle_pairs( *grid, rrPotential, *ghost, optional, force_buf, force_op , compute_force_field_set, parallel_execution_context() );	  	  
	    }
	  else
	    {
	      LinearXForm cp_xform { domain->xform() };
	      auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter  );
	      compute_cell_particle_pairs( *grid, rrPotential, *ghost, optional, force_buf, force_op , compute_force_field_set, parallel_execution_context() );	  	  
	    }

	  lout << "Bispectrum computed" << std::endl;
	  
	}

      // Computation of per-atom bispectrum only if asked TODO
      if (std::binary_search(per_atom_data.begin(), per_atom_data.end(), "steinhardt"))
	{
	  lerr << "Steinhardt parameters not yet implemented (ongoing: Paul Lafourcade)" << std::endl;
	}
      
    }

    //  private:
        
    std::vector<PerThreadContext> m_thread_ctx;
    std::shared_ptr<SnapLegacyCG> m_cg = nullptr;
    double m_rcut = 0.0;
    int m_cg_nt = 2;    
};

  template<class GridT> using ComputeLocalStructuralMetricsOperatorTmpl = ComputeLocalStructuralMetricsOperator<GridT>;
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(compute_local_structural_metrics)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_local_structural_metrics", make_grid_variant_operator< ComputeLocalStructuralMetricsOperatorTmpl > );
  }

}
