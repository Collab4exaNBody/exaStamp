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
#include <onika/cpp_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/file_utils.h>

#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>
#include <exaStamp/mechanical/compute_local_centrosymmetry.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include <vector>
#include <memory>
#include <iostream>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;
  
  template< class GridT,
	    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
	    >
  class ComputeLocalCentrosymmetryOperator : public OperatorNode
  {    

    // ========= I/O slots =======================
    ADD_SLOT( double                , rcut_max          , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost             , INPUT , false );
    ADD_SLOT( GridT                 , grid              , INPUT );
    ADD_SLOT( Domain                , domain            , INPUT , REQUIRED );
    
    ADD_SLOT( GridParticleLocalStructuralMetrics, local_structural_data , INPUT_OUTPUT );
    ADD_SLOT( double                            , rcut                  , INPUT , REQUIRED);
    ADD_SLOT( uint64_t                          , nnn                   , INPUT , 12);

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

      double rrPotential = *(this->rcut_max);
      
      size_t n_cells = grid->number_of_cells();
      auto cells = grid->cells();
      GridParticleLocalStructuralMetrics& local_structural_data = *(this->local_structural_data);
      local_structural_data.resize( n_cells );
      
      lout << "\t- Computing per-atom CSP" << std::endl;
      for(size_t c=0; c<n_cells;++c)
        {
          local_structural_data[c].csp.resize( cells[c].size() );
        }

      CentroSymmetryOp force_op { *rcut, *nnn, local_structural_data};
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
      
    }
    
  };
  
  template<class GridT> using ComputeLocalCentrosymmetryOperatorTmpl = ComputeLocalCentrosymmetryOperator<GridT>;
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(compute_local_centrosymmetry)
  {
    OperatorNodeFactory::instance()->register_factory( "compute_local_centrosymmetry", make_grid_variant_operator< ComputeLocalCentrosymmetryOperatorTmpl > );
  }

}
