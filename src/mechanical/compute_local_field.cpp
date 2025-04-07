#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>


#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>
#include <exaStamp/mechanical/average_local_field.h>
#include <exaStamp/mechanical/compute_local_field.h>

#include <iostream>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;
  
  template< class GridT,
	    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
	    >
  class ComputeLocalFIeldOperator : public OperatorNode
  {    

    // ========= I/O slots =======================
    ADD_SLOT( double                , rcut_max            , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors                   , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                  , ghost               , INPUT , false );
    ADD_SLOT( GridT                 , grid                , INPUT );
    ADD_SLOT( Domain                , domain              , INPUT , REQUIRED );

    ADD_SLOT( double                , rcut                , INPUT , 0.0);
    
    using ComputeBuffer = ComputePairBuffer2<false,false>;
    using ComputeFields = FieldSet< > ;
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

      lout << "\t- Computing per-atom local field" << std::endl;

       // using ForceCPBuf = SimpleNbhComputeBuffer< FieldSet<field::_local_field> >; /* we want extra neighbor storage space to store these fields */
       // ComputePairBufferFactory< ForceCPBuf > local_buf;  
      
      auto local_field = grid->field_accessor( field::local_field );
      auto local_op_fields = make_field_tuple_from_field_set( FieldSet<>{}, local_field );
        
      double rrPotential = *(this->rcut_max);
      LocalFieldOp local_op { *rcut};
      ComputePairNullWeightIterator cp_weight{};
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto local_op_buf = make_compute_pair_buffer<ComputeBuffer>();
      ComputePairTrivialCellFiltering cpu_cell_filter = {}; 	  
     
      if( domain->xform_is_identity() )
        {
          NullXForm cp_xform;
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter );
          compute_cell_particle_pairs( *grid, rrPotential, *ghost, optional, local_op_buf, local_op , local_op_fields, parallel_execution_context() );	  	  
        }
      else
        {
          LinearXForm cp_xform { domain->xform() };
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter  );
          compute_cell_particle_pairs( *grid, rrPotential, *ghost, optional, local_op_buf, local_op , local_op_fields, parallel_execution_context() );	  	  
        }
      lout << "\t- Computing per-atom local field END" << std::endl;

    }
    
  };
  
  template<class GridT> using ComputeLocalFIeldOperatorTmpl = ComputeLocalFIeldOperator<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_local_field) {
    OperatorNodeFactory::instance()->register_factory("compute_local_field",
                                                      make_grid_variant_operator<ComputeLocalFIeldOperatorTmpl>);
  }
}
