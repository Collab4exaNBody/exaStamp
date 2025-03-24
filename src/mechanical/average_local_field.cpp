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
#include <exaStamp/mechanical/average_local_field.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include <exaStamp/compute/thermodynamic_state.h>

#include <exanb/mpi/grid_update_ghosts.h>

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
  class AverageLocalFieldOperator : public OperatorNode
  {    

    using UpdateGhostsScratch = typename UpdateGhostsUtils::UpdateGhostsScratch;
    
    // ========= I/O slots =======================
    ADD_SLOT( double                , rcut_max            , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors                   , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    //    ADD_SLOT( bool                  , ghost               , INPUT , false );
    ADD_SLOT( GridT                 , grid                , INPUT );
    ADD_SLOT( Domain                , domain              , INPUT , REQUIRED );
    ADD_SLOT( double                , rcut                , INPUT , REQUIRED);
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme , INPUT_OUTPUT , OPTIONAL );
    
    ADD_SLOT( MPI_Comm                 , mpi               , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( long                     , mpi_tag           , INPUT , 0 );
    ADD_SLOT( bool                     , gpu_buffer_pack   , INPUT , false );
    ADD_SLOT( bool                     , async_buffer_pack , INPUT , false );
    ADD_SLOT( bool                     , staging_buffer    , INPUT , false );
    ADD_SLOT( bool                     , serialize_pack_send , INPUT , false );
    ADD_SLOT( bool                     , wait_all          , INPUT , false );

    ADD_SLOT( UpdateGhostsScratch      , ghost_comm_buffers, PRIVATE );
    
    using ComputeBuffer = ComputePairBuffer2<false,false>;
    using ComputeFields = FieldSet< > ;
    static constexpr ComputeFields compute_force_field_set{};
    
  public:
  
    inline void execute () override final
    {

      if( ! ghost_comm_scheme.has_value() ) return;
      if( grid->number_of_particles() == 0 ) return;

      auto pecfunc = [self=this](auto ... args) { return self->parallel_execution_context(args...); };
      auto pesfunc = [self=this](unsigned int i) { return self->parallel_execution_stream(i); };
      
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      bool has_chunk_neighbors = chunk_neighbors.has_value();
      if( !has_chunk_neighbors )
      {
        lerr << "No neighbors input data available" << std::endl;
        std::abort();
      }

      using ForceCPBuf = SimpleNbhComputeBuffer< FieldSet<field::_local_field> >; /* we want extra neighbor storage space to store these fields */
      ComputePairBufferFactory< ForceCPBuf > local_op_buf;  


      auto local_field = grid->field_accessor( field::local_field );
      auto ghost_update_fields = onika::make_flat_tuple( local_field );

      grid_update_ghosts( ldbg, *mpi, *ghost_comm_scheme, *grid, *domain, nullptr,
                          *ghost_comm_buffers, pecfunc, pesfunc, ghost_update_fields,
                          *mpi_tag, *gpu_buffer_pack, *async_buffer_pack, *staging_buffer,
                          *serialize_pack_send, *wait_all, std::false_type{} );
      
      auto local_op_fields = make_field_tuple_from_field_set( FieldSet<>{}, local_field );
        
      double rrPotential = *(this->rcut_max);
      AverageOp local_op { *rcut};
      ComputePairNullWeightIterator cp_weight{};
      ComputePairOptionalLocks<false> cp_locks {};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      //      auto local_op_buf = make_compute_pair_buffer<ComputeBuffer>();
      
      ComputePairTrivialCellFiltering cpu_cell_filter = {}; 	  
      ComputePairTrivialParticleFiltering cpu_particle_filter = {}; 	  

      auto c_local_field = grid->field_const_accessor( field::local_field );
      auto nbh_fields = onika::make_flat_tuple(c_local_field); // we want external field c_emb_field to be populated for each neighbor
      
      if( domain->xform_is_identity() )
        {
          NullXForm cp_xform;
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter, cpu_particle_filter, nbh_fields );
          compute_cell_particle_pairs( *grid, rrPotential, false, optional, local_op_buf, local_op , local_op_fields, parallel_execution_context() );	  	  
        }
      else
        {
          LinearXForm cp_xform { domain->xform() };
          auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter, cpu_particle_filter, nbh_fields  );
          compute_cell_particle_pairs( *grid, rrPotential, false, optional, local_op_buf, local_op , local_op_fields, parallel_execution_context() );	  	  
        }

    }
    
  };
  
  template<class GridT> using AverageLocalFieldOperatorTmpl = AverageLocalFieldOperator<GridT>;
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(average_local_field)
  {
    OperatorNodeFactory::instance()->register_factory( "average_local_field", make_grid_variant_operator< AverageLocalFieldOperatorTmpl > );
  }

}
