// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

#include <exaStamp/potential/eam/eam_buffer.h>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/mpi/grid_update_ghosts.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace UpdateGhostsUtils;

  template<class GridT>
  class EAMUpdateGhostsEMB : public OperatorNode
  {
    using UpdateGhostsScratch = typename UpdateGhostsUtils::UpdateGhostsScratch;

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( MPI_Comm                 , mpi               , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme , INPUT_OUTPUT , OPTIONAL );
    ADD_SLOT( GridT                    , grid              , INPUT_OUTPUT);
    ADD_SLOT( Domain                   , domain            , INPUT );
    ADD_SLOT( long                     , mpi_tag           , INPUT , 0 );

    ADD_SLOT( bool                     , gpu_buffer_pack   , INPUT , false );
    ADD_SLOT( bool                     , async_buffer_pack , INPUT , false );
    ADD_SLOT( bool                     , staging_buffer    , INPUT , false );
    ADD_SLOT( bool                     , serialize_pack_send , INPUT , false );
    ADD_SLOT( bool                     , wait_all          , INPUT , false );

    ADD_SLOT( UpdateGhostsScratch      , ghost_comm_buffers, PRIVATE );

  public:
    inline void execute() override final
    {
      if( ! ghost_comm_scheme.has_value() ) return;
      if( grid->number_of_particles() == 0 ) return;
    
      auto pecfunc = [self=this](auto ... args) { return self->parallel_execution_context(args...); };
      auto pesfunc = [self=this](unsigned int i) { return self->parallel_execution_stream(i); };

      auto rho_emb_field = grid->field_accessor( field::rho_dEmb );
      auto ghost_update_fields = onika::make_flat_tuple( rho_emb_field );

      grid_update_ghosts( ldbg, *mpi, *ghost_comm_scheme, *grid, *domain, nullptr,
                          *ghost_comm_buffers, pecfunc, pesfunc, ghost_update_fields,
                          *mpi_tag, *gpu_buffer_pack, *async_buffer_pack, *staging_buffer,
                          *serialize_pack_send, *wait_all, std::false_type{} );
    }

  };

  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "ghost_update_emb", make_grid_variant_operator<EAMUpdateGhostsEMB> );
  }

}

