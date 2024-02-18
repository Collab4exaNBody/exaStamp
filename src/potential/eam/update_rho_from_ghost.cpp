// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <exanb/mpi/grid_update_from_ghosts.h>
#include <exaStamp/potential/eam/eam_buffer.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace UpdateGhostsUtils;

  template<class GridT>
  class EAMUpdateRhoFromGhosts : public OperatorNode
  {
    using UpdateGhostsScratch = typename UpdateGhostsUtils::UpdateGhostsScratch;

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( MPI_Comm                 , mpi               , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( GhostCommunicationScheme , ghost_comm_scheme , INPUT_OUTPUT , OPTIONAL );
    ADD_SLOT( GridT                    , grid              , INPUT_OUTPUT);
    ADD_SLOT( long                     , mpi_tag           , INPUT , 0 );

    ADD_SLOT( bool                     , gpu_buffer_pack   , INPUT , false );
    ADD_SLOT( bool                     , async_buffer_pack , INPUT , false );
    ADD_SLOT( bool                     , staging_buffer    , INPUT , false );
    ADD_SLOT( bool                     , serialize_pack_send , INPUT , false );
    ADD_SLOT( bool                     , wait_all          , INPUT , false );

    ADD_SLOT( EamAdditionalFields      , eam_extra_fields  , INPUT_OUTPUT );
    ADD_SLOT( UpdateGhostsScratch      , ghost_comm_buffers, PRIVATE );

  public:
    inline void execute() override final
    {
      if( ! ghost_comm_scheme.has_value() ) return;
      if( grid->number_of_particles() == 0 ) return;

      auto pecfunc = [self=this]() { return self->parallel_execution_context(); };
      auto pesfunc = [self=this](unsigned int i) { return self->parallel_execution_stream(i); };

      eam_extra_fields->m_rho.resize( grid->number_of_particles() );
      double * rho_ptr = eam_extra_fields->m_rho.data();
      auto rho_field = make_external_field_flat_array_accessor( *grid , rho_ptr , field::rho );
      auto update_fields = onika::make_flat_tuple( rho_field );

      grid_update_from_ghosts( ldbg, *mpi, *ghost_comm_scheme, *grid, nullptr,
                          *ghost_comm_buffers, pecfunc,pesfunc, update_fields,
                          *mpi_tag, *gpu_buffer_pack, *async_buffer_pack, *staging_buffer,
                          *serialize_pack_send, *wait_all, UpdateValueAdd{} );
    }

  };

  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "update_rho_from_ghost", make_grid_variant_operator<EAMUpdateRhoFromGhosts> );
  }

}

