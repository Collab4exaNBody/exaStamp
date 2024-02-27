// #pragma xstamp_cuda_enable // DO NOT REMOVE THIS LINE !!

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/field_sets.h>
#include <exanb/core/check_particles_inside_cell.h>

#include <exanb/mpi/data_types.h>
#include <exanb/mpi/ghosts_comm_scheme.h>

#include <mpi.h>
#include <vector>
#include <string>
#include <list>
#include <algorithm>

#include <exanb/mpi/update_from_ghost_utils.h>
#include <exanb/mpi/ghosts_comm_scheme.h>
#include <exanb/mpi/update_from_ghosts.h>
#include <exanb/grid_cell_particles/cell_particle_update_functor.h>

namespace exaStamp
{
  using namespace exanb;
    
  // === register factory ===
  template<typename GridT> using UpdateForceEnergyFromGhosts = UpdateFromGhosts< GridT , FieldSet<field::_fx,field::_fy,field::_fz, field::_ep>, UpdateValueAdd >;
  template<typename GridT> using UpdateFlatForceEnergyFromGhosts = UpdateFromGhosts< GridT , FieldSet<field::_flat_fx,field::_flat_fy,field::_flat_fz, field::_flat_ep>, UpdateValueAdd >;
  template<typename GridT> using UpdateVirialForceEnergyFromGhosts = UpdateFromGhosts< GridT , FieldSet<field::_fx,field::_fy,field::_fz, field::_ep, field::_virial>, UpdateValueAdd >; //NICO
  template<typename GridT> using UpdateFromGhostsTestId = UpdateFromGhosts< GridT , FieldSet<field::_id>, UpdateValueAssertEqual >;

  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "update_force_energy_from_ghost", make_grid_variant_operator<UpdateForceEnergyFromGhosts> );
    OperatorNodeFactory::instance()->register_factory( "flat_force_energy_from_ghost", make_grid_variant_operator<UpdateFlatForceEnergyFromGhosts> );
    OperatorNodeFactory::instance()->register_factory( "update_virial_force_energy_from_ghost", make_grid_variant_operator<UpdateVirialForceEnergyFromGhosts> );
    OperatorNodeFactory::instance()->register_factory( "update_from_ghost_check_id", make_grid_variant_operator<UpdateFromGhostsTestId> );
  }

}

