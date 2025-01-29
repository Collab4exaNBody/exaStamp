//  // DO NOT REMOVE THIS LINE !!

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
#include <onika/math/basic_types_stream.h>
#include <exanb/core/grid.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/check_particles_inside_cell.h>

#include <onika/soatl/field_tuple.h>

#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include <tuple>

#include <mpi.h>
#include <exanb/mpi/update_ghost_utils.h>
#include <exanb/mpi/ghosts_comm_scheme.h>
#include <exanb/mpi/update_ghosts.h>
#include <onika/mpi/data_types.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace UpdateGhostsUtils;

  // === register factory ===
  template<typename GridT> using UpdateGhostsAllFieldsNoFV = UpdateGhostsNode< GridT , RemoveFields< typename GridT::Fields , FieldSet<field::_fx,field::_fy,field::_fz,field::_ep, field::_vx, field::_vy, field::_vz > > , true >;
  template<typename GridT> using UpdateGhostsRandVandVir = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz,  field::_vx, field::_vy, field::_vz, field::_virial > , false >;
  template<typename GridT> using UpdateGhostsRandV = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz,  field::_vx, field::_vy, field::_vz > , false >;

  template<typename GridT> using UpdateGhostsRandRf = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz, field::_rxf, field::_ryf, field::_rzf > , false >;  
  template<typename GridT> using UpdateGhostsRandRfandV = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz, field::_rxf, field::_ryf, field::_rzf, field::_vx, field::_vy, field::_vz > , false >;
  
  template<typename GridT> using UpdateGhostsRQ = UpdateGhostsNode< GridT , FieldSet<field::_rx, field::_ry, field::_rz , field::_orient > , false >;
  template<typename GridT> using UpdateGhostsIdMol = UpdateGhostsNode< GridT , FieldSet<field::_idmol> , false >;

  ONIKA_AUTORUN_INIT(update_ghosts)
  {
    OperatorNodeFactory::instance()->register_factory( "ghost_update_r_v_vir",   make_grid_variant_operator<UpdateGhostsRandVandVir> );
    OperatorNodeFactory::instance()->register_factory( "ghost_update_r_v",       make_grid_variant_operator<UpdateGhostsRandV> );
    OperatorNodeFactory::instance()->register_factory( "ghost_update_r_rf",      make_grid_variant_operator<UpdateGhostsRandRf> );    
    OperatorNodeFactory::instance()->register_factory( "ghost_update_r_rf_v",    make_grid_variant_operator<UpdateGhostsRandRfandV> );    
    OperatorNodeFactory::instance()->register_factory( "ghost_update_rq",        make_grid_variant_operator<UpdateGhostsRQ> );
    OperatorNodeFactory::instance()->register_factory( "ghost_update_idmol",     make_grid_variant_operator<UpdateGhostsIdMol> );
  }

}

