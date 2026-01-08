/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

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

#include <exanb/compute/compute_cell_particles.h>
#include <exanb/compute/compute_pair_traits.h>
#include <exaStamp/potential/coulombic/wolf.h>
#include <exaStamp/unit_system.h>

#include <exanb/core/config.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/core/concurent_add_contributions.h>

namespace exaStamp
{
  using namespace exanb;

  struct WolfSelfForceOp
  {
    const ParticleSpecie * m_species = nullptr;
    const WolfParameters m_params = {};
    double m_defaultCharge = 0.0;
    
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () (double& ep, unsigned int type) const
    {
      const double c = m_species[type].m_charge;
      double e_self = -(m_params.e_shift / 2.0 + m_params.alpha / sqrt(M_PI) ) * c * c * m_params.qqrd2e;
      double _ep=0.;
      _ep = EXASTAMP_QUANTITY( e_self * eV );
      ep += _ep;
    }
    inline void operator () (double& ep) const
    {
      ep += 0.;
    }
  };

}

namespace exanb
{
  template<> struct ComputeCellParticlesTraits<exaStamp::WolfSelfForceOp>
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{

  using namespace exanb;
  using namespace onika;
  
  template<class GridT, class _Ep, class _Type>
  class WolfSelf : public OperatorNode
  {      
    ADD_SLOT( GridT            , grid              , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies  , species           , INPUT , REQUIRED );
    ADD_SLOT( WolfParameters    , parameters        , INPUT , REQUIRED );

  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute () override final
    {
      WolfSelfForceOp func = { species->data() , *parameters , species->at(0).m_charge };
      auto ep = grid->field_accessor( onika::soatl::FieldId<_Ep>{} );
      
      onika::soatl::FieldId<_Type> type_field = {};
      if( grid->has_allocated_field( type_field ) )
      {
        auto type = grid->field_accessor( type_field );          
        auto cp_fields = onika::make_flat_tuple( ep, type );
        compute_cell_particles( *grid , false , func , cp_fields , parallel_execution_context() );
      }
      else
      {
        auto cp_fields = onika::make_flat_tuple( ep );
        compute_cell_particles( *grid , false , func , cp_fields , parallel_execution_context() );
      }
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Self contribution for Wolf potential.
)EOF";
    }
    
  };

  template<class GridT> using WolfSelfTmpl = WolfSelf<GridT, field::_ep, field::_type >;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_wolf_self)
  {  
    OperatorNodeFactory::instance()->register_factory( "coulombic_wolf_self" , make_grid_variant_operator<WolfSelfTmpl> );
  }

}
