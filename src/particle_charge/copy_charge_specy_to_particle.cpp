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

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_stream.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <vector>

namespace exaStamp
{

  using namespace exanb;

  template<typename GridT, class = AssertGridHasFields< GridT, field::_type, field::_charge> >
  struct CopyChargeFromSpecyToParticleNode : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies , species , INPUT , REQUIRED );
    ADD_SLOT( GridT           , grid    , INPUT_OUTPUT );

    inline void execute () override final
    {
      const ParticleSpecie* __restrict__ sp = species->data();

      auto cells = grid->cells();
      IJK dims = grid->dimension();
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,i,loc)
        {
          const uint8_t * __restrict__ ptype = cells[i][field::type];
          double* __restrict__ charge = cells[i][field::charge];
          size_t n = cells[i].size();
#         pragma omp simd
          for(size_t j=0;j<n;j++)
          {
            assert( ptype[j]>=0 && ptype[j]<species->size() );
            charge[j] = sp[ ptype[j] ].m_charge;
          }
        }
        GRID_OMP_FOR_END
      }
    }

  };
  
  namespace tmplhelper { template<class GridT> using CopyChargeFromSpecyToParticleNode = ::exaStamp::CopyChargeFromSpecyToParticleNode<GridT>; }

  // === register factories ===  
  ONIKA_AUTORUN_INIT(copy_charge_specy_to_particle)
  {
    OperatorNodeFactory::instance()->register_factory( "copy_charge_specy_to_particle", make_grid_variant_operator< tmplhelper::CopyChargeFromSpecyToParticleNode > );
  }

}

