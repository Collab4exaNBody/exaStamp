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
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/geometry.h>
#include <onika/math/basic_types_stream.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/log.h>
#include <onika/thread.h>
#include <exanb/grid_cell_particles/replicate_domain.h>

#include <vector>
#include <iostream>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  struct ShiftParticleIdAndCMolFunctor
  {
    template<class FieldTupleT>
    inline void operator () ( FieldTupleT& tp , int64_t id_offset ) const
    {
      if constexpr ( tp.has_field( field::id ) )
      {
        tp[field::id] += id_offset;
      }
      if constexpr ( tp.has_field( field::cmol ) )
      {
        auto cmol = tp[field::cmol];
        for(size_t i=0;i<cmol.size();i++)
        {
          if( cmol[i] != std::numeric_limits<uint64_t>::max() ) cmol[i] += id_offset;
        }
        tp[field::cmol] = cmol;
      }
    }
  };

  template<class GridT> using ReplicateDomainCMol = ReplicateDomain<GridT,ShiftParticleIdAndCMolFunctor>;

   // === register factories ===
  ONIKA_AUTORUN_INIT(replicate_domain_cmol)
  {
    OperatorNodeFactory::instance()->register_factory( "replicate_domain_cmol", make_grid_variant_operator< ReplicateDomainCMol > );
  }

}

