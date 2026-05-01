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
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>
#include <exanb/core/particle_type_id.h>

#include <omp.h>

#include "airebo_params.h"

namespace exaStamp
{

  using namespace exanb;

  class AireboInit : public OperatorNode
  {
    ADD_SLOT( AireboParams       , parameters , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double             , rcut_max   , INPUT_OUTPUT , 0.0 );  
    ADD_SLOT( double             , bondorder_cutoff    , OUTPUT);    
    ADD_SLOT( double             , rebo_cutoff         , OUTPUT );   
    ADD_SLOT( double             , lj_cutoff           , OUTPUT );    
    ADD_SLOT( double             , torsion_cutoff      , OUTPUT );
    ADD_SLOT( ParticleSpecies    , species    , INPUT        , REQUIRED );

  public:

    inline void execute() override final
    {
      ldbg << "Initializing AIREBO potential" << std::endl;

      const auto& p = *parameters;

      // --- REBO max range: max(rcmax[i][j]) over all species pairs ---
      double rcmax_max = 0.0;
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
          rcmax_max = std::max(rcmax_max, p.rcmax[i][j]);

      // bondorder_cutoff: BondOrderOp + NconjOp — REBO range only
      *bondorder_cutoff = rcmax_max;

      // rebo_cutoff: REBOForceOp — 2×rcmax so that j's REBO neighbours
      //   are always within i's buffer (r_il ≤ r_ij + r_jl ≤ 2×rcmax)
      *rebo_cutoff = 3.0 * rcmax_max;

      // lj_cutoff: LJForceOp — full LJ range max(cutlj×sigma[i][j])
      //   (also covers the 3/4-body path search: cut3rebo=3×rcmax ≤ lj_cutoff)
      if (p.ljflag) {
        double lj_max = 0.0;
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 2; ++j)
            lj_max = std::max(lj_max, p.cutlj * p.sigma[i][j]);
        *lj_cutoff = lj_max;
      } else {
        *lj_cutoff = 0.0;
      }

      // torsion_cutoff: TorsionForceOp — same geometry as REBO (k and l
      //   neighbours of i and j both within rcmax, so buf needs 2×rcmax)
      *torsion_cutoff = p.torflag ? 2.0 * rcmax_max : 0.0;

      // rcut_max: neighbour list construction cutoff = max of all active cutoffs
      *rcut_max = std::max({ *rcut_max,
                             *bondorder_cutoff,
                             *rebo_cutoff,
                             *lj_cutoff,
                             *torsion_cutoff });

      lout << "AIREBO bondorder_cutoff = " << *bondorder_cutoff << std::endl;
      lout << "AIREBO rebo_cutoff      = " << *rebo_cutoff      << std::endl;
      lout << "AIREBO lj_cutoff        = " << *lj_cutoff        << std::endl;
      lout << "AIREBO torsion_cutoff   = " << *torsion_cutoff   << std::endl;
      lout << "AIREBO rcut_max         = " << *rcut_max         << std::endl;
    }

  };

  ONIKA_AUTORUN_INIT(airebo_init)
  {
    OperatorNodeFactory::instance()->register_factory("airebo_init", make_simple_operator<AireboInit>);
  }

}
