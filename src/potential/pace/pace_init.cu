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

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"
#include "ace/ace_b_basis.h"

#include "pace_params.h"
#include "pace_config.h"

namespace exaStamp
{

  using namespace exanb;

  class PaceInit : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( PaceParams      , parameters , INPUT        , REQUIRED );
    ADD_SLOT( double          , rcut_max   , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies , species    , INPUT        , REQUIRED );
    ADD_SLOT( PaceContext     , pace_ctx   , OUTPUT );

  public:

    inline void execute () override final
    {
      ldbg << "Initializing ACE potential" << std::endl;

      pace_ctx->recursive = parameters->recursive;
      pace_ctx->aceimpl   = std::make_unique<ACEImpl>();
      pace_ctx->aceimpl->basis_set = std::make_unique<ACECTildeBasisSet>();
      pace_ctx->aceimpl->ace       = std::make_unique<ACERecursiveEvaluator>();

      const auto& potential_file_name = parameters->pace_coef;
      const auto& sp = *species;
      const int nspecies = sp.size();
      pace_ctx->nspecies = nspecies;

      if (hasExtension(potential_file_name, ".yaml")) {
        ACEBBasisSet bBasisSet(potential_file_name);
        *pace_ctx->aceimpl->basis_set = bBasisSet.to_ACECTildeBasisSet();
      } else {
        *pace_ctx->aceimpl->basis_set = ACECTildeBasisSet(potential_file_name);
      }

      pace_ctx->aceimpl->ace->set_recursive(pace_ctx->recursive);
      pace_ctx->aceimpl->ace->element_type_mapping.init(nspecies + 1);

      for (int i = 1; i <= nspecies; i++) {
        const char *elemname = sp[i-1].m_name;
        if (atomicNumberByName(elemname) == -1) {
          lerr << "Element " << elemname << " is not a valid chemical element" << std::endl;
        }
        SPECIES_TYPE mu = pace_ctx->aceimpl->basis_set->get_species_index_by_name(elemname);
        if (mu != -1) {
          ldbg << "Mapping atom type #" << i << " (" << elemname << ") -> ACE species type #" << mu << std::endl;
          pace_ctx->aceimpl->ace->element_type_mapping(i) = mu;
        } else {
          fatal_error() << "Element " << elemname << " is not supported by ACE potential from file " << potential_file_name << std::endl;
        }
      }
      pace_ctx->aceimpl->ace->set_basis(*pace_ctx->aceimpl->basis_set, 1);

      for (int i = 0; i < nspecies; i++)
        for (int j = 0; j < nspecies; j++)
          *rcut_max = std::max( *rcut_max, pace_ctx->aceimpl->basis_set->radial_functions->cut(i,j) );
    }

  };

  ONIKA_AUTORUN_INIT(pace_init)
  {
    OperatorNodeFactory::instance()->register_factory( "pace_init", make_simple_operator<PaceInit> );
  }

}
