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
#include <exanb/core/config.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>

#include <omp.h>

#include "include/mtp_params.h"
#include "include/mtp_config.h"

namespace exaStamp
{

  using namespace exanb;

  class MtpInit : public OperatorNode
  {
    ADD_SLOT( MtpParams        , parameters , INPUT        , REQUIRED );
    ADD_SLOT( double           , rcut_max   , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies  , species    , INPUT        , REQUIRED );
    ADD_SLOT( MtpContext       , mtp_ctx    , OUTPUT );

  public:

    inline void execute() override final
    {
      ldbg << "Initializing MTP potential" << std::endl;

      const auto& mtp_file = parameters->mtp_file;

      // Eagerly construct one EMTP per OpenMP thread, scratch fixed to MAX_PARTICLE_NEIGHBORS
      // (not grown on demand) -- see emtp.h.
      const int nt = omp_get_max_threads();
      ldbg << "Initializing " << nt << " MTP thread context(s)" << std::endl;
      mtp_ctx->m_emtp.resize(nt);
      for (int t = 0; t < nt; t++) {
        mtp_ctx->m_emtp[t] = std::make_shared<EMTP>(mtp_file, static_cast<int>(exanb::MAX_PARTICLE_NEIGHBORS));
      }

      auto& emtp = *mtp_ctx->m_emtp[0];

      *rcut_max = std::max(*rcut_max, emtp.rcut);
      ldbg << "MTP cutoff radius: " << emtp.rcut << std::endl;

      // Build exaStamp-type -> MTP-species mapping. PURELY POSITIONAL (see mtp_config.h): the
      // .mtp file carries no species names, so type i (in the declared species: block order)
      // maps to MTP species index i. Only the count can be validated -- log the resulting
      // mapping so a misordered species: block is at least visible.
      const auto& sp = *species;
      const int nspecies = sp.size();
      if (nspecies != emtp.species_count) {
        fatal_error() << "mtp_init: declared species count (" << nspecies << ") does not match "
                      << "species_count in MTP file " << mtp_file << " (" << emtp.species_count << ")" << std::endl;
      }
      mtp_ctx->nspecies = nspecies;
      mtp_ctx->type_map.resize(nspecies);
      for (int i = 0; i < nspecies; i++) {
        mtp_ctx->type_map[i] = i;
        ldbg << "Mapping atom type #" << i << " (" << sp[i].m_name
             << ") -> MTP species index #" << i
             << " (positional -- MTP files carry no species names; verify this order matches training)"
             << std::endl;
      }
    }

  };

  ONIKA_AUTORUN_INIT(mtp_init)
  {
    OperatorNodeFactory::instance()->register_factory("mtp_init", make_simple_operator<MtpInit>);
  }

}
