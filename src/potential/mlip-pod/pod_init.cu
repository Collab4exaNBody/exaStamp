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

#ifdef EAPOD_USE_OPENBLAS
extern "C" void openblas_set_num_threads(int);
#endif

#include "pod_params.h"
#include "pod_config.h"

namespace exaStamp
{

  using namespace exanb;

  class PodInit : public OperatorNode
  {
    ADD_SLOT( PodParams       , parameters , INPUT        , REQUIRED );
    ADD_SLOT( double          , rcut_max   , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies , species    , INPUT        , REQUIRED );
    ADD_SLOT( PodContext      , pod_ctx    , OUTPUT );

  public:

    inline void execute() override final
    {
      ldbg << "Initializing POD potential" << std::endl;

      const auto& pod_file   = parameters->pod_file;
      const auto& coeff_file = parameters->coeff_file;

      // Store paths for lazy per-thread EAPOD construction in pod_force.
      pod_ctx->pod_file   = pod_file;
      pod_ctx->coeff_file = coeff_file;

#ifdef EAPOD_USE_OPENBLAS
      // Prevent OpenBLAS from spawning its own threads inside the OMP parallel
      // region — each exaStamp thread calls BLAS sequentially for its own atom.
      openblas_set_num_threads(1);
      ldbg << "OpenBLAS: set to single-threaded mode (OMP handles parallelism)" << std::endl;
#endif

      // Eagerly construct one EAPOD per OpenMP thread.
      // Done here at init time so execute() carries zero thread-setup overhead.
      const int nt = omp_get_max_threads();
      ldbg << "Initializing " << nt << " POD thread context(s)" << std::endl;
      pod_ctx->m_eapod.resize(nt);
      for (int t = 0; t < nt; t++) {
        auto ep = std::make_shared<EAPOD>(pod_file, coeff_file);
        ep->allocate_temp_memory(ep->Njmax);
        pod_ctx->m_eapod[t] = ep;
      }

      // Use thread-0 instance for metadata queries.
      auto& eapod = *pod_ctx->m_eapod[0];

      // Set cutoff radius.
      *rcut_max = std::max(*rcut_max, eapod.rcut);
      ldbg << "POD cutoff radius: " << eapod.rcut << std::endl;

      // Build exaStamp-type → POD-type mapping.
      const auto& sp = *species;
      const int nspecies = sp.size();
      pod_ctx->nspecies = nspecies;
      pod_ctx->type_map.assign(nspecies, -1);

      for (int i = 0; i < nspecies; i++) {
        const std::string name(sp[i].m_name);
        bool found = false;
        for (int j = 0; j < (int)eapod.species.size(); j++) {
          if (eapod.species[j] == name) {
            pod_ctx->type_map[i] = j + 1;
            ldbg << "Mapping atom type #" << i << " (" << name
                 << ") -> POD species type #" << (j+1) << std::endl;
            found = true;
            break;
          }
        }
        if (!found) {
          fatal_error() << "Element " << name
                        << " is not supported by POD potential from file "
                        << pod_file << std::endl;
        }
      }
    }

  };

  ONIKA_AUTORUN_INIT(pod_init)
  {
    OperatorNodeFactory::instance()->register_factory("pod_init", make_simple_operator<PodInit>);
  }

}
