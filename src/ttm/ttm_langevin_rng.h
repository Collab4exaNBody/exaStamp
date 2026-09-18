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
#pragma once

#include <exanb/grid_cell_particles/particle_random_selection.h> // splitmix64, particle_random_key
#include <onika/cuda/cuda.h>
#include <cmath>

namespace exaStamp
{
  // Stateless Langevin coupling noise, replacing the old thread-indexed onika::parallel::random_engine()
  // (which was never even reproducible run-to-run, since draw order followed schedule(dynamic) thread
  // timing): same (particle id, MD step, axis) always yields the same draw, independent of thread count
  // or scheduling -- and, being a pure function of its arguments, this is GPU-portable as-is (reused
  // unchanged by the Pass B GPU port).
  //
  // salt selects which of a particle's 3 independent per-axis draws (0=x, 1=y, 2=z) this call produces.

  // uniform draw in [-0.5,0.5), matching LAMMPS fix-ttm's own noise convention (lammps_noise: true).
  ONIKA_HOST_DEVICE_FUNC inline double ttm_langevin_uniform_rand( uint64_t p_id, uint64_t md_step, int salt )
  {
    return exanb::particle_random_key( p_id, md_step * 8 + uint64_t(salt) ) - 0.5;
  }

  // standard-normal draw (mean 0, variance 1), this codebase's own default noise convention.
  // Box-Muller from 2 independent uniform draws (cos branch only -- the sin branch is simply unused,
  // no need for the second gaussian a full Box-Muller pair would give).
  ONIKA_HOST_DEVICE_FUNC inline double ttm_langevin_gauss_rand( uint64_t p_id, uint64_t md_step, int salt )
  {
    const double u1 = exanb::particle_random_key( p_id, md_step * 8 + uint64_t(salt) * 2     );
    const double u2 = exanb::particle_random_key( p_id, md_step * 8 + uint64_t(salt) * 2 + 1 );
    const double r = std::sqrt( -2.0 * std::log( u1 > 1.e-300 ? u1 : 1.e-300 ) );
    return r * std::cos( 2.0 * M_PI * u2 );
  }
}
