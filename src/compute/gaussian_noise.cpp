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
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <onika/parallel/random.h>

#include <exanb/grid_cell_particles/particle_region.h>
#include <exanb/compute/gaussian_noise.h>

namespace exaStamp
{
  using namespace exanb;

  template<class GridT> using GaussianNoiseR = GaussianNoise < GridT , field::_id , FieldSet<field::_rx,field::_ry,field::_rz> >;
  template<class GridT> using GaussianNoiseV = GaussianNoise < GridT , field::_id , FieldSet<field::_vx,field::_vy,field::_vz> >;
  template<class GridT> using GaussianNoiseF = GaussianNoise < GridT , field::_id , FieldSet<field::_fx,field::_fy,field::_fz> >;

  // === register factories ===
  ONIKA_AUTORUN_INIT(gaussian_noise)
  {
   OperatorNodeFactory::instance()->register_factory( "gaussian_noise_r", make_grid_variant_operator< GaussianNoiseR > );
   OperatorNodeFactory::instance()->register_factory( "gaussian_noise_v", make_grid_variant_operator< GaussianNoiseV > );
   OperatorNodeFactory::instance()->register_factory( "gaussian_noise_f", make_grid_variant_operator< GaussianNoiseF > );
  }

}

