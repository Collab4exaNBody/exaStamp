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

#include <onika/math/basic_types.h>
#include <vector>
#include <cstdint>

namespace exaStamp
{
  using namespace exanb;

  // Output of grain_segmentation_nps (see grain_segmentation_algo.cpp) -- one entry per LOCAL atom
  // (owned+ghost, same indexing the caller's flattened struct_type/orientation/bond arrays use).
  struct GrainSegmentationResult
  {
    std::vector<int32_t> atom_grain_id;     // 1-based, 0 = not part of any grain (dissolved/never merged/unadopted orphan)
    std::vector<long>    grain_size;        // index 0 -> grain id 1, atom count
    std::vector<double>  grain_orientation; // 4 doubles/grain (w,x,y,z), normalized mean quaternion
    std::vector<int>     grain_structure_type; // PTM_MATCH_* code, same for every member atom by construction
    std::vector<Vec3d>   grain_color;       // RGB in [0,1]
    long n_grains = 0;
    double merge_threshold_log = 0.0; // the (possibly auto-selected) log-distance cutoff actually applied
    long n_orphans_adopted = 0;
  };
}
