#
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements. See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership. The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#
# ======================================
# === DeePMD (native, no deepmd-kit) ===
# ======================================

# Per-atom scratch bounds (block shared memory on GPU, about 43 KB with the defaults)
set(XSTAMP_DEEPMD_MAX_NEIGHBORS "256" CACHE STRING "DeePMD max particle neighbors within the model cutoff")
set(XSTAMP_DEEPMD_MAX_TYPES "16" CACHE STRING "DeePMD max number of model types")
set(XSTAMP_DEEPMD_MAX_EMB_WIDTH "128" CACHE STRING "DeePMD max embedding network width")
set(XSTAMP_DEEPMD_MAX_EMB_LAYERS "4" CACHE STRING "DeePMD max embedding network layers")
set(XSTAMP_DEEPMD_MAX_FIT_WIDTH "256" CACHE STRING "DeePMD max fitting network width")
set(XSTAMP_DEEPMD_MAX_FIT_LAYERS "4" CACHE STRING "DeePMD max fitting network layers, output layer included")
set(XSTAMP_DEEPMD_MAX_DESCRIPTOR "2048" CACHE STRING "DeePMD max descriptor size (embedding width x axis_neuron)")
set(exaStampMlipDeepmd_COMPILE_DEFINITIONS
    -DXSTAMP_DEEPMD_MAX_NEIGHBORS=${XSTAMP_DEEPMD_MAX_NEIGHBORS}
    -DXSTAMP_DEEPMD_MAX_TYPES=${XSTAMP_DEEPMD_MAX_TYPES}
    -DXSTAMP_DEEPMD_MAX_EMB_WIDTH=${XSTAMP_DEEPMD_MAX_EMB_WIDTH}
    -DXSTAMP_DEEPMD_MAX_EMB_LAYERS=${XSTAMP_DEEPMD_MAX_EMB_LAYERS}
    -DXSTAMP_DEEPMD_MAX_FIT_WIDTH=${XSTAMP_DEEPMD_MAX_FIT_WIDTH}
    -DXSTAMP_DEEPMD_MAX_FIT_LAYERS=${XSTAMP_DEEPMD_MAX_FIT_LAYERS}
    -DXSTAMP_DEEPMD_MAX_DESCRIPTOR=${XSTAMP_DEEPMD_MAX_DESCRIPTOR})

set(exaStampMlipDeepmd_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)

# DPA3 host matrix products through cblas_dgemm (Accelerate on macOS, else a BLAS found by CMake), else built-in loops
option(EXASTAMP_DEEPMD_BLAS "DeePMD DPA3 host matrix products through cblas_dgemm" ${APPLE})
if(EXASTAMP_DEEPMD_BLAS)
  if(APPLE)
    list(APPEND exaStampMlipDeepmd_LINK_LIBRARIES "-framework Accelerate")
  else()
    find_package(BLAS REQUIRED)
    list(APPEND exaStampMlipDeepmd_LINK_LIBRARIES ${BLAS_LIBRARIES})
  endif()
  list(APPEND exaStampMlipDeepmd_COMPILE_DEFINITIONS -DXSTAMP_DEEPMD_BLAS)
else()
  message(WARNING "EXASTAMP_DEEPMD_BLAS is OFF: DeePMD DPA2, DPA3 and DPA4 matrix products use built-in loops, about 7 times slower for DPA4. "
                  "Enable it with a BLAS for production runs.")
endif()
onika_add_plugin(exaStampMlipDeepmd ${CMAKE_CURRENT_SOURCE_DIR})
