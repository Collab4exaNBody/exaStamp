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
# ====================================
# ===== exaStamp MTP potential   =====
# ====================================

option(EXASTAMP_MLIP_MTP_BUILD "Build MTP potential" OFF)

# Mirror our own explicit toggle into onika_add_plugin's auto-generated one so it doesn't show up
# as a second, redundant, independently-settable ccmake entry.
set(ONIKA_BUILD_exaStampMTP_DEFAULT ${EXASTAMP_MLIP_MTP_BUILD})
set(ONIKA_BUILD_exaStampMTP_FORCED ON)

# Unlike POD, the core MTP math needs no BLAS/LAPACK (confirmed against the LAMMPS ML-MTP
# reference implementation -- zero such dependency; the out-of-scope MaxVol active-learning
# feature would need one, but even that's precomputed offline by MLIP-3, never inverted here) --
# so this static lib is cheap and safe to build unconditionally, unlike POD's eapod.
add_library(emtp STATIC emtp.cpp)
target_include_directories(emtp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(emtp PUBLIC cxx_std_17)
set_target_properties(emtp PROPERTIES POSITION_INDEPENDENT_CODE ON)

set(exaStampMTP_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors emtp)
onika_add_plugin(exaStampMTP ${CMAKE_CURRENT_SOURCE_DIR})
