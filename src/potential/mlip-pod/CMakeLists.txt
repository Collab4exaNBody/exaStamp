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
# ===== exaStamp POD potential   =====
# ====================================

# BLAS and LAPACK are found in the top-level CMakeLists when EXASTAMP_BUILD_POD=ON.
add_library(eapod STATIC eapod.cpp)
target_include_directories(eapod PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(eapod PUBLIC cxx_std_17)
set_target_properties(eapod PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(eapod PUBLIC BLAS::BLAS LAPACK::LAPACK)

set(exaStampPOD_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors eapod)
onika_add_plugin(exaStampPOD ${CMAKE_CURRENT_SOURCE_DIR})
