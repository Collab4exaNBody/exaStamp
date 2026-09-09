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

option(EXASTAMP_MLIP_POD_BUILD "Build POD potential" OFF)
set(EXASTAMP_MLIP_POD_BLAS_VENDOR "" CACHE STRING
    "BLA_VENDOR hint for the POD potential (e.g. OpenBLAS, Intel10_64lp, Generic). Empty = CMake auto-selects.")

# Mirror our own explicit toggle into onika_add_plugin's auto-generated one so it doesn't show up
# as a second, redundant, independently-settable ccmake entry.
set(ONIKA_BUILD_exaStampMlipPOD_DEFAULT ${EXASTAMP_MLIP_POD_BUILD})
set(ONIKA_BUILD_exaStampMlipPOD_FORCED ON)

if(EXASTAMP_MLIP_POD_BUILD)
  if(EXASTAMP_MLIP_POD_BLAS_VENDOR)
    set(BLA_VENDOR ${EXASTAMP_MLIP_POD_BLAS_VENDOR})
    message(STATUS "POD BLAS vendor: ${EXASTAMP_MLIP_POD_BLAS_VENDOR}")
  endif()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  add_library(eapod STATIC eapod.cpp)
  target_include_directories(eapod PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(eapod PUBLIC cxx_std_17)
  set_target_properties(eapod PROPERTIES POSITION_INDEPENDENT_CODE ON)
  target_link_libraries(eapod PUBLIC BLAS::BLAS LAPACK::LAPACK)

  set(exaStampMlipPOD_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors eapod)
endif()

onika_add_plugin(exaStampMlipPOD ${CMAKE_CURRENT_SOURCE_DIR})
