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
# ===== exaStamp PACE potential  =====
# ====================================

set(PACE_ROOT_DIR /usr/local/pace CACHE STRING "PACE installation directory")

if(IS_DIRECTORY ${PACE_ROOT_DIR})
  message(STATUS "PACE directory found in ${PACE_ROOT_DIR}, PACE is forced ON")
  set(ONIKA_BUILD_exaStampPace_FORCED ON)
else()
  message(STATUS "PACE directory not found in ${PACE_ROOT_DIR}, PACE is forced OFF")
  set(ONIKA_BUILD_exaStampPace_FORCED OFF)
endif()

set(ONIKA_BUILD_exaStampPace ${ONIKA_BUILD_exaStampPace_FORCED} CACHE BOOL "Build exaStampPace")

if(ONIKA_BUILD_exaStampPace)
  set(WIGNER_PATH ${PACE_ROOT_DIR}/wigner-cpp)
  set(WIGNER_INCLUDE_PATH ${WIGNER_PATH}/include)
  
  add_library(cnpy SHARED ${PACE_ROOT_DIR}/cnpy/cnpy.cpp)  
  target_include_directories(cnpy PUBLIC ${PACE_ROOT_DIR}/cnpy)

  file(GLOB exaStampPace_EXTERNAL_SRCS
    ${PACE_ROOT_DIR}/ML-PACE/ace-evaluator/*.cpp
    ${PACE_ROOT_DIR}/ML-PACE/ace/*.cpp
  )
  add_library(pace STATIC ${exaStampPace_EXTERNAL_SRCS})
  target_include_directories(pace PUBLIC ${WIGNER_INCLUDE_PATH})
  target_include_directories(pace PUBLIC ${PACE_ROOT_DIR}/ML-PACE)
  target_compile_definitions(pace PUBLIC EXTRA_C_PROJECTIONS)
  target_compile_definitions(pace PUBLIC COMPUTE_B_GRAD)
  target_link_libraries(pace PUBLIC cnpy)

  set_target_properties(pace PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
  )
  
  set(exaStampPace_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors cnpy pace)
  onika_add_plugin(exaStampPace ${CMAKE_CURRENT_SOURCE_DIR})
endif()
