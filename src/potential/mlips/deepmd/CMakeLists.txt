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
# ===================================
# === micro-stamp meam potential  ===
# ===================================

# set(DeepMD_ROOT_DIR /usr/local/deepmd CACHE STRING "DEEPMD installation directory")

# # # default is disabled, user must explicitly enable it
# # set(ONIKA_BUILD_exaStampDeepMD_DEFAULT OFF)

# # if(NOT IS_DIRECTORY ${DeepMD_ROOT_DIR})
# # 	message(STATUS "DeepMD directory not found in ${DeepMD_ROOT_DIR}, DeepMD is forced OFF")
# #   set(ONIKA_BUILD_exaStampDeepMD_FORCED OFF)
# # else()
# #   set(ONIKA_BUILD_exaStampDeepMD_FORCED ON)
# # endif()


# set(exaStampDeepMD_LINK_LIBRARIES exaStampPotential exaStampParticleSpecies exanbGridCellParticles exanbCompute ${DeepMD_ROOT_DIR}/lib/libdeepmd_cc_low.so)
# onika_add_plugin(exaStampDeepMD ${CMAKE_CURRENT_SOURCE_DIR})

# # set(exaStampDEEPMD_LINK_LIBRARIES exaStampPotential exaStampParticleSpecies)
# # onika_add_plugin(exaStampDEEPMD ${CMAKE_CURRENT_SOURCE_DIR})

# TODO with the following

# Add the following to CMake 
# -DDeePMD_DIR=/home/lafourcadep/CODES/ATOMISTIC/FFs/deepmd-kit/lib/cmake/DeePMD

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

#set(exaStampDeepMD_COMPILE_DEFINITIONS TORCH_STABLE_ONLY=1)
find_package(DeePMD REQUIRED)
set(exaStampDeepMD_INCLUDE_DIRS ${DeepMD_ROOT_DIR}/include/deepmd ${TORCH_ROOT_DIR}/include/torch)
set(exaStampDeepMD_LINK_LIBRARIES exaStampPotential exaStampParticleSpecies exanbGridCellParticles exanbCompute PRIVATE DeePMD::deepmd_cc "${TORCH_LIBRARIES}")
set(CMAKE_CXX_STANDARD 17)
  
onika_add_plugin(exaStampDeepMD ${CMAKE_CURRENT_SOURCE_DIR})
