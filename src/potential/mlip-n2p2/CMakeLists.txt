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

option(EXASTAMP_MLIP_N2P2_BUILD "Build N2P2 potential" OFF)
set(EXASTAMP_MLIP_N2P2_ROOT_DIR /usr/local/n2p2 CACHE STRING "N2P2 installation directory")

# Mirror our own explicit toggle into onika_add_plugin's auto-generated one so it doesn't show up
# as a second, redundant, independently-settable ccmake entry -- FORCED in both branches, so it
# always tracks EXASTAMP_MLIP_N2P2_BUILD (further forced OFF if the SDK directory isn't found).
if(NOT IS_DIRECTORY ${EXASTAMP_MLIP_N2P2_ROOT_DIR})
  if(EXASTAMP_MLIP_N2P2_BUILD)
    message(STATUS "N2P2 directory not found in ${EXASTAMP_MLIP_N2P2_ROOT_DIR}, N2P2 is forced OFF")
  endif()
  set(ONIKA_BUILD_exaStampMlipN2P2_FORCED ON)
  set(ONIKA_BUILD_exaStampMlipN2P2_DEFAULT OFF)
else()
  set(ONIKA_BUILD_exaStampMlipN2P2_FORCED ON)
  set(ONIKA_BUILD_exaStampMlipN2P2_DEFAULT ${EXASTAMP_MLIP_N2P2_BUILD})
endif()

# ${MPI_CXX_LIBRARIES}
set(microStampN2P2_INCLUDE_DIRS ${EXASTAMP_MLIP_N2P2_ROOT_DIR}/include)
set(microStampN2P2_LINK_LIBRARIES exaStampPotential ${EXASTAMP_MLIP_N2P2_ROOT_DIR}/lib/libnnp.so ${EXASTAMP_MLIP_N2P2_ROOT_DIR}/lib/libnnpif.so)
onika_add_plugin(exaStampMlipN2P2 ${CMAKE_CURRENT_SOURCE_DIR})

