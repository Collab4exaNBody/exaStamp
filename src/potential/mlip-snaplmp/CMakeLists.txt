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

option(EXASTAMP_MLIP_SNAPLMP_BUILD "Build SnapLMP potential" OFF)
set(EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR "${XSTAMP_THIRD_PARTY_TOOLS_ROOT}/lammps" CACHE STRING "LAMMPS source directory")

# Mirror our own explicit toggle into onika_add_plugin's auto-generated one so it doesn't show up
# as a second, redundant, independently-settable ccmake entry -- FORCED in both branches, so it
# always tracks EXASTAMP_MLIP_SNAPLMP_BUILD (further forced OFF if the LAMMPS source isn't found,
# or if exaNBody's mdSnap contrib -- which exaStampSnapLMP_LINK_LIBRARIES links against below --
# isn't available).
if(NOT IS_DIRECTORY ${EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR} OR NOT XNB_BUILD_mdSnap)
  if(EXASTAMP_MLIP_SNAPLMP_BUILD)
    if(NOT IS_DIRECTORY ${EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR})
      message(STATUS "LAMMPS directory not found in ${EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR}, SnapLMP is forced OFF")
    endif()
    if(NOT XNB_BUILD_mdSnap)
      message(STATUS "exaNBody's mdSnap contrib not available, SnapLMP is forced OFF")
    endif()
  endif()
  set(ONIKA_BUILD_exaStampSnapLMP_FORCED ON)
  set(ONIKA_BUILD_exaStampSnapLMP_DEFAULT OFF)
else()
  set(ONIKA_BUILD_exaStampSnapLMP_FORCED ON)
  set(ONIKA_BUILD_exaStampSnapLMP_DEFAULT ${EXASTAMP_MLIP_SNAPLMP_BUILD})
endif()

set(exaStampSnapLMP_EXTERNAL_SRCS ${EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR}/src/memory.cpp ${EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR}/src/math_special.cpp ${EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR}/src/ML-SNAP/sna.cpp)
set(exaStampSnapLMP_LINK_LIBRARIES mdSnap exaStampPotential exanbParticleNeighbors)
set(exaStampSnapLMP_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR}/src/ML-SNAP ${EXASTAMP_MLIP_SNAPLMP_LMP_SRC_DIR}/src)
onika_add_plugin(exaStampSnapLMP ${CMAKE_CURRENT_SOURCE_DIR})

