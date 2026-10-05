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

option(EXASTAMP_MLIP_PACE_BUILD "Build PACE potential" OFF)
set(EXASTAMP_MLIP_PACE_GIT_REPO "git@github.com:Collab4exaNBody/exaStamp_mlip_pace.git" CACHE STRING "Git repository URL for PACE compute libraries")
set(EXASTAMP_MLIP_PACE_GIT_TAG "main" CACHE STRING "Git tag/branch for PACE compute libraries")
set(EXASTAMP_MLIP_PACE_EXTERNAL_DIR "${CMAKE_BINARY_DIR}/external/pace")

# Mirror our own explicit toggle into onika_add_plugin's auto-generated one so it doesn't show up
# as a second, redundant, independently-settable ccmake entry.
set(ONIKA_BUILD_exaStampMlipPace_DEFAULT ${EXASTAMP_MLIP_PACE_BUILD})
set(ONIKA_BUILD_exaStampMlipPace_FORCED ON)

if(EXASTAMP_MLIP_PACE_BUILD)
  if(NOT EXISTS "${EXASTAMP_MLIP_PACE_EXTERNAL_DIR}/.git")
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/external")
    message(STATUS "Cloning PACE libraries from ${EXASTAMP_MLIP_PACE_GIT_REPO} (${EXASTAMP_MLIP_PACE_GIT_TAG}) into ${EXASTAMP_MLIP_PACE_EXTERNAL_DIR}")
    execute_process(
      COMMAND git clone --branch ${EXASTAMP_MLIP_PACE_GIT_TAG} ${EXASTAMP_MLIP_PACE_GIT_REPO} pace
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/external"
      RESULT_VARIABLE _pace_clone_result
    )
    if(NOT _pace_clone_result EQUAL 0)
      message(FATAL_ERROR "Failed to clone PACE repository (exit code: ${_pace_clone_result})")
    endif()
  else()
    message(STATUS "PACE libraries already present in ${EXASTAMP_MLIP_PACE_EXTERNAL_DIR}")
  endif()

  # YAML_CPP_INSTALL_DIR is a general onika-ecosystem variable, set at the top level, forwarded
  # implicitly to this external add_subdirectory (not PACE-specific -- do not rename/relocate it).
  add_subdirectory(${EXASTAMP_MLIP_PACE_EXTERNAL_DIR} ${CMAKE_CURRENT_BINARY_DIR}/pace-mlip)
  set(exaStampMlipPace_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors pace::pace)
endif()

onika_add_plugin(exaStampMlipPace ${CMAKE_CURRENT_SOURCE_DIR})
