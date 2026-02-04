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
# message(STATUS "XNB_BUILD_mdSnap = ${XNB_BUILD_mdSnap}")
if(XNB_BUILD_mdSnap)
  set(exaStampSnapLegacy_COMPILE_DEFINITIONS -DLAMMPS=1 ${KMP_ALIGNED_ALLOCATOR_DEFINITIONS})
  set(exaStampSnapLegacy_INCLUDE_DIRS ${SNAP_INCLUDE_DIRS})
  set(exaStampSnapLegacy_LINK_LIBRARIES mdSnap exanbParticleNeighbors ${SNAP_LIBRARIES})
  onika_add_plugin(exaStampSnapLegacy ${CMAKE_CURRENT_SOURCE_DIR})
endif()

