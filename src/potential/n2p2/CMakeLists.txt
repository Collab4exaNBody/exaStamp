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

set(N2P2_ROOT_DIR /usr/local/n2p2 CACHE STRING "N2P2 installation directory")

# default is disabled, user must explicitly enable it
set(ONIKA_BUILD_exaStampN2P2_DEFAULT OFF)

if(NOT IS_DIRECTORY ${N2P2_ROOT_DIR})
	message(STATUS "N2P2 directory not found in ${N2P2_ROOT_DIR}, N2P2 is forced OFF")
  set(ONIKA_BUILD_exaStampN2P2_FORCED ON)
else()
  set(ONIKA_BUILD_exaStampN2P2_FORCED OFF)
endif()

# ${MPI_CXX_LIBRARIES}
set(microStampN2P2_INCLUDE_DIRS ${N2P2_ROOT_DIR}/include)
set(microStampN2P2_LINK_LIBRARIES exaStampPotential ${N2P2_ROOT_DIR}/lib/libnnp.so ${N2P2_ROOT_DIR}/lib/libnnpif.so)
onika_add_plugin(exaStampN2P2 ${CMAKE_CURRENT_SOURCE_DIR})

