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
# =================================
# === micro-stamp core library  ===
# =================================
find_package(ZLIB REQUIRED)

set(exaStampIO_LINK_LIBRARIES exaStampCompute exanbMPI exanbIO exanbDefBox exaStampParticleSpecies exanbGridCellParticles exanbParticleNeighbors ${ZLIB_LIBRARIES} stdc++fs)
set(exaStampIO_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR} ${BASEN_INCLUDE_DIRS}) 
onika_add_plugin(exaStampIO ${CMAKE_CURRENT_SOURCE_DIR})

