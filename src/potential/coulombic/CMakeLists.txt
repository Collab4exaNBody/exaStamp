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

set(exaStampCoulombic_LINK_LIBRARIES exaStampCompute exaStampPotential exanbParticleNeighbors)
# pair kernels are header only and defined once, next to the pair potentials using them
set(exaStampCoulombic_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}/../pair_potentials/coul_wolf_pair/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../pair_potentials/coul_dsf/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../reaction_field/include)
onika_add_plugin(exaStampCoulombic ${CMAKE_CURRENT_SOURCE_DIR})

