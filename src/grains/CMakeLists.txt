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
# Grain segmentation (misorientation-based clustering of PTM-classified atoms into grains),
# following OVITO's crystalanalysis/modifier/grains modifier: Node-Pair-Sampling graph clustering +
# Dijkstra-by-distance orphan atom adoption. Consumes ptm_type/ptm_orientation fields from the ptm
# plugin -- needs its include dir for ptm_constants.h's PTM_MATCH_* enum only (no link dependency,
# the disorientation math itself is duplicated in device-compatible form, see
# compute_grain_bond_misorientation.cu's own header comment).

set(exaStampGrains_LINK_LIBRARIES exaStampCompute exaStampThermoState exanbDefBox exanbIO exanbParticleNeighbors exaStampParticleSpecies)
set(exaStampGrains_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../ptm/include)
onika_add_plugin(exaStampGrains ${CMAKE_CURRENT_SOURCE_DIR})
