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
# ====================================================================
# === YE3T linear ACE (.yace), native onika kernels (CPU and GPU) ===
# ====================================================================

option(EXASTAMP_MLIP_YE3T_BUILD "Build the YE3T linear ACE potential (native kernels, no external library)" ON)
set(ONIKA_BUILD_exaStampMlipYe3t_DEFAULT ${EXASTAMP_MLIP_YE3T_BUILD})
set(ONIKA_BUILD_exaStampMlipYe3t_FORCED ON)

set(exaStampMlipYe3t_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../mlip-utils/include)
set(exaStampMlipYe3t_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
onika_add_plugin(exaStampMlipYe3t ${CMAKE_CURRENT_SOURCE_DIR})
