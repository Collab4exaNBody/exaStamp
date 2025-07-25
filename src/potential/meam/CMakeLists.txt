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

set(XSTAMP_MEAM_MAX_NEIGHBORS "32" CACHE STRING "MEAM max particle neighbors")
set(XSTAMP_MEAM_MULTIMAT_MAX_TYPES "4" CACHE STRING "MEAM/LJ max different particle types")
option(XSTAMP_MEAM_ENFORCE_OVERFLOW_CHECK "MEAM compute buffer overflow checking" OFF)
option(XSTAMP_MEAM_MERGE_FORCE_UPDATES "MEAM fusion force updates" ON)
option(XSTAMP_MEAM_USE_SHARED_MEM "MEAM uses more shared memory and less local variables" OFF)

set(exaStampMEAM_COMPILE_DEFINITIONS -DXSTAMP_MEAM_MAX_NEIGHBORS=${XSTAMP_MEAM_MAX_NEIGHBORS} -DXSTAMP_MEAM_MULTIMAT_MAX_TYPES=${XSTAMP_MEAM_MULTIMAT_MAX_TYPES})
  
if(XSTAMP_MEAM_MERGE_FORCE_UPDATES)
  set(exaStampMEAM_COMPILE_DEFINITIONS ${exaStampMEAM_COMPILE_DEFINITIONS} -DXSTAMP_MEAM_MERGE_FORCE_UPDATES=1)
endif()

if(XSTAMP_MEAM_ENFORCE_OVERFLOW_CHECK)
  set(exaStampMEAM_COMPILE_DEFINITIONS ${exaStampMEAM_COMPILE_DEFINITIONS} -DXSTAMP_MEAM_ENFORCE_OVERFLOW_CHECK=1)
endif()

if(XSTAMP_MEAM_USE_SHARED_MEM)
  set(exaStampMEAM_COMPILE_DEFINITIONS ${exaStampMEAM_COMPILE_DEFINITIONS} -DXSTAMP_MEAM_USE_SHARED_MEM=1)
endif()

set(exaStampMEAM_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors lennard_jones)

onika_add_plugin(exaStampMEAM ${CMAKE_CURRENT_SOURCE_DIR})

