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

# ======================================================= #
# ==================== Pair Potentials ================== #
add_subdirectory(pair_potential_template)
add_subdirectory(pair_potentials)
add_subdirectory(multimat_pair_potential)
# ======================================================= #

# ======================================================= #
# ============ Embedded Atom Model Potentials =========== #
add_subdirectory(eam)
add_subdirectory(eam_potential_template)
add_subdirectory(eam_potentials)
# ======================================================= #

# ======================================================= #
# ================== Reactive Potentials ================ #
add_subdirectory(LCHBOP)
# ======================================================= #

# ======================================================= #
# ======== Modified Embedded Atom Model Potentials ====== #
add_subdirectory(meam)
# ======================================================= #

# ======================================================= #
# ================= Coulombic Potentials ================ #
add_subdirectory(coul_dsf)
add_subdirectory(coul_wolf)
add_subdirectory(ewald)
add_subdirectory(reaction_field)
# ======================================================= #

# ======================================================= #
# ================== Combined Potentials ================ #
add_subdirectory(ljexp6rf)
# ======================================================= #

# ======================================================= #
# ======= Machine Learning Interatomic Potentials ======= #
if(EXASTAMP_BUILD_PACE)
  add_subdirectory(mlip-pace)
endif()

if(EXASTAMP_BUILD_POD)
  add_subdirectory(mlip-pod)
endif()

add_subdirectory(n2p2)

add_subdirectory(snap)
add_subdirectory(snaplegacy)
add_subdirectory(snaplmp)
# ======================================================= #

# ======================================================= #
# =============== Experimental Potentials =============== #
add_subdirectory(range_potential)
# ======================================================= #

