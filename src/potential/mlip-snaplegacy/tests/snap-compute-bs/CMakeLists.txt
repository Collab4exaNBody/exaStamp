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

add_executable(snap-compute-bs-cpp snap-compute-bs.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyCG.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyGSH.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyBS.cpp)
target_include_directories(snap-compute-bs-cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_compile_options(snap-compute-bs-cpp PRIVATE -DDEBUG_SNAP_BS)

add_executable(snap-compute-bs-cpp-lammps snap-compute-bs-lammps.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyCG.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyGSH.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyBS.cpp)
target_include_directories(snap-compute-bs-cpp-lammps PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_compile_options(snap-compute-bs-cpp-lammps PRIVATE -DDEBUG_SNAP_BS -DLAMMPS=1)

file(COPY atom_positions.ref dbs.ref bs.ref2 Mo_conf.xyz bs_lammps.ref DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/snap-compute-bs.sh" "${CMAKE_CURRENT_BINARY_DIR}/snap-compute-bs.sh")
onika_add_test(utests-snap-compute-bs ${CMAKE_CURRENT_BINARY_DIR}/snap-compute-bs.sh)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/snap-compute-bs-lammps.sh" "${CMAKE_CURRENT_BINARY_DIR}/snap-compute-bs-lammps.sh")
onika_add_test(utests-snap-compute-bs-lammps ${CMAKE_CURRENT_BINARY_DIR}/snap-compute-bs-lammps.sh)

