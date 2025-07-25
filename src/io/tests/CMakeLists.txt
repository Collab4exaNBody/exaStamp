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
# ==========================================
# === xyz2legacy : xyz to StampV3 dump   ===
# ==========================================
add_executable(xyz2legacy xyz2legacy.cpp)
target_include_directories(xyz2legacy PUBLIC ${PROJECT_SOURCE_DIR}/src/exaStamp/io/include ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(xyz2legacy exaStampIO ${MPI_CXX_LIBRARIES})
install(TARGETS xyz2legacy DESTINATION bin)

# =============================================
# === readDumpStampV3 : StampV3 dump viewer ===
# =============================================
add_executable(readDumpStampV3 readDumpStampV3.cpp)
target_include_directories(readDumpStampV3 PUBLIC ${PROJECT_SOURCE_DIR}/src/exaStamp/io/include ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(readDumpStampV3 exaStampIO ${MPI_CXX_LIBRARIES})
install(TARGETS readDumpStampV3 DESTINATION bin)

# ==========================================
# === XSP2LAMMPS : a data file converter ===
# ==========================================
if(XNB_BUILD_CONTRIBS) # because we need saru random generator by now
  add_executable(genDroplet genDroplet.cpp)
  target_include_directories(genDroplet PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/src/exaStamp/io/include ${MPI_CXX_INCLUDE_PATH})
  target_link_libraries(genDroplet exaStampIO SARU ${MPI_CXX_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
  target_compile_options(genDroplet PUBLIC ${OpenMP_CXX_FLAGS})
  install(TARGETS genDroplet DESTINATION bin)
endif()


# ==========================================
# === XSP2LAMMPS : a data file converter ===
# ==========================================
add_executable(stampv3_to_lammps stampv3_to_lammps.cpp)
target_include_directories(stampv3_to_lammps PUBLIC ${PROJECT_SOURCE_DIR}/src/exaStamp/io/include ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(stampv3_to_lammps exaStampIO ${MPI_CXX_LIBRARIES})
install(TARGETS stampv3_to_lammps DESTINATION bin)

add_test(stampv3_to_lammps_par ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/stampv3_to_lammps ${CMAKE_CURRENT_SOURCE_DIR}/StampV3.MpiIO ${CMAKE_CURRENT_BINARY_DIR}/MpiIO2lammps_par_lammps.dat ${MPIEXEC_POSTFLAGS} )
add_test(stampv3_to_lammps_mono ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/stampv3_to_lammps ${CMAKE_CURRENT_SOURCE_DIR}/StampV3.MpiIO ${CMAKE_CURRENT_BINARY_DIR}/MpiIO2lammps_mono_lammps.dat ${MPIEXEC_POSTFLAGS} )

