# licensed to the apache software foundation (asf) under one
# or more contributor license agreements.  see the notice file
# distributed with this work for additional information
# regarding copyright ownership.  the asf licenses this file
# to you under the apache license, version 2.0 (the
# "license"); you may not use this file except in compliance
# with the license.  you may obtain a copy of the license at
# 
#   http://www.apache.org/licenses/license-2.0
# 
# unless required by applicable law or agreed to in writing,
# software distributed under the license is distributed on an
# "as is" basis, without warranties or conditions of any
# kind, either express or implied.  see the license for the
# specific language governing permissions and limitations
# under the license.

# ==========================================
# === Rigid molecule droplet generator   ===
# ==========================================

add_executable(genDropletMolRigid genDropletMolRigid.cpp StampV3LegacyIOStructures.cpp)
target_include_directories(genDropletMolRigid PUBLIC ${MPI_CXX_INCLUDE_PATH})
target_link_libraries(genDropletMolRigid ${MPI_CXX_LIBRARIES} ${OpenMP_CXX_LIBRARIES})
target_compile_options(genDropletMolRigid PUBLIC ${OpenMP_CXX_FLAGS})

install(TARGETS genDropletMolRigid DESTINATION bin)



