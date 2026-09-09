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

add_executable(snap-compute-cg-cpp snap-compute-cg.cpp ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/SnapLegacyCG.cpp)
target_include_directories(snap-compute-cg-cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../include)
target_compile_options(snap-compute-cg-cpp PRIVATE -DDEBUG_SNAP_CG)

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/snap-compute-cg.sh" "${CMAKE_CURRENT_BINARY_DIR}/snap-compute-cg.sh")
onika_add_test(unit-snap-compute-cg ${CMAKE_CURRENT_BINARY_DIR}/snap-compute-cg.sh)

