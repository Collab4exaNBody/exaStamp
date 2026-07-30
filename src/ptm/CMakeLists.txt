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
# Vendored Polyhedral Template Matching (PTM, MIT license, see PTM_LICENSE) from
# https://github.com/pmla/polyhedral-template-matching (as bundled by OVITO, which is where these
# sources were copied from). lib/ptm_*.cpp compile into an internal exaStampPTM shared library
# that compute_ptm.cpp (this plugin's own operator) links against automatically.

set(exaStampPTM_LINK_LIBRARIES exanbCore)
# compute_dxa_edge_vectors.cpp needs DelaunayTessellation/DXAEdgeVectors's struct definitions
# (header-only, no link dependency on Geogram or the delaunay plugin's own .so) -- src/delaunay's
# include dir isn't on any other plugin's include path by default, so add it explicitly.
set(exaStampPTM_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../delaunay/include)
onika_add_plugin(exaStampPTM ${CMAKE_CURRENT_SOURCE_DIR})
