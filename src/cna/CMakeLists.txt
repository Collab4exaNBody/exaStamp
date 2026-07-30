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
# Common Neighbor Analysis (CNA), reimplemented from scratch (not vendored) following OVITO's
# own adaptive-cutoff algorithm (src/ovito/crystalanalysis/modifier/structureanalysis/
# StructureAnalysis.cpp -- the exact structure-identification engine DXA itself uses, per OVITO's
# own user manual: "This is done with the help of the Common Neighbor Analysis (CNA) method").
# No PTM dependency at all -- self-contained.

set(exaStampCNA_LINK_LIBRARIES exanbCore)
onika_add_plugin(exaStampCNA ${CMAKE_CURRENT_SOURCE_DIR})
