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
# =================================
# === micro-stamp core library  ===
# =================================

set(XSTAMP_MAX_MOLECULE_ATOMS "24" CACHE STRING "Maximum number of atoms in a molecule")
set(XSTAMP_MAX_MOLECULE_BONDS "24" CACHE STRING "Maximum number of atom bonds in a molecule")
set(XSTAMP_MAX_MOLECULE_BENDS "24" CACHE STRING "Maximum number of angles in a molecule")
set(XSTAMP_MAX_MOLECULE_TORSIONS "24" CACHE STRING "Maximum number of torsions in a molecule")
set(XSTAMP_MAX_MOLECULE_IMPROPERS "24" CACHE STRING "Maximum number of imprompers in a molecule")
set(XSTAMP_MAX_MOLECULE_PAIRS "32" CACHE STRING "Maximum number of distant atom pairs in a molecule")

set(exaStampMolecule_COMPILE_DEFINITIONS
  -DXSTAMP_MAX_MOLECULE_ATOMS=${XSTAMP_MAX_MOLECULE_ATOMS}
  -DXSTAMP_MAX_MOLECULE_BONDS=${XSTAMP_MAX_MOLECULE_BONDS}
  -DXSTAMP_MAX_MOLECULE_BENDS=${XSTAMP_MAX_MOLECULE_BENDS}
  -DXSTAMP_MAX_MOLECULE_TORSIONS=${XSTAMP_MAX_MOLECULE_TORSIONS}
  -DXSTAMP_MAX_MOLECULE_IMPROPERS=${XSTAMP_MAX_MOLECULE_IMPROPERS}
  -DXSTAMP_MAX_MOLECULE_PAIRS=${XSTAMP_MAX_MOLECULE_PAIRS}
 )

set(exaStampMolecule_LINK_LIBRARIES exanbCompute exaStampCompute exanbParticleNeighbors exaStampParticleSpecies exp6 lennard_jones exaStampLJExp6RFParticleCharge exaStampReactionFieldParticleCharge exaStampPotential exanbIO exaStampIO)
onika_add_plugin(exaStampMolecule ${CMAKE_CURRENT_SOURCE_DIR})

