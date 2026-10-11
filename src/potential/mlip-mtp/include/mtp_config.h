/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "../emtp.h"

struct MtpContext
{
  // Per-thread EMTP instances (index == OpenMP thread id).
  std::vector<std::shared_ptr<EMTP>> m_emtp;

  // Maps exaStamp 0-indexed particle type -> MTP 0-indexed species. PURELY POSITIONAL: unlike
  // POD, the .mtp file carries no species names at all (confirmed against both the real file
  // format and LAMMPS's own PairMTP::coeff(), which only accepts "pair_coeff * *" -- there is
  // nothing else to name-match against). type_map[i] = i by construction; mtp_init only
  // validates that species_count (file) matches the declared species: block's length and logs
  // the resulting mapping, since a misordered species: block cannot be caught structurally.
  std::vector<int> type_map;

  int nspecies = 1;
};
