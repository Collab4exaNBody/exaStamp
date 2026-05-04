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

#include "eapod.h"

struct PodContext
{
  // Potential file paths — kept for lazy per-thread EAPOD construction.
  std::string pod_file;
  std::string coeff_file;

  // Per-thread EAPOD instances (index == OpenMP thread id).
  // Each owns its own working memory (tmpmem, bd, bdd, …).
  std::vector<std::shared_ptr<EAPOD>> m_eapod;

  // Maps exaStamp 0-indexed particle type → POD 1-indexed element type.
  std::vector<int> type_map;

  int nspecies = 1;
};
