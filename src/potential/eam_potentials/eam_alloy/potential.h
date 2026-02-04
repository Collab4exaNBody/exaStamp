/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include "eam_alloy.h"

#define USTAMP_POTENTIAL_NAME eam_alloy
#define USTAMP_POTENTIAL_PARMS EamAlloyParameters

#define USTAMP_POTENTIAL_ENABLE_CUDA 1
#define USTAMP_POTENTIAL_EAM_MM 1 // set if potential supports multimaterial and thus needs type numbers as additional arguments

#define USTAMP_POTENTIAL_EAM_MM_INIT_TYPES(p,nt,pe) (p).initialize_types_table( nt , pe )
#define USTAMP_POTENTIAL_EAM_MM_FORCE eam_alloy_mm_force
#define USTAMP_POTENTIAL_EAM_EMB eam_alloy_fEmbed
#define USTAMP_POTENTIAL_EAM_EMB_DERIV eam_alloy_fEmbed_deriv
#define USTAMP_POTENTIAL_EAM_RHO_NODERIV eam_alloy_rho_noderiv

