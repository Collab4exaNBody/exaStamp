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

// rigid molecule
#define USTAMP_POTENTIAL_KIND USTAMP_POTENTIAL_KIND_RIGIDMOL
#define USTAMP_POTENTIAL_RIGIDMOL 1
#define OPERATOR_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_rigidmol_force)
#include "rigidmol_force_op_singlemat.h"
#define POTENTIAL_ADDITIONAL_FIELDS ,field::_orient,field::_couple

