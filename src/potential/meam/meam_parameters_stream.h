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

#include "meam_parameters.h"

#include <ostream>

static inline std::ostream& operator << (std::ostream& out, const MeamParameters& p)
{
# define MEAM_PARAMETER(x) out<< #x << "=" << p.x << " "
  MEAM_PARAMETER(rmax);
  MEAM_PARAMETER(rmin);
  MEAM_PARAMETER(Ecoh);
  MEAM_PARAMETER(E0);
  MEAM_PARAMETER(A);
  MEAM_PARAMETER(r0);
  MEAM_PARAMETER(alpha);
  MEAM_PARAMETER(delta);
  MEAM_PARAMETER(beta0);
  MEAM_PARAMETER(beta1);
  MEAM_PARAMETER(beta2);
  MEAM_PARAMETER(beta3);
  MEAM_PARAMETER(t0);
  MEAM_PARAMETER(t1);
  MEAM_PARAMETER(t2);
  MEAM_PARAMETER(t3);
  MEAM_PARAMETER(s0);
  MEAM_PARAMETER(s1);
  MEAM_PARAMETER(s2);
  MEAM_PARAMETER(s3);
  MEAM_PARAMETER(Cmin);
  MEAM_PARAMETER(Cmax);
  MEAM_PARAMETER(Z);
  MEAM_PARAMETER(rc);
  MEAM_PARAMETER(rp);
# undef MEAM_PARAMETER
  return out;
}

