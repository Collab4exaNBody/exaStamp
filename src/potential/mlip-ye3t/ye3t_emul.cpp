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

// Debug build of the YE3T kernels with the GPU semantics, run on the host (mlip_onika_emul.h): kernel namespace and
// driver renamed, no clash with the regular build of the same headers.

#include <exaStamp/mlip_utils/mlip_onika_emul.h>

#define ye3t_onika ye3t_onika_emulated
#define ye3t_onika_compute ye3t_onika_compute_emulated
#include <exaStamp/mlip_ye3t/ye3t_driver.h>

namespace exaStamp
{
  namespace mlip_ye3t
  {
    void ye3t_onika_emul( const Ye3tTables& T, Ye3tBuffers& buf, const mlip_utils::BondsView& bv, mlip_utils::OutputsView& out,
                          int nthreads, double mem_gb, int& n_rho_limit )
    {
      if( nthreads > ONIKA_CU_MAX_THREADS_PER_BLOCK ) nthreads = ONIKA_CU_MAX_THREADS_PER_BLOCK;
      ye3t_onika_compute_emulated( T, buf, bv, out, true, mem_gb, ::mlip_onika_emul::EmulLauncher{ nthreads }, n_rho_limit );
    }
  }
}
