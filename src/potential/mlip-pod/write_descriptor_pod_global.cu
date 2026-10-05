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

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/memory/allocator.h>

#include <string>
#include <mpi.h>

#include "../mlip-utils/include/descriptor_writer_common.h"

// Plain export of compute_descriptor_pod_global's output array. compute_descriptor_pod_global is
// multi-rank capable and already identically MPI_Allreduce'd on every rank by the time this runs,
// so this writer just needs to pick one rank (0) to actually write, to avoid every rank racing to
// write the same file -- same pattern as write_descriptor_snap_global.cu.
namespace exaStamp
{
  using namespace exanb;

  class WriteDescriptorPodGlobal : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , pod_global , INPUT , REQUIRED , DocString{"see compute_descriptor_pod_global"} );
    ADD_SLOT( long , ncoeff_all , INPUT , REQUIRED , DocString{"see compute_descriptor_pod_global"} );
    ADD_SLOT( std::string , format , INPUT , std::string("text") , DocString{"Output format: 'text' (default, one line per row, space-separated) or 'npy' (single .npy v1.0 file, shape (rows,ncoeff_all))."} );
    ADD_SLOT( std::string , filename , INPUT , std::string("pod_global.txt") , DocString{"Output file. In 'npy' format a trailing '.txt' is stripped and '.npy' appended."} );

  public:
    inline void execute() override final
    {
      write_descriptor_global_common( "write_descriptor_pod_global", *mpi,
          pod_global->data(), pod_global->size(), *ncoeff_all, *format, *filename );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes compute_descriptor_pod_global's (1+3*natoms+6) x ncoeff_all array to a single file, from
rank 0 only (the array is already identically MPI_Allreduce'd on every rank) -- row 0 = global
per-configuration descriptor vector, rows 1..3*natoms = its gradient w.r.t. each atom id's x/y/z,
rows 3*natoms+1..+6 = virial (Voigt order) (see compute_descriptor_pod_global's documentation for
the exact layout and how to use it for linear-potential fitting).

'format: npy' writes a single .npy v1.0 file, shape (1+3*natoms+6, ncoeff_all), directly loadable
with numpy.load() -- the whole output is already one homogeneous 2-D array, no per-field splitting
needed.

Usage example:

pod_init: { parameters: { pod_file: "Ta_param.pod", coeff_file: "Ta_coefficients.pod" } }
compute_descriptor_pod_global
write_descriptor_pod_global: { filename: "pod_global.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(write_descriptor_pod_global)
  {
    OperatorNodeFactory::instance()->register_factory( "write_descriptor_pod_global", make_simple_operator< WriteDescriptorPodGlobal > );
  }

}
