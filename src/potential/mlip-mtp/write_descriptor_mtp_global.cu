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

// Writes compute_descriptor_mtp_global's array from rank 0 (the array is identical on every rank).
namespace exaStamp
{
  using namespace exanb;

  class WriteDescriptorMtpGlobal : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , mtp_global , INPUT , REQUIRED , DocString{"see compute_descriptor_mtp_global"} );
    ADD_SLOT( long , ncoeff_all , INPUT , REQUIRED , DocString{"see compute_descriptor_mtp_global"} );
    ADD_SLOT( std::string , format , INPUT , std::string("text") , DocString{"Output format: 'text' (default, one line per row, space-separated) or 'npy' (single .npy v1.0 file, shape (rows,ncoeff_all))."} );
    ADD_SLOT( std::string , filename , INPUT , std::string("mtp_global.txt") , DocString{"Output file. In 'npy' format a trailing '.txt' is stripped and '.npy' appended."} );

  public:
    inline void execute() override final
    {
      write_descriptor_global_common( "write_descriptor_mtp_global", *mpi,
          mtp_global->data(), mtp_global->size(), *ncoeff_all, *format, *filename );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes the (1+3*natoms+6) x (species_count+ncoeff) array of compute_descriptor_mtp_global to a single file,
from rank 0: row 0 = summed descriptor, rows 1..3*natoms = force-signed gradient, last 6 rows =
virial (Voigt order). See compute_descriptor_mtp_global for the layout.

'format: npy' writes a .npy v1.0 file of shape (1+3*natoms+6, species_count+ncoeff), loadable with numpy.load().

Usage example:

init_parameters:
  - species
  - mtp_init: { parameters: { mtp_file: "pot.almtp" } }

compute_descriptor_mtp: { compute_derivative: true }
compute_descriptor_mtp_global
write_descriptor_mtp_global: { filename: "mtp_global.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(write_descriptor_mtp_global)
  {
    OperatorNodeFactory::instance()->register_factory( "write_descriptor_mtp_global", make_simple_operator< WriteDescriptorMtpGlobal > );
  }

}
