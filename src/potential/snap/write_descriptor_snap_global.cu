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
#include <onika/file_utils.h>
#include <onika/log.h>

#include <fstream>
#include <iomanip>
#include <string>
#include <mpi.h>

#include "../mlip-pod/include/npy_writer.h"

// Plain export of compute_descriptor_snap_global's output array. Unlike
// write_descriptor_pod_global (single-MPI-rank only), compute_descriptor_snap_global is multi-rank
// capable -- its output array is already identically MPI_Allreduce'd on every rank by the time this
// runs, so this writer just needs to pick one rank (0) to actually write, to avoid every rank racing
// to write the same file.
namespace exaStamp
{
  using namespace exanb;

  class WriteDescriptorSnapGlobal : public OperatorNode
  {
    ADD_SLOT( MPI_Comm , mpi , INPUT , REQUIRED );
    ADD_SLOT( onika::memory::CudaMMVector<double> , snap_global , INPUT , REQUIRED , DocString{"see compute_descriptor_snap_global"} );
    ADD_SLOT( long , ncoeff_all , INPUT , REQUIRED , DocString{"see compute_descriptor_snap_global"} );
    ADD_SLOT( std::string , format , INPUT , std::string("text") , DocString{"Output format: 'text' (default, one line per row, space-separated) or 'npy' (single .npy v1.0 file, shape (rows,ncoeff_all))."} );
    ADD_SLOT( std::string , filename , INPUT , std::string("snap_global.txt") , DocString{"Output file. In 'npy' format a trailing '.txt' is stripped and '.npy' appended."} );

  public:
    inline void execute() override final
    {
      int rank = 0;
      MPI_Comm_rank( *mpi, &rank );
      if( rank != 0 ) return;

      if( *format != "text" && *format != "npy" )
      {
        fatal_error() << "write_descriptor_snap_global: unknown format '"<<*format<<"' (choices: text, npy)" << std::endl;
      }

      const long nc = *ncoeff_all;
      const long rows = ( nc > 0 ) ? static_cast<long>(snap_global->size() / static_cast<size_t>(nc)) : 0;

      if( *format == "npy" )
      {
        std::string prefix = *filename;
        static constexpr const char * TXT_SUFFIX = ".txt";
        if( prefix.size() >= 4 && prefix.compare(prefix.size()-4, 4, TXT_SUFFIX) == 0 ) prefix.resize(prefix.size()-4);
        write_npy( onika::data_file_path(prefix+".npy"), {static_cast<size_t>(rows), static_cast<size_t>(nc)}, "<f8", snap_global->data(), sizeof(double) );
        return;
      }

      std::ofstream fout( onika::data_file_path(*filename) );
      fout << std::setprecision(17);
      for( long r=0; r<rows; r++ )
      {
        const double * const prow = snap_global->data() + static_cast<size_t>(r)*nc;
        for( long c=0; c<nc; c++ ) { if(c>0) fout << " "; fout << prow[c]; }
        fout << "\n";
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(

Writes compute_descriptor_snap_global's (1+3*natoms+6) x ncoeff array to a single file, from rank 0
only (the array is already identically MPI_Allreduce'd on every rank) -- row 0 = summed descriptor,
rows 1..3*natoms = per-atom gradient, rows 3*natoms+1..+6 = virial (Voigt order). See
compute_descriptor_snap_global's own documentation for the exact layout and how to use it for
linear-potential fitting.

'format: npy' writes a single .npy v1.0 file, shape (1+3*natoms+6, ncoeff), directly loadable with
numpy.load().

Usage example:

compute_descriptor_snap: { compute_derivative: true, parameters: { param: "W.snapparam", coef: "W.snapcoeff" } }
compute_descriptor_snap_global
write_descriptor_snap_global: { filename: "snap_global.txt" }

)EOF";
    }
  };

  ONIKA_AUTORUN_INIT(write_descriptor_snap_global)
  {
    OperatorNodeFactory::instance()->register_factory( "write_descriptor_snap_global", make_simple_operator< WriteDescriptorSnapGlobal > );
  }

}
