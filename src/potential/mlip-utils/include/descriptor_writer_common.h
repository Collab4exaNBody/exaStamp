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

#include "npy_writer.h"
#include <onika/file_utils.h>
#include <onika/log.h>

#include <fstream>
#include <iomanip>
#include <string>
#include <mpi.h>

// Shared _global writer body for every write_descriptor_<family>_global operator (POD/SNAP/k2b/MTP)
// -- these were ~90% byte-identical hand-rolled copies before this header existed (confirmed by
// diff), differing only by slot/type names. Header-only, no link dependency, no CMake registration
// needed for this directory (headers-only, found via relative #include, not add_subdirectory).
// Each per-family write_descriptor_<family>_global.cu keeps its own thin OperatorNode (ADD_SLOT
// types genuinely differ per family), and its execute() just calls this one function.
//
// The per-atom counterpart (write_descriptor_common) that used to live alongside this one was
// removed entirely -- write_descriptor_pod.cu/write_descriptor_snap.cu/write_descriptor_mtp.cu/
// write_descriptor_k2b.cu are all gone now, only write_descriptor_<family>_global is built for any
// family. Per-atom classification use cases consume compute_descriptor_<family>'s output directly
// in-graph instead (see compute_slcsa.msp/compute_bispectrum.msp), no writer involved.
namespace exaStamp
{
  using namespace exanb;

  // compute_descriptor_<family>_global's (1+3*natoms+6) x ncoeff_all array is already identically
  // MPI_Allreduce'd on every rank by the time this runs, so this just picks one rank (0) to
  // actually write, no gather needed.
  inline void write_descriptor_global_common(
      const char * op_name,
      MPI_Comm mpi_comm,
      const double * arr,
      size_t arr_size,
      long ncoeff_all,
      const std::string & format,
      const std::string & filename )
  {
    int rank = 0;
    MPI_Comm_rank( mpi_comm, &rank );
    if( rank != 0 ) return;

    if( format != "text" && format != "npy" )
    {
      fatal_error() << op_name << ": unknown format '"<<format<<"' (choices: text, npy)" << std::endl;
    }

    const long nc = ncoeff_all;
    const long rows = ( nc > 0 ) ? static_cast<long>(arr_size / static_cast<size_t>(nc)) : 0;

    if( format == "npy" )
    {
      std::string prefix = filename;
      static constexpr const char * TXT_SUFFIX = ".txt";
      if( prefix.size() >= 4 && prefix.compare(prefix.size()-4, 4, TXT_SUFFIX) == 0 ) prefix.resize(prefix.size()-4);
      write_npy( onika::data_file_path(prefix+".npy"), {static_cast<size_t>(rows), static_cast<size_t>(nc)}, "<f8", arr, sizeof(double) );
      return;
    }

    std::ofstream fout( onika::data_file_path(filename) );
    fout << std::setprecision(17);
    for( long r=0; r<rows; r++ )
    {
      const double * const prow = arr + static_cast<size_t>(r)*nc;
      for( long c=0; c<nc; c++ ) { if(c>0) fout << " "; fout << prow[c]; }
      fout << "\n";
    }
  }
}
