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
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/histogram.h>

#include <mpi.h>
#include <fstream>
#include <string>

namespace exaStamp
{
  using namespace exanb;

  class WriteHistogram : public OperatorNode
  {
    ADD_SLOT( MPI_Comm    , mpi       , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( Histogram<> , histogram , INPUT , REQUIRED );
    ADD_SLOT( std::string , filename  , INPUT , std::string("histogram.csv") );
    ADD_SLOT( std::string , separator , INPUT , std::string(" ; ") );

  public:
    inline void execute () override final
    {
      int rank = 0;
      MPI_Comm_rank( *mpi, &rank );
      if( rank != 0 ) return;

      const Histogram<>& h = *histogram;
      size_t n = h.m_data.size();
      double bin_width = ( n > 0 ) ? ( h.m_max_val - h.m_min_val ) / n : 0.0;

      std::ofstream fout( *filename );
      fout << "# bin_center" << (*separator) << "count" << std::endl;
      for(size_t i=0;i<n;i++)
      {
        double center = h.m_min_val + ( i + 0.5 ) * bin_width;
        fout << center << (*separator) << h.m_data[i] << std::endl;
      }
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Writes a Histogram<> (as produced by histogram_energy, histogram_charge, histogram_vx/vy/vz,
histogram_rx/ry/rz, histogram_vnorm, histogram_fnorm, histogram_cell_particles, ...)
to a .csv file, one row per bin: bin center and count.

Usage example:
  - histogram_energy
  - write_histogram:
      filename: "energy_histogram.csv"
      separator: " "
)EOF";
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(write_histogram)
  {
    OperatorNodeFactory::instance()->register_factory( "write_histogram" , make_simple_operator< WriteHistogram > );
  }

}
