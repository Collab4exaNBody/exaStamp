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

#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
#include <onika/file_utils.h>

#include <exaStamp/io/read_stamp_v3.h>

#include <mpi.h>
#include <fstream>

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT>
  struct ReadStampV3Node : public OperatorNode
  {
    ADD_SLOT(MPI_Comm               , mpi              , INPUT , MPI_COMM_WORLD                       , DocString{"MPI communicator"} );
    ADD_SLOT(std::string            , filename         , INPUT , REQUIRED                             , DocString{"File name"} );
    ADD_SLOT(ReadBoundsSelectionMode, bounds_mode      , INPUT , ReadBoundsSelectionMode::FILE_BOUNDS , DocString{"Bounds mode"} );
    ADD_SLOT(double                 , enlarge_bounds   , INPUT , 0.0                                  , DocString{"Bounds enlargement"} );
    ADD_SLOT(bool                   , pbc_adjust_xform , INPUT , false                                , DocString{"Adjust xform to preserve sizes in directions with periodic boundary conditions"} );
    ADD_SLOT(Domain                 , domain           , INPUT_OUTPUT                                 , DocString{"Domain name"} );
    ADD_SLOT(GridT                  , grid             , INPUT_OUTPUT                                 , DocString{"Particle grid"} );
    ADD_SLOT(long                   , timestep         , OUTPUT                                       , DocString{"Iteration number"} );
    ADD_SLOT(double                 , physical_time    , INPUT_OUTPUT                                 , DocString{"Physical time, modified by reader"} );

    // -----------------------------------------------
    // ----------- Operator documentation ------------
    inline std::string documentation() const override final
    {
      return R"EOF(
        Reads a protection file in StampV3 format:
        this format is fully compatible with the Stamp3 MD code (read/write)
        )EOF";
    }

    inline void execute () override final
    {
      std::string file_name = onika::data_file_path( *filename );
      *physical_time = 0.0;
      read_stamp_v3(*mpi,file_name,*enlarge_bounds,*bounds_mode,*grid,*domain,*timestep,*physical_time, *pbc_adjust_xform);
    }

    inline void yaml_initialize(const YAML::Node& node) override final
    {
      YAML::Node tmp;
      if( node.IsScalar() )
      {
        tmp["filename"] = node;
      }
      else
      {
        tmp = node;
        if( tmp["file"] && ! tmp["filename"] )
        {
          lerr << "Warning, read_stamp_v3 now uses 'filename' for file to read, and no more 'file'. Please adjust your input file." << std::endl;
          tmp["file"] = node["filename"];
        }
      }
      this->OperatorNode::yaml_initialize( tmp );
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(read_stamp_v3)
  {
    OperatorNodeFactory::instance()->register_factory( "read_stamp_v3", make_grid_variant_operator< ReadStampV3Node > );
  }

}
