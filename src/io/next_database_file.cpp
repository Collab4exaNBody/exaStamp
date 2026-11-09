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
#include <onika/log.h>

#include <filesystem>
#include <string>
#include <vector>

// Loop cursor over list_file_directory's file list (see create_descriptor_database_<family>.msp).
namespace exaStamp
{
  using namespace exanb;

  class NextDatabaseFile : public OperatorNode
  {
    ADD_SLOT( std::vector<std::string> , file_list           , INPUT , REQUIRED , DocString{"File list from list_file_directory"} );
    ADD_SLOT( long                     , n_total_files       , INPUT , REQUIRED , DocString{"Number of files from list_file_directory"} );
    ADD_SLOT( std::string              , desc_database        , INPUT , REQUIRED , DocString{"Output directory (created if needed); output_filename = desc_database/<source file stem>, without extension"} );
    ADD_SLOT( long                     , cursor               , INPUT_OUTPUT , 0 , DocString{"Index of the next file"} );
    ADD_SLOT( std::string              , filename             , OUTPUT , DocString{"Next file to read"} );
    ADD_SLOT( std::string              , output_filename      , OUTPUT , DocString{"Output path prefix for that file"} );
    // must be INPUT_OUTPUT to be seen by the loop condition, and must be seeded in `global:`
    // (compute_desc_continue: true) since the condition is evaluated before the first iteration
    ADD_SLOT( bool                     , compute_desc_continue, INPUT_OUTPUT );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Loop cursor over the file list of list_file_directory. At each call, outputs the next file as
`filename` and an output prefix `output_filename` = desc_database/<file stem>, and sets
compute_desc_continue to false once every file has been processed. compute_desc_continue must be
set to true in `global:`. See list_file_directory for an example.
)EOF";
    }

    inline void execute() override final
    {
      if( *cursor < *n_total_files )
      {
        const std::string & src = (*file_list)[*cursor];
        *filename = src;
        // the descriptor writers do not create missing directories
        std::filesystem::create_directories( *desc_database );
        // no extension: the npy writer appends ".npy" itself
        *output_filename = *desc_database + "/" + std::filesystem::path(src).stem().string();
        ++(*cursor);
        *compute_desc_continue = true;
      }
      else
      {
        *compute_desc_continue = false;
      }
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(next_database_file)
  {
    OperatorNodeFactory::instance()->register_factory( "next_database_file", make_simple_operator< NextDatabaseFile > );
  }

}
