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

// Cursor operator for process_files_loop's batch{loop:true} body (see
// create_descriptor_database.msp): each iteration emits the next entry of list_file_directory's
// file_list (as `filename`, auto-wired into read_xyz_file_with_xform), a matching per-file output
// path re-using that file's own stem (as `output_filename`, auto-wired into write_descriptor_snap
// -- the input file's own name IS the label, no separate file_id column/sidecar manifest needed),
// and the loop-continue boolean `compute_desc_continue` the batch's condition: watches.
namespace exaStamp
{
  using namespace exanb;

  class NextDatabaseFile : public OperatorNode
  {
    ADD_SLOT( std::vector<std::string> , file_list           , INPUT , REQUIRED );
    ADD_SLOT( long                     , n_total_files       , INPUT , REQUIRED );
    ADD_SLOT( std::string              , desc_database        , INPUT , REQUIRED , DocString{"Output directory; output_filename = desc_database/<source file stem> (extension-less prefix, matching write_descriptor_*'s own npy-prefix convention)"} );
    ADD_SLOT( long                     , cursor               , INPUT_OUTPUT , 0 );
    ADD_SLOT( std::string              , filename             , OUTPUT );
    ADD_SLOT( std::string              , output_filename      , OUTPUT );
    // INPUT_OUTPUT, no literal default here (matches onika's sim_continue's own loop-condition
    // slot `result`) -- a batch{loop:true} condition slot is fed via the framework's loop_output/
    // loop_input propagation, which only wires up correctly against an INPUT_OUTPUT slot; a pure
    // OUTPUT slot's write is invisible to eval_condition() (confirmed by an infinite-loop repro),
    // and adding a literal default here (instead of seeding via `global:`) caused a resource-cycle
    // crash on startup -- seed the initial value via `global: { compute_desc_continue: true }`
    // in the .msp instead (needed since eval_condition() runs once before the body's first pass).
    ADD_SLOT( bool                     , compute_desc_continue, INPUT_OUTPUT );

  public:
    inline void execute() override final
    {
      if( *cursor < *n_total_files )
      {
        const std::string & src = (*file_list)[*cursor];
        *filename = src;
        // std::ofstream on a missing parent dir silently no-ops (write_descriptor_snap does not
        // check is_open()) -- create it here so a missing desc_database doesn't silently drop output
        std::filesystem::create_directories( *desc_database );
        // no extension here: write_descriptor_snap's npy path uses `filename` as a bare prefix and
        // appends ".npy" itself (it only strips a trailing ".txt", so passing "*.npy" here would
        // double up into "*.npy.npy")
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
