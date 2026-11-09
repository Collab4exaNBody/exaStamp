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

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

// Lists the files with a given suffix under a directory tree (recursively), sorted by path, for
// next_database_file to iterate over (see create_descriptor_database_<family>.msp).
namespace exaStamp
{
  using namespace exanb;

  class ListFileDirectory : public OperatorNode
  {
    ADD_SLOT( std::string              , xyz_database  , INPUT , REQUIRED , DocString{"Directory to search recursively"} );
    ADD_SLOT( std::string              , pattern       , INPUT , std::string(".xyz") , DocString{"Filename suffix to match"} );
    ADD_SLOT( std::vector<std::string> , file_list     , OUTPUT , DocString{"Sorted list of matching file paths"} );
    ADD_SLOT( long                     , n_total_files , OUTPUT , DocString{"Number of matching files"} );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Lists the files whose name ends with `pattern` under the directory `xyz_database`, recursively
(one subdirectory per configuration category is typical), sorted by path. Aborts if no file
matches. Used with next_database_file to compute descriptors over a training database:

process_files_loop:
  loop: true
  condition: compute_desc_continue
  body:
    - next_database_file
    - grid_clear
    - read_xyz_file_with_xform      # reads `filename`
    - ...                           # init, compute_descriptor_<family>, compute_descriptor_<family>_global
    - write_descriptor_<family>_global: { format: npy }   # rebind filename to output_filename

global:
  xyz_database: "Ta_DB_XYZ"
  desc_database: "Ta_Desc_SNAP"
  compute_desc_continue: true

simulation body: ... init_prolog, list_file_directory, process_files_loop ...

See data/regression_new/mlip_training/create_descriptor_database_*.msp for complete examples.
)EOF";
    }

    inline void execute() override final
    {
      namespace fs = std::filesystem;
      std::vector<std::string> files;
      for( const auto & entry : fs::recursive_directory_iterator( *xyz_database ) )
      {
        if( ! entry.is_regular_file() ) continue;
        const std::string name = entry.path().filename().string();
        const std::string & suf = *pattern;
        if( name.size() >= suf.size() && name.compare( name.size()-suf.size(), suf.size(), suf ) == 0 )
        {
          files.push_back( entry.path().string() );
        }
      }
      std::sort( files.begin(), files.end() );
      *n_total_files = static_cast<long>( files.size() );
      *file_list = std::move( files );
      lout << "list_file_directory: found "<< *n_total_files <<" file(s) matching '"<< *pattern <<"' under "<< *xyz_database << std::endl;
      if( *n_total_files == 0 )
      {
        fatal_error() << "list_file_directory: no '"<<*pattern<<"' files found under '"<<*xyz_database<<"' -- nothing to process" << std::endl;
      }
    }
  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(list_file_directory)
  {
    OperatorNodeFactory::instance()->register_factory( "list_file_directory", make_simple_operator< ListFileDirectory > );
  }

}
