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

#include <exaStamp/mlip_ye3t/ye3t_tables.h>
#include <exaStamp/mlip_utils/mlip_flat_bonds.h>
#include <exaStamp/yaml_check_keys.h>

#include <yaml-cpp/yaml.h>

#include <memory>
#include <string>
#include <vector>

namespace exaStamp
{
  struct Ye3tParams
  {
    std::string file;
    bool analytic = false;      // radial functions evaluated directly instead of the splines (debug, not the reference)
    double mem_gb = 1.0;        // budget of the per-chunk buffers
    int emulate = 0;            // debug: > 0 runs the GPU code paths on the host, blocks of that many threads
  };

  // model built by ye3t_init, with the flat arrays the force operator gathers into
  struct Ye3tContext
  {
    std::shared_ptr<mlip_ye3t::Ye3tTables> tables;
    std::vector<int> species_type;    // exaStamp species index -> model element (-1: not in the model)
    double mem_gb = 1.0;
    int emulate = 0;
    mlip_utils::FlatBonds flat;
    mlip_ye3t::Ye3tBuffers buf;
  };
}

namespace YAML
{
  template<> struct convert<exaStamp::Ye3tParams>
  {
    static bool decode(const Node& node, exaStamp::Ye3tParams& v)
    {
      if( !node.IsMap() ) { return false; }
      exaStamp::check_yaml_keys( node , { "file", "analytic", "mem_gb", "emulate" } , "ye3t parameters" );
      if( !node["file"] ) { return false; }
      v.file = node["file"].as<std::string>();
      if( node["analytic"] ) v.analytic = node["analytic"].as<bool>();
      if( node["mem_gb"] ) v.mem_gb = node["mem_gb"].as<double>();
      if( node["emulate"] ) v.emulate = node["emulate"].as<int>();
      return true;
    }
  };
}
