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

#include <exaStamp/potential/eam/eam_buffer.h>
#include <yaml-cpp/yaml.h>

namespace YAML
{
  using exaStamp::EamMultimatParameters;
  template<class EamParametersT> struct convert< EamMultimatParameters<EamParametersT> >
  {
    static bool decode(const Node& node, EamMultimatParameters<EamParametersT> & v)
    {
      if( ! node.IsMap() ) return false;
      v.m_type_a = node["type_a"].as< std::string >();
      v.m_type_b = node["type_b"].as< std::string >();
      v.m_parameters = node["parameters"].as< EamParametersT >();
      return true;
    }
  };
}

