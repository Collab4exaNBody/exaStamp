#pragma once

#include "eam_buffer.h"
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

