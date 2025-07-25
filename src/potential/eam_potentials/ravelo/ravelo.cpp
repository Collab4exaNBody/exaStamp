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

#include "ravelo.h"
#include <onika/file_utils.h>
#include <onika/log.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{

  bool convert<exaStamp::EamRaveloParameters>::decode(const Node& _node, exaStamp::EamRaveloParameters& v)
  {
    using onika::physics::Quantity;
    using exanb::lout;
    using exanb::lerr;
    using exanb::ldbg;
    
    Node node = _node;
    std::string file_to_load;
    
    if( node.IsScalar() )
    {
      file_to_load = node.as<std::string>();
    }
    else
    {
      if( !node.IsMap() ) { return false; }
      if( node["file"] )
      {
        file_to_load = node["file"].as<std::string>();
      }
    }

    if( ! file_to_load.empty() )
    {
      file_to_load = onika::data_file_path(file_to_load);
      ldbg << "parameters and tabulated EAM data from "<<file_to_load<<std::endl;
      node = LoadFile(file_to_load);
    }

    double fs = 1.0; // length scale factor
    double dfs = 1.0; // length scale factor
    if( node["format"] )
    {
      if (node["format"].as<std::string>() == "exastampv1")
      {
        fs = onika::physics::avogadro * 1.e-1 ; // .... <--
        dfs = onika::physics::avogadro * 1.e-3 ;
      }
    }

    if( !node.IsMap() ) { return false; }
    std::vector<double> support;
    std::vector<double> value;
    std::vector<double> derivative;

    for(const Quantity& q : node["rho"].as< std::vector<Quantity> >() ) { support.push_back(q.convert()); }
    for(const Quantity& q : node["f"].as< std::vector<Quantity> >() ) { value.push_back(q.convert()*fs); }
    for(const Quantity& q : node["df"].as< std::vector<Quantity> >() ) { derivative.push_back(q.convert()*dfs); }

    if( support.size() != value.size() || support.size() != derivative.size() || support.empty() )
    {
      lerr << "bad data length: rhof as "<<support.size()<<" values, f as "<<value.size()<<" values, df as "<<derivative.size()<<" values"<<std::endl;
      return false;
    }
    
    v.m_f.set_points( support, value );
    v.m_df.set_points( support, derivative );

    support.clear();
    value.clear();
    derivative.clear();

    //    if( !node.IsMap() ) { return false; }
#     define CONVERT_EAMRAVELO_PARAM(x) v.x = node[#x].as<Quantity>().convert() //; std::cout << #x << " = " << v.x << std::endl
      CONVERT_EAMRAVELO_PARAM(Ec);
      CONVERT_EAMRAVELO_PARAM(alpha);
      CONVERT_EAMRAVELO_PARAM(a0);
      CONVERT_EAMRAVELO_PARAM(f3);
      CONVERT_EAMRAVELO_PARAM(f4);
      CONVERT_EAMRAVELO_PARAM(U0);
      CONVERT_EAMRAVELO_PARAM(r1);
      CONVERT_EAMRAVELO_PARAM(alphap);
      CONVERT_EAMRAVELO_PARAM(beta3);
      CONVERT_EAMRAVELO_PARAM(beta4);
      CONVERT_EAMRAVELO_PARAM(rs);
      CONVERT_EAMRAVELO_PARAM(s);
      CONVERT_EAMRAVELO_PARAM(rc);
      CONVERT_EAMRAVELO_PARAM(a1);
      CONVERT_EAMRAVELO_PARAM(a2);
      CONVERT_EAMRAVELO_PARAM(a3);
      CONVERT_EAMRAVELO_PARAM(a4);
      CONVERT_EAMRAVELO_PARAM(p);
      CONVERT_EAMRAVELO_PARAM(q);
      CONVERT_EAMRAVELO_PARAM(rho0);
#     undef CONVERT_EAMRAVELO_PARAM

    return true;
  }

}

