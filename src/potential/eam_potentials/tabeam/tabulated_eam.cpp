#include "tabulated_eam.h"
#include <exanb/core/file_utils.h>
#include <onika/log.h>
#include <onika/physics/constants.h>

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{

  bool convert<exaStamp::TabEAMPotentialParms>::decode(const Node& _node, exaStamp::TabEAMPotentialParms& v)
  {
    using exanb::UnityConverterHelper;
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
      file_to_load = exanb::data_file_path(file_to_load);
      ldbg << "tabulated EAM data from "<<file_to_load<<std::endl;
      node = LoadFile(file_to_load);
    }

    double dphis = 1.0;
    double phis = 1.0;
    double drhos = 1.0; // drho scale factor
    double ls = 1.0; // length scale factor
    double fs = 1.0; // length scale factor
    double dfs = 1.0; // length scale factor
    double rhos = 1.0;
    if( node["format"] )
    {
      if (node["format"].as<std::string>() == "exastampv1")
      {
        phis = onika::physics::avogadro * 1.e-1 ;
        dphis = onika::physics::avogadro * 1.e-11;
        fs = onika::physics::avogadro * 1.e-1 ; // .... <--
        dfs = onika::physics::avogadro * 1.e-3 ;
        drhos = 1.e-8; // drho scale factor
        rhos = 1.0;
        ls = exanb::make_quantity(1.0,"m").convert();
      }
    }

    if( !node.IsMap() ) { return false; }
    std::vector<double> support;
    std::vector<double> value;
    std::vector<double> derivative;
    for(const Quantity& q : node["r"].as< std::vector<Quantity> >() ) { support.push_back(q.convert()*ls); }
    for(const Quantity& q : node["phi"].as< std::vector<Quantity> >() ) { value.push_back(q.convert()*phis); }
    for(const Quantity& q : node["dphi"].as< std::vector<Quantity> >() ) { derivative.push_back(q.convert()*dphis); }
    v.m_phi.set_points( support, value );
    v.m_dphi.set_points( support, derivative );

    if( support.size() != value.size() || support.size() != derivative.size() || support.empty() )
    {
      lerr << "bad data length: r as "<<support.size()<<" values, phi as "<<value.size()<<" values, dphi as "<<derivative.size()<<" values"<<std::endl;
      return false;
    }

//      support.clear();
    value.clear();
    derivative.clear();
//      for(const Quantity& q : node["r"].as< std::vector<Quantity> >() ) { support.push_back(q.convert()); }
    for(const Quantity& q : node["rho"].as< std::vector<Quantity> >() ) { value.push_back(q.convert()*rhos); }
    for(const Quantity& q : node["drho"].as< std::vector<Quantity> >() ) { derivative.push_back(q.convert()*drhos); }
    v.m_rho.set_points( support, value );
    v.m_drho.set_points( support, derivative );

    if( support.size() != value.size() || support.size() != derivative.size() || support.empty() )
    {
      lerr << "bad data length: r as "<<support.size()<<" values, rho as "<<value.size()<<" values, drho as "<<derivative.size()<<" values"<<std::endl;
      return false;
    }

    support.clear();
    value.clear();
    derivative.clear();
    for(const Quantity& q : node["rhof"].as< std::vector<Quantity> >() ) { support.push_back(q.convert()); }
    for(const Quantity& q : node["f"].as< std::vector<Quantity> >() ) { value.push_back(q.convert()*fs); }
    for(const Quantity& q : node["df"].as< std::vector<Quantity> >() ) { derivative.push_back(q.convert()*dfs); }
    v.m_f.set_points( support, value );
    v.m_df.set_points( support, derivative );

    if( support.size() != value.size() || support.size() != derivative.size() || support.empty() )
    {
      lerr << "bad data length: rhof as "<<support.size()<<" values, f as "<<value.size()<<" values, df as "<<derivative.size()<<" values"<<std::endl;
      return false;
    }

    return true;
  }

}

