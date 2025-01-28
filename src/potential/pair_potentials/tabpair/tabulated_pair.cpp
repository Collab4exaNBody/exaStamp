#include "tabulated_pair.h"
#include <exanb/core/file_utils.h>
#include <onika/yaml/yaml_utils.h>
#include <onika/log.h>
#include <onika/physics/constants.h>

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{

  bool convert<exaStamp::TabPairPotentialParms>::decode(const Node& in_node, exaStamp::TabPairPotentialParms& v)
  {
    using exanb::UnityConverterHelper;
    using onika::physics::Quantity;
    using exanb::lout;
    using exanb::lerr;
    using exanb::fatal_error;
    using exanb::yaml_load_file_abort_on_except;
    auto& ldbg = exanb::lout; //using exanb::ldbg;
    
    std::string file_to_load;
    
    if( in_node.IsScalar() )
    {
      file_to_load = in_node.as<std::string>();
    }
    else
    {
      if( ! in_node.IsMap() )
      {
        lerr << "Decoding exaStamp::TabPairPotentialParms : YAML node is not a map as expected"<<std::endl;
        return false;
      }
      if( in_node["file"] )
      {
        file_to_load = in_node["file"].as<std::string>();
      }
    }

    Node node;
    if( ! file_to_load.empty() )
    {
      file_to_load = exanb::data_file_path(file_to_load);
      ldbg << "tabulated pair data from "<<file_to_load<<std::endl;
      node = yaml_load_file_abort_on_except(file_to_load);
    }
    else
    {
      ldbg << "interpreting node as inline raw data"<<std::endl;
      node = in_node;
    }

    if( ! node.IsMap() )
    {
      lerr << "Decoding exaStamp::TabPairPotentialParms : YAML data node is not a map as expected"<<std::endl;
      return false;
    }


    double ds = 1.0; // derivative scale factor
    double es = 1.0; // energy scale factor
    double ls = 1.0; // length scale factor
    if( node["format"] )
    {
      if (node["format"].as<std::string>() == "exastampv1")
      {
        ds = onika::physics::avogadro * 1.e-11 ;
        ls = exanb::make_quantity(1.0,"m").convert();
        es = exanb::make_quantity(1.0,"J").convert();
        ldbg << "exastampv1 format : ds=" <<ds<<", ls="<<ls<<", es="<<es<< std::endl;
      }
    }

    if( ! node["r"] )
    {
      lerr << "missing 'r' key in tabulated file "<<file_to_load<<std::endl;
      return false;
    }
    if( ! node["e"] )
    {
      lerr << "missing 'e' key in tabulated file "<<file_to_load<<std::endl;
      return false;
    }
    if( ! node["de"] )
    {
      lerr << "missing 'de' key in tabulated file "<<file_to_load<<std::endl;
      return false;
    }

    std::vector<double> support;
    std::vector<double> value;
    std::vector<double> derivative;    
    for(const Quantity& q : node["r"].as< std::vector<Quantity> >() ) { support.push_back( ls * q.convert() ); }
    for(const Quantity& q : node["e"].as< std::vector<Quantity> >() ) { value.push_back( es * q.convert() ); }
    for(const Quantity& q : node["de"].as< std::vector<Quantity> >() ) { derivative.push_back( ds * q.convert() ); }
    if( support.size() != value.size() || support.size() != derivative.size() || support.empty() )
    {
      lerr << "bad data length: r as "<<support.size()<<" values, e as "<<value.size()<<" values, de as "<<derivative.size()<<" values"<<std::endl;
      return false;
    }
    ldbg << "tabulated pair read "<<support.size()<<" values"<<std::endl;
#   ifndef NDEBUG
    size_t n = support.size();
    for(size_t i=0;i<n;i++)
    {
      ldbg << support[i]<< " ; " << value[i] << " ; " << derivative[i] << std::endl;
    }
#   endif
    v.m_e.set_points( support, value );
    v.m_de.set_points( support, derivative );
    return true;
  }

}

