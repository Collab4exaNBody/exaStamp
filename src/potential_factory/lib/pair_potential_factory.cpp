#include <exaStamp/potential_factory/pair_potential_factory.h>

#include <exanb/core/grid.h>
#include <onika/log.h>
#include <exanb/core/plugin.h>
#include <onika/yaml/yaml_utils.h>

#include <cassert>
#include <iostream>
#include <memory>

namespace exaStamp
{
  using namespace exanb;

  using PairPotentialCreator = std::function<std::shared_ptr<PairPotential>(YAML::Node)>;

  bool PairPotentialFactory::register_factory( const std::string& name, PairPotentialCreator creator )
  {
    if( ! exanb::quiet_plugin_register() ) { lout<<"  potential   "<<name<<std::endl; }
    plugin_db_register( "potential" , name );
    s_factories[ name ] = creator;
    return true;
  }

  std::shared_ptr<PairPotential> PairPotentialFactory::make_instance( YAML::Node node )
  {
    /*
    ldbg << "PairPotentialFactory::make_instance, node content :" << std::endl;
    dump_node_to_stream( ldbg , node );
    ldbg << std::endl;
    */
    
    if( ! node["potential"] )
    {
      lerr << "Fatal: Yaml node does not have a 'potential' entry." << std::endl;
      std::abort();
    }
    
    std::string pot_name = node["potential"].as<std::string>();
    auto it = s_factories.find( pot_name );
    if( it == s_factories.end() )
    {
      std::string suggested_plugin = suggest_plugin_for( "potential" , pot_name );
      if( ! suggested_plugin.empty() )
      {
        ldbg << "auto loading "<< suggested_plugin<<" to find potential "<<pot_name<< std::endl;
        load_plugins( { suggested_plugin } );
        it = s_factories.find( pot_name );
      }
    }
        
    if( it == s_factories.end() )
    {
      lerr<<"Could not find a potential factory for name '"<<pot_name<<"'"<<std::endl;
      std::abort();
    }
    YAML::Node tmp;
    if( node["parameters"] )
    {
      tmp = node["parameters"];
    }
    else
    {
      tmp = "null";
    }
    return it->second ( tmp );
  }

  std::list<std::string> PairPotentialFactory::available_potentials()
  {
    std::list<std::string> result;
    for(const auto& f : s_factories)
    {
	    result.push_back( f.first );
    }
    return result;
  }

  std::map< std::string, PairPotentialCreator > PairPotentialFactory::s_factories;
}

