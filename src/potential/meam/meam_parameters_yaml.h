#include <yaml-cpp/yaml.h>
#include "meam_parameters.h"

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  template<> struct convert<MeamParameters>
  {
    static bool decode(const Node& node, MeamParameters& v);
  };
}

