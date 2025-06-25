#pragma once

#include <onika/math/basic_types.h>
#include <vector>
#include <string>
#include <memory>

namespace exaStamp
{

  struct ScalarSourceTerm
  {
    virtual inline double S( exanb::Vec3d r, double t ) const { return 0; }
    virtual inline ~ScalarSourceTerm() = default;
  };

  std::shared_ptr<ScalarSourceTerm> make_source_term( const std::string& stype, const std::vector<double>& p );
}

