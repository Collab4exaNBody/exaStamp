#pragma once

#include <exaStamp/npt/npt.h>
#include <onika/math/basic_types_stream.h>
#include <iostream>

namespace exaStamp
{

  inline std::ostream& operator << (std::ostream& out, const NPTConfig& data)
  {
    out<<"{ Tstart="<<data.m_Tstart<<", Tdamp="<<data.m_Tdamp<<" }";
    return out;
  }

}
