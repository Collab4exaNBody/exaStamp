#pragma once

#include <exaStamp/parrinellorahman/parrinellorahman.h>
#include <exanb/core/basic_types_stream.h>
#include <iostream>

namespace exaStamp
{

  inline std::ostream& operator << (std::ostream& out, const ParrinelloRahmanConfig& data)
  {
    out<<"{ Text="<<data.m_Text<<", masseNVT="<<data.m_masseNVT<<" , Pext="<<data.m_Pext<<", m_masseB="<<data.m_masseB<<" }";
    return out;
  }

}
