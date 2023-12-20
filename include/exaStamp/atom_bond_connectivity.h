#pragma once

#include <array>
#include <cstdint>

namespace exaStamp
{  

  // an atom can have up to 4 chemical bonds
  using AtomBondConnectivity = std::array<uint64_t,4>;

}

