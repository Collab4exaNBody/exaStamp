#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <limits>

namespace exaStamp
{  

  // an atom can have up to 4 chemical bonds
  using AtomBondConnectivity = std::array<uint64_t,4>;
}

