#pragma once

#include <cstdint>
#include <exanb/core//basic_types_def.h>

namespace exaStamp
{
  using namespace exanb;
  struct SlipLinkField
  {
    Vec3d anchor_displ; // anchor relative position from SL position
    double xj_frac;
    uint64_t left_bead_id;
  };

}

