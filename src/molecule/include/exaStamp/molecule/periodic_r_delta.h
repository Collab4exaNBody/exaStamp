#pragma once

#include <exanb/core/basic_types.h>

namespace exaStamp
{

  static inline Vec3d periodic_r_delta(const Vec3d& r1, const Vec3d& r2, const Vec3d& size_box, double half_min_size_box)
  {
    Vec3d rij = r2 - r1;

    if(rij.x>  half_min_size_box) rij.x -= size_box.x;
    if(rij.x< -half_min_size_box) rij.x += size_box.x;
    if(rij.y>  half_min_size_box) rij.y -= size_box.y;
    if(rij.y< -half_min_size_box) rij.y += size_box.y;
    if(rij.z>  half_min_size_box) rij.z -= size_box.z;
    if(rij.z< -half_min_size_box) rij.z += size_box.z;
    assert(norm(rij)<half_min_size_box);
    assert( (rij != Vec3d{0,0,0}) );
    
    return rij;
  }

}
