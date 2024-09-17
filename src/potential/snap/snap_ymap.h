#pragma once

namespace exaStamp
{
# include "ymap_jmax2.hxx"
# include "ymap_jmax3.hxx"
# include "ymap_jmax4.hxx"

  static inline constexpr int snap_force_Y_count(int twojmax, int idxu_max)
  {
    const int jmax = twojmax / 2;
    if( jmax == 2 ) return Y_jmax2_jju_count;
    else if( jmax == 3 ) return Y_jmax3_jju_count;
    else if( jmax == 4 ) return Y_jmax4_jju_count;
    else return idxu_max;
  }
  
  static inline constexpr bool snap_force_use_Y(int twojmax , int jju)
  {
    const int jmax = twojmax / 2;
    if( jmax == 2 ) return Y_jmax2_jju_map[jju] != -1;
    else if( jmax == 3 ) return Y_jmax3_jju_map[jju] != -1;
    else if( jmax == 4 ) return Y_jmax4_jju_map[jju] != -1;
    else return true;
  }

  static inline constexpr int snap_force_Y_map(int twojmax , int jju)
  {
    const int jmax = twojmax / 2;
    if( jmax == 2 ) return Y_jmax2_jju_map[jju];
    else if( jmax == 3 ) return Y_jmax3_jju_map[jju];
    else if( jmax == 4 ) return Y_jmax4_jju_map[jju];
    else return jju;
  }

}

