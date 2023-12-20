#pragma once
#include <onika/cuda/cuda.h>
#include "snap_ext.h"
#include <cstdlib>

namespace SnapExt
{
  template<int JMAX> struct SnapConstants;

  template<int N> struct ConstOneIndices{ int idx[N]; };

  template<> struct SnapConstants<3>
  {
    static constexpr int compact_gsh_size = 140;
    static constexpr int gsh_start = 1;
    static constexpr int m_jmax = 3;
    static constexpr int m_two_jmax = 6;
    static constexpr int gsh_data_size = 343;
    static constexpr int n_bs = 140;
    static constexpr size_t main_data_size = 140;
    static constexpr size_t cgprod_data_size = 5361;
    static constexpr size_t subit_data_size = 5361;
    static constexpr size_t n_cmm_ones = 28;
    static constexpr ConstOneIndices<n_cmm_ones> cmm_init_ones = { { 0,1,4,5,9,13,14,19,24,29,30,36,42,48,54,55,62,69,76,83,90,91,99,107,115,123,131,139 } };

    static constexpr double gsh_sqrt_frac[7][7] = {
     { 0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0 }
    ,{ 0x0p+0,0x1p+0,0x1.6a09e667f3bcdp-1,0x1.279a74590331cp-1,0x1p-1,0x1.c9f25c5bfedd9p-2,0x1.a20bd700c2c3ep-2 }
    ,{ 0x0p+0,0x1.6a09e667f3bcdp+0,0x1p+0,0x1.a20bd700c2c3ep-1,0x1.6a09e667f3bcdp-1,0x1.43d136248490fp-1,0x1.279a74590331cp-1 }
    ,{ 0x0p+0,0x1.bb67ae8584caap+0,0x1.3988e1409212ep+0,0x1p+0,0x1.bb67ae8584caap-1,0x1.8c97ef43f7248p-1,0x1.6a09e667f3bcdp-1 }
    ,{ 0x0p+0,0x1p+1,0x1.6a09e667f3bcdp+0,0x1.279a74590331cp+0,0x1p+0,0x1.c9f25c5bfedd9p-1,0x1.a20bd700c2c3ep-1 }
    ,{ 0x0p+0,0x1.1e3779b97f4a8p+1,0x1.94c583ada5b53p+0,0x1.4a7e9cb8a3491p+0,0x1.1e3779b97f4a8p+0,0x1p+0,0x1.d363d1848dcbfp-1 }
    ,{ 0x0p+0,0x1.3988e1409212ep+1,0x1.bb67ae8584caap+0,0x1.6a09e667f3bcdp+0,0x1.3988e1409212ep+0,0x1.186f174f88472p+0,0x1p+0 }
    };
  };

template<> struct SnapConstants<4>
{
  static constexpr int compact_gsh_size = 285;
  static constexpr int gsh_start = 1;
  static constexpr int m_jmax = 4;
  static constexpr int m_two_jmax = 8;
  static constexpr int gsh_data_size = 729;
  static constexpr int n_bs = 285;
  static constexpr size_t main_data_size = 285;
  static constexpr size_t cgprod_data_size = 27100;
  static constexpr size_t subit_data_size = 27100;
  static constexpr size_t n_cmm_ones = 45;
  static constexpr ConstOneIndices<n_cmm_ones> cmm_init_ones = 
  { { 0,1,4,5,9,13,14,19,24,29,30,36,42,48,54,55,62,69,76,83,90,91,99,107,115,123,131,139,140,149,158,167,176,185,194,203,204,214,224,234,244,254,264,274,284 } };

  static constexpr double gsh_sqrt_frac[9][9] = {
   { 0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0,0x0p+0 }
  ,{ 0x0p+0,0x1p+0,0x1.6a09e667f3bcdp-1,0x1.279a74590331cp-1,0x1p-1,0x1.c9f25c5bfedd9p-2,0x1.a20bd700c2c3ep-2,0x1.83091e6a7f7e6p-2,0x1.6a09e667f3bcdp-2 }
  ,{ 0x0p+0,0x1.6a09e667f3bcdp+0,0x1p+0,0x1.a20bd700c2c3ep-1,0x1.6a09e667f3bcdp-1,0x1.43d136248490fp-1,0x1.279a74590331cp-1,0x1.11acee560242ap-1,0x1p-1 }
  ,{ 0x0p+0,0x1.bb67ae8584caap+0,0x1.3988e1409212ep+0,0x1p+0,0x1.bb67ae8584caap-1,0x1.8c97ef43f7248p-1,0x1.6a09e667f3bcdp-1,0x1.4f2ec413cb52ap-1,0x1.3988e1409212ep-1 }
  ,{ 0x0p+0,0x1p+1,0x1.6a09e667f3bcdp+0,0x1.279a74590331cp+0,0x1p+0,0x1.c9f25c5bfedd9p-1,0x1.a20bd700c2c3ep-1,0x1.83091e6a7f7e6p-1,0x1.6a09e667f3bcdp-1 }
  ,{ 0x0p+0,0x1.1e3779b97f4a8p+1,0x1.94c583ada5b53p+0,0x1.4a7e9cb8a3491p+0,0x1.1e3779b97f4a8p+0,0x1p+0,0x1.d363d1848dcbfp-1,0x1.b0b80ef844ba1p-1,0x1.94c583ada5b53p-1 }
  ,{ 0x0p+0,0x1.3988e1409212ep+1,0x1.bb67ae8584caap+0,0x1.6a09e667f3bcdp+0,0x1.3988e1409212ep+0,0x1.186f174f88472p+0,0x1p+0,0x1.da05179501504p-1,0x1.bb67ae8584caap-1 }
  ,{ 0x0p+0,0x1.52a7fa9d2f8eap+1,0x1.deeea11683f49p+0,0x1.870be4c1c28b2p+0,0x1.52a7fa9d2f8eap+0,0x1.2ee73dadc9b57p+0,0x1.1482f86c40c43p+0,0x1p+0,0x1.deeea11683f49p-1 }
  ,{ 0x0p+0,0x1.6a09e667f3bcdp+1,0x1p+1,0x1.a20bd700c2c3ep+0,0x1.6a09e667f3bcdp+0,0x1.43d136248490fp+0,0x1.279a74590331cp+0,0x1.11acee560242ap+0,0x1p+0 }
  };
};

}

