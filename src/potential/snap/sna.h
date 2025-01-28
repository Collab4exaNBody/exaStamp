/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Aidan Thompson, Christian Trott, SNL
------------------------------------------------------------------------- */

#pragma once

#include <onika/log.h>
#include <onika/cuda/cuda.h>
#include "pointers.h"
#include "sna_defs.h"
#include "snap_ymap.h"

#include <onika/integral_constant.h>

namespace SnapInternal
{

  struct SNA_ZINDICES {
    int16_t j1;
    int16_t j2;
    int16_t j;
    int16_t ma1min;
    int16_t ma2max;
    int16_t mb1min;
    int16_t mb2max;
    int16_t na;
    int16_t nb;
    int16_t jju;
  };

  struct alignas(8) SNA_ZINDICES_ALT
  {
    int16_t j1;
    int16_t j2;
    int16_t j;
    int16_t ma1min;
    int16_t ma2max;
    int16_t mb1min;
    int16_t mb2max;
    int16_t na;
    int16_t nb;
    int16_t jju;
    int16_t jjb;
    int16_t idx_cgblock;
  };

  static_assert( sizeof(SNA_ZINDICES_ALT)==24 && alignof(SNA_ZINDICES_ALT)==8 );

  struct JJUPairBlend
  {
    int32_t jju1;
    int32_t jju2;
    double coeff;
  };

  struct SNA_BINDICES {
    int16_t j1;
    int16_t j2;
    int16_t j;
  };

  struct SNA
  {
    SNA(Memory *, double, int, double, int, int, int, int, int, int, int);
    
    Memory* memory = nullptr;

    //SNA(LAMMPS *lmp) : Pointers(lmp){};
    ~SNA();
    void build_indexlist();
    void init();

    void create_twojmax_arrays();
    void destroy_twojmax_arrays();
    void init_clebsch_gordan();
    void print_clebsch_gordan();
    void init_rootpqarray();
    double deltacg(int, int, int);
    void compute_ncoeff();

    double rmin0 = 0.0;
    double rfac0 = 0.0;
    double wself = 0.0;

    int nelements = 0;        // number of elements
    int ndoubles = 0;         // number of multi-element pairs
    int ntriples = 0;         // number of multi-element triplets
    int ncoeff = 0;
    int twojmax = 0;
    int idxcg_max = 0;
    int idxu_max = 0;
    int idxz_max = 0;
    int idxb_max = 0;
    int idxz_max_alt = 0;
    int idxu_max_alt = 0;

    // Sets the style for the switching function
    // 0 = none
    // 1 = cosine
    int switch_flag = 0;

    // Sets the style for the inner switching function
    // 0 = none
    // 1 = cosine
    int switch_inner_flag = 0;

    int bzero_flag = 0;       // 1 if bzero subtracted from barray
    int bnorm_flag = 0;       // 1 if barray divided by j+1
    int chem_flag = 0;        // 1 for multi-element bispectrum components
    int wselfall_flag = 0;    // 1 for adding wself to all element labelings

    double * __restrict__ bzero = nullptr;        // array of B values for isolated atoms
    SNA_ZINDICES * __restrict__ idxz = nullptr;
    SNA_ZINDICES_ALT * __restrict__ idxz_alt = nullptr;
    SNA_BINDICES * __restrict__ idxb = nullptr;
    double * __restrict__ rootpqarray = nullptr;
    double *  __restrict__ cglist = nullptr;
    int *  __restrict__ idxcg_block = nullptr;
    int *  __restrict__ idxz_block = nullptr;
    int *  __restrict__ idxb_block = nullptr;
    int *  __restrict__ y_jju_map = nullptr; // size is idxu_max, stored indices are in [0;idxu_max_alt[
  };

  static inline constexpr int compute_idxcg_max(int twojmax)
  {
    int idxcg_count = 0;
    for (int j1 = 0; j1 <= twojmax; j1++)
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
          for (int m1 = 0; m1 <= j1; m1++)
            for (int m2 = 0; m2 <= j2; m2++)
              idxcg_count++;
        }
    return idxcg_count;
  }
  template<int twojmax> static inline constexpr auto compute_idxcg_max( onika::IntConst<twojmax> ) { return onika::IntConst<compute_idxcg_max(twojmax)>{}; }

  static inline constexpr int compute_idxb_max(int twojmax)
  {
    int idxb_count = 0;
    for (int j1 = 0; j1 <= twojmax; j1++)
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
          if (j >= j1) idxb_count++;
    return idxb_count;
  }
  template<int twojmax> static inline constexpr auto compute_idxb_max( onika::IntConst<twojmax> ) { return onika::IntConst<compute_idxb_max(twojmax)>{}; }

  static inline constexpr int compute_idxz_max(int twojmax)
  {
    int idxz_count = 0;
    for (int j1 = 0; j1 <= twojmax; j1++)
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2)
          for (int mb = 0; 2*mb <= j; mb++)
            for (int ma = 0; ma <= j; ma++)
              idxz_count++;
    return idxz_count;
  }
  template<int twojmax> static inline constexpr auto compute_idxz_max( onika::IntConst<twojmax> ) { return onika::IntConst<compute_idxz_max(twojmax)>{}; }

  static inline constexpr int compute_idxz_max_alt(int twojmax)
  {
    int idxz_count = 0;
    for (int j1 = 0; j1 <= twojmax; j1++)
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2) {
          for (int mb = 0; 2*mb <= j; mb++)
            for (int ma = 0; ma <= j; ma++) {
              const int jju = IDXU_BLOCK(j) + (j+1)*mb + ma;
              if( exaStamp::snap_force_use_Y(twojmax,jju) ) idxz_count++;
            }
        }
    return idxz_count;
  }  
  template<int twojmax>
  static inline constexpr auto compute_idxz_max_alt( onika::IntConst<twojmax> )
  {
    return onika::IntConst<compute_idxz_max_alt(twojmax)>{};
  }

  static inline constexpr int compute_idxu_max_alt(int twojmax, int idxu_max)
  {
    return exaStamp::snap_force_Y_count(twojmax,idxu_max);
  }
  template<int twojmax, int idxu_max>
  static inline constexpr auto compute_idxu_max_alt( onika::IntConst<twojmax> , onika::IntConst<idxu_max> )
  {
    return onika::IntConst<compute_idxu_max_alt(twojmax,idxu_max)>{};
  }


  static inline constexpr int compute_coeff_count(int twojmax, int ntriples)
  {
    int ncount = 0;
    for (int j1 = 0; j1 <= twojmax; j1++)
      for (int j2 = 0; j2 <= j1; j2++)
        for (int j = j1 - j2;
             j <= std::min(twojmax, j1 + j2); j += 2)
          if (j >= j1) ncount++;
    return ncount * ntriples;
  }
  template<int twojmax,int ntriples> static inline constexpr auto compute_coeff_count( onika::IntConst<twojmax> , onika::IntConst<ntriples> ) { return onika::IntConst<compute_coeff_count(twojmax,ntriples)>{}; }

  static inline constexpr int sum_int_sqr(int X) { return SUM_INT_SQR(X); }
  template<int X> static inline constexpr auto sum_int_sqr( onika::IntConst<X> ) { return onika::IntConst< SUM_INT_SQR(X) >{}; }


  template<class T> struct InitialValue
  {
    static inline T check( const T& x ) { return x; }
  };
  template<class T, T value> struct InitialValue< onika::IntegralConst<T,value> >
  {
    static inline onika::IntegralConst<T,value> check( const T& chkval ) { if( value != chkval ) exanb::fatal_error() << "bad initialization value" << std::endl; return {}; }
  };

  template< class JMaxT , class NElementsT >
  struct ReadOnlySnapParameters
  {
#   define DEDUCED_AUTO_CONST_MEMBER(name,expr) const decltype(expr) name = expr

    const JMaxT jmax = {};
    DEDUCED_AUTO_CONST_MEMBER( twojmax , jmax * onika::IntConst<2>{} );
    DEDUCED_AUTO_CONST_MEMBER( idxcg_max , compute_idxcg_max(twojmax) );
    DEDUCED_AUTO_CONST_MEMBER( idxu_max , sum_int_sqr( twojmax + onika::IntConst<1>{} ) );
    DEDUCED_AUTO_CONST_MEMBER( idxz_max , compute_idxz_max(twojmax) );
    DEDUCED_AUTO_CONST_MEMBER( idxb_max , compute_idxb_max(twojmax) );
    DEDUCED_AUTO_CONST_MEMBER( idxz_max_alt , compute_idxz_max_alt(twojmax) );
    DEDUCED_AUTO_CONST_MEMBER( idxu_max_alt , compute_idxu_max_alt(twojmax,idxu_max) );
    
    const NElementsT nelements = {};
    DEDUCED_AUTO_CONST_MEMBER( ndoubles , nelements * nelements );
    DEDUCED_AUTO_CONST_MEMBER( ntriples , nelements * nelements * nelements );
    DEDUCED_AUTO_CONST_MEMBER( ncoeff , compute_coeff_count( twojmax , ntriples ) );

#   undef DEDUCED_AUTO_CONST_MEMBER
  
    ReadOnlySnapParameters() = default;
  
    inline ReadOnlySnapParameters( const SNA * sna )
      : jmax( InitialValue<JMaxT>::check( sna->twojmax / 2 ) )
      , nelements( InitialValue<NElementsT>::check( sna->nelements ) )
      , bzero ( sna->bzero )
      , idxz ( sna->idxz )
      , idxz_alt ( sna->idxz_alt )
      , idxb ( sna->idxb )
      , rootpqarray ( sna->rootpqarray )
      , cglist ( sna->cglist )
      , idxcg_block ( sna->idxcg_block )
      , idxz_block ( sna->idxz_block )
      , idxb_block ( sna->idxb_block )
      , y_jju_map( sna->y_jju_map )
      , rmin0 ( sna->rmin0 )
      , rfac0 ( sna->rfac0 )
      , wself ( sna->wself )
      , switch_flag ( sna->switch_flag )
      , switch_inner_flag ( sna->switch_inner_flag )
      , bzero_flag ( sna->bzero_flag )
      , bnorm_flag ( sna->bnorm_flag )
      , chem_flag ( sna->chem_flag )
      , wselfall_flag ( sna->wselfall_flag )
    {
      if( twojmax   != sna->twojmax   ) exanb::fatal_error() << "twojmax : sna value ("<<sna->twojmax<<") != "<<twojmax<<std::endl;
      if( idxcg_max != sna->idxcg_max ) exanb::fatal_error() << "idxcg_max : sna value ("<<sna->idxcg_max<<") != "<<idxcg_max<<std::endl;
      if( idxu_max  != sna->idxu_max  ) exanb::fatal_error() << "idxu_max : sna value ("<<sna->idxu_max<<") != "<<idxu_max<<std::endl;
      if( idxz_max  != sna->idxz_max  ) exanb::fatal_error() << "idxz_max : sna value ("<<sna->idxz_max<<") != "<<idxz_max<<std::endl;
      if( idxb_max  != sna->idxb_max  ) exanb::fatal_error() << "idxb_max : sna value ("<<sna->idxb_max<<") != "<<idxb_max<<std::endl;
      if( idxz_max_alt != sna->idxz_max_alt ) exanb::fatal_error() << "idxz_max_alt : sna value ("<<sna->idxz_max_alt<<") != "<<idxz_max_alt<<std::endl;
      if( idxu_max_alt != sna->idxu_max_alt ) exanb::fatal_error() << "idxu_max_alt : sna value ("<<sna->idxu_max_alt<<") != "<<idxu_max_alt<<std::endl;
      if( nelements != sna->nelements ) exanb::fatal_error() << "nelements : sna value ("<<sna->nelements<<") != "<<nelements<<std::endl;
      if( ndoubles  != sna->ndoubles  ) exanb::fatal_error() << "ndoubles : sna value ("<<sna->ndoubles<<") != "<<ndoubles<<std::endl;
      if( ntriples  != sna->ntriples  ) exanb::fatal_error() << "nelements : sna value ("<<sna->ntriples<<") != "<<ntriples<<std::endl;
      if( ncoeff    != sna->ncoeff    ) exanb::fatal_error() << "ncoeff : sna value ("<<sna->ncoeff<<") != "<<ncoeff<<std::endl;
    }
  
    double const * const __restrict__ bzero = nullptr;
    SNA_ZINDICES const * const __restrict__ idxz = nullptr;
    SNA_ZINDICES_ALT const * const __restrict__ idxz_alt = nullptr;
    SNA_BINDICES const * const __restrict__ idxb = nullptr;
    double const * const __restrict__ rootpqarray = nullptr;
    double const * const __restrict__ cglist = nullptr;
    int const * const __restrict__ idxcg_block = nullptr;
    int const * const __restrict__ idxz_block = nullptr;
    int const * const __restrict__ idxb_block = nullptr;
    int const * const __restrict__ y_jju_map = nullptr;

    double const rmin0 = 0.0;
    double const rfac0 = 0.0;
    double const wself = 0.0;

    bool const switch_flag =false;
    bool const switch_inner_flag = false;
    bool const bzero_flag = false; 
    bool const bnorm_flag = false; 
    bool const chem_flag = false; 
    bool const wselfall_flag = false;
    
    template<class StreamT>
    inline StreamT& to_stream(StreamT& out) const
    {
      out << "flags: switch="<<std::boolalpha<<switch_flag<<" switch_inner="<<switch_inner_flag<<" bzero="<<bzero_flag<<" bnorm="<<bnorm_flag<<" chem="<<chem_flag<<" wselfall="<<wselfall_flag<<std::endl;
      out << "const: rmin0="<<rmin0<<" rfac0="<<rfac0<<" wself="<<wself<<std::endl;
      out << "sizes: jmax="<<jmax<<" idxcg_max="<<idxcg_max<<" idxu_max="<<idxu_max<<" idxz_max="<<idxz_max<<" idxb_max="
          << idxb_max<<" nelements="<<nelements<<" ndoubles="<<ndoubles<<" ntriples="<<ntriples<<" ncoeff="<<ncoeff<<std::endl;

      std::unordered_set<int> unique_idxcg;
      std::unordered_set<int> unique_jju;
      out << "IDXZ:" << std::endl;
      for(int jjz=0;jjz<idxz_max;jjz++)
      {
        const int j1 = IDXZ_ALT(jjz).j1;
        const int j2 = IDXZ_ALT(jjz).j2;
        const int j = IDXZ_ALT(jjz).j;
        //const int ma1min = IDXZ_ALT(jjz).ma1min;
        //const int ma2max = IDXZ_ALT(jjz).ma2max;
        //const int na = IDXZ_ALT(jjz).na;
        //const int mb1min = IDXZ_ALT(jjz).mb1min;
        //const int mb2max = IDXZ_ALT(jjz).mb2max;
        //const int nb = IDXZ_ALT(jjz).nb;
        const int jju = IDXZ_ALT(jjz).jju;

        const int idxcg_block_j1_j2_j = IDXCG_BLOCK(j1,j2,j);
        int jjb = -1;
        int betaj_scale = 0;
        if (j >= j1) {
          jjb = IDXB_BLOCK(j1,j2,j);
          if (j1 == j) {
            if (j2 == j) betaj_scale = 3;
            else betaj_scale = 2;
          } else betaj_scale = 1;
        } else if (j >= j2) {
          jjb = IDXB_BLOCK(j,j2,j1);
          if (j2 == j) betaj_scale = 2;
          else betaj_scale = 1;
        } else {
          jjb = IDXB_BLOCK(j2,j,j1);
          betaj_scale = 1;
        }
        unique_idxcg.insert( idxcg_block_j1_j2_j );
        unique_jju.insert( jju );
        out << "  jju=" <<jju<<" jjb="<<jjb<<" betascale="<<betaj_scale<<std::endl;
      }
      out << "unique idxcg values = "<<unique_idxcg.size()<<std::endl;
      out << "unique jju values = "<<unique_jju.size()<<std::endl;
      return out;
    }
  };

}    // namespace SnapInternal

