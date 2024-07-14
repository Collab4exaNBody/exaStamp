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

#include "pointers.h"

#define SNA_CU_ARRAY(a,i) a [ ( (i) * ONIKA_CU_BLOCK_SIZE ) + ONIKA_CU_THREAD_IDX ]

#define WJ(i)           SNA_CU_ARRAY(wj,i)
#define RCUTIJ(i)       SNA_CU_ARRAY(rcutij,i)
#define SINNERIJ(i)     SNA_CU_ARRAY(sinnerij,i)
#define DINNERIJ(i)     SNA_CU_ARRAY(dinnerij,i)
#define ELEMENT(i)      SNA_CU_ARRAY(element,i)
#define ULIST_IDIM      idxu_max
#define ULIST_R_IJ(j,i) SNA_CU_ARRAY(ulist_r_ij, (j) * ULIST_IDIM + (i) )
#define ULIST_I_IJ(j,i) SNA_CU_ARRAY(ulist_i_ij, (j) * ULIST_IDIM + (i) )
#define ULISTTOT_R(i)   SNA_CU_ARRAY(ulisttot_r,i)
#define ULISTTOT_I(i)   SNA_CU_ARRAY(ulisttot_i,i)
#define DULIST_IDIM     3
#define DULIST_R(j,i)   SNA_CU_ARRAY(dulist_r, (j) * DULIST_IDIM + (i) )
#define DULIST_I(j,i)   SNA_CU_ARRAY(dulist_i, (j) * DULIST_IDIM + (i) )
#define ZLIST_R(i)      SNA_CU_ARRAY(zlist_r,i)
#define ZLIST_I(i)      SNA_CU_ARRAY(zlist_i,i)
#define BLIST(i)        SNA_CU_ARRAY(blist,i)
#define YLIST_R(i)      SNA_CU_ARRAY(ylist_r,i)
#define YLIST_I(i)      SNA_CU_ARRAY(ylist_i,i)

namespace LAMMPS_NS
{

  struct SNA_ZINDICES {
    int j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, jju;
  };

  struct SNA_BINDICES {
    int j1, j2, j;
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
    SNA_BINDICES * __restrict__ idxb = nullptr;
    double * __restrict__ *  __restrict__ rootpqarray = nullptr;
    double *  __restrict__ cglist = nullptr;
    int *  __restrict__ *  __restrict__ *  __restrict__ idxcg_block = nullptr;
    int *  __restrict__ idxu_block = nullptr;
    int *  __restrict__ *  __restrict__ *  __restrict__ idxz_block = nullptr;
    int *  __restrict__ *  __restrict__ *  __restrict__ idxb_block = nullptr;
  };

  struct ReadOnlySnapParameters
  {
    inline ReadOnlySnapParameters( const SNA * sna )
      : bzero ( sna->bzero )
      , idxz ( sna->idxz )
      , idxb ( sna->idxb )
      , rootpqarray ( sna->rootpqarray )
      , cglist ( sna->cglist )
      , idxcg_block ( sna->idxcg_block )
      , idxu_block ( sna->idxu_block )
      , idxz_block ( sna->idxz_block )
      , idxb_block ( sna->idxb_block )
      , rmin0 ( sna->rmin0 )
      , rfac0 ( sna->rfac0 )
      , wself ( sna->wself )
      , nelements ( sna->nelements )
      , ndoubles ( sna->ndoubles )
      , ntriples ( sna->ntriples )
      , ncoeff ( sna->ncoeff )
      , twojmax ( sna->twojmax )
      , idxcg_max ( sna->idxcg_max )
      , idxu_max ( sna->idxu_max )
      , idxz_max ( sna->idxz_max )
      , idxb_max ( sna->idxb_max )
      , switch_flag ( sna->switch_flag )
      , switch_inner_flag ( sna->switch_inner_flag )
      , bzero_flag ( sna->bzero_flag )
      , bnorm_flag ( sna->bnorm_flag )
      , chem_flag ( sna->chem_flag )
      , wselfall_flag ( sna->wselfall_flag )
    {}
  
    double const * const __restrict__ bzero = nullptr;
    SNA_ZINDICES const * const __restrict__ idxz = nullptr;
    SNA_BINDICES const * const __restrict__ idxb = nullptr;
    double const * const __restrict__ * const __restrict__ rootpqarray = nullptr;
    double const * const __restrict__ cglist = nullptr;
    int const * const __restrict__ * const __restrict__ * const __restrict__ idxcg_block = nullptr;
    int const * const __restrict__ idxu_block = nullptr;
    int const * const __restrict__ * const __restrict__ * const __restrict__ idxz_block = nullptr;
    int const * const __restrict__ * const __restrict__ * const __restrict__ idxb_block = nullptr;

    double const rmin0 = 0.0;
    double const rfac0 = 0.0;
    double const wself = 0.0;

    int const nelements = 0;
    int const ndoubles = 0;
    int const ntriples = 0;
    int const ncoeff = 0;
    int const twojmax = 0;
    int const idxcg_max = 0;
    int const idxu_max = 0;
    int const idxz_max = 0;
    int const idxb_max = 0;

    bool const switch_flag =false;
    bool const switch_inner_flag = false;
    bool const bzero_flag = false; 
    bool const bnorm_flag = false; 
    bool const chem_flag = false; 
    bool const wselfall_flag = false; 
  };

  struct SnapScratchBuffers
  {    
    void initialize( const SNA * conf, size_t block_size=1 );
    void finalize();
    
    void grow_rij(int);
    void create_twojmax_arrays();
    void destroy_twojmax_arrays();

    const SNA* config = nullptr;
    Memory* memory = nullptr;
    
    int nmax = 0;    // allocated size of short lists

    int * __restrict__ element = nullptr;    // short element list [0,nmax)
    double * __restrict__ blist = nullptr;
    double * __restrict__  wj = nullptr;        // short weight list
    double * __restrict__ rcutij = nullptr;    // short cutoff list
    double * __restrict__ sinnerij = nullptr;    // short inner cutoff midpoint list
    double * __restrict__ dinnerij = nullptr;    // short inner half-width list

    double * __restrict__ ulisttot_r = nullptr;
    double * __restrict__ ulisttot_i = nullptr;
    double * __restrict__ ulist_r_ij = nullptr;
    double * __restrict__ ulist_i_ij = nullptr;
    double * __restrict__ zlist_r = nullptr;
    double * __restrict__ zlist_i = nullptr;

    double * __restrict__ dulist_r = nullptr;
    double * __restrict__ dulist_i = nullptr;
    double * __restrict__ ylist_r = nullptr;
    double * __restrict__ ylist_i = nullptr;
  };


}    // namespace LAMMPS_NS

