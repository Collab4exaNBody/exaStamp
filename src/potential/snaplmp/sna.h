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

#define WJ(i)           wj[i]
#define RCUTIJ(i)       rcutij[i]
#define SINNERIJ(i)     sinnerij[i]
#define DINNERIJ(i)     dinnerij[i]
#define ELEMENT(i)      element[i]
#define ULIST_IDIM      idxu_max
#define ULIST_R_IJ(j,i) ulist_r_ij[j][i] // ulist_r_ij[ (j)*ULIST_IDIM + (i) ]
#define ULIST_I_IJ(j,i) ulist_i_ij[j][i] // ulist_i_ij[ (j)*ULIST_IDIM + (i) ]
#define ULISTTOT_R(i)   ulisttot_r[i]
#define ULISTTOT_I(i)   ulisttot_i[i]
#define DULIST_IDIM     3
#define DULIST_R(j,i)   dulist_r[j][i] // dulist_r[ (j) * DULIST_IDIM + (i) ]
#define DULIST_I(j,i)   dulist_i[j][i] // dulist_i[ (j) * DULIST_IDIM + (i) ]
#define ZLIST_R(i)      zlist_r[i]
#define ZLIST_I(i)      zlist_i[i]
#define BLIST(i)        blist[i]
#define YLIST_R(i)      ylist_r[i]
#define YLIST_I(i)      ylist_i[i]

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

    double *bzero = nullptr;        // array of B values for isolated atoms
    SNA_ZINDICES *idxz = nullptr;
    SNA_BINDICES *idxb = nullptr;
    double **rootpqarray = nullptr;
    double *cglist = nullptr;
    int ***idxcg_block = nullptr;
    int *idxu_block = nullptr;
    int ***idxz_block = nullptr;
    int ***idxb_block = nullptr;
  };

  struct SnapScratchBuffers
  {
    SnapScratchBuffers( const SNA * conf, Memory* mem );
    ~SnapScratchBuffers();
    void grow_rij(int);
    void create_twojmax_arrays();
    void destroy_twojmax_arrays();

    const SNA* config = nullptr;
    Memory* memory = nullptr;
    
    int nmax = 0;    // allocated size of short lists
    int elem_duarray = 0;    // element of j in derivative

    int *element = nullptr;    // short element list [0,nmax)
    double *blist = nullptr;
//    double **rij = nullptr;      // short rij list
//    int *inside = nullptr;       // short neighbor list
    double *wj = nullptr;        // short weight list
    double *rcutij = nullptr;    // short cutoff list

    // only allocated for switch_inner_flag=1
    double *sinnerij = nullptr;    // short inner cutoff midpoint list
    double *dinnerij = nullptr;    // short inner half-width list

    double *ulisttot_r = nullptr;
    double *ulisttot_i = nullptr;
    double **ulist_r_ij = nullptr;
    double **ulist_i_ij = nullptr;
    double *zlist_r = nullptr;
    double *zlist_i = nullptr;

    double **dulist_r = nullptr;
    double **dulist_i = nullptr;
    double *ylist_r = nullptr;
    double *ylist_i = nullptr;
  };


}    // namespace LAMMPS_NS

