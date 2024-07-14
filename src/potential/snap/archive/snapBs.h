#pragma once

#include <vector>
#include <map>
#include <mutex>
#include <fstream>

#include "snap_3Dtypes.h"
#include "snapCg.h"
#include "snapGsh.h"
#include "snap_vector.h"
#include "snap_ext.h"

/*
 * References:
 * - Bartok (PhD):
 *   Gaussian Approximation Potential: an interatomic potential derived from first principles Quantum Mechanics.
 * - Thompson (JCP):
 *   Spectral neighbor analysis method for automated generation of quantum-accurate interatomic potentials.
 */

//#include <iostream> // DEBUG
namespace SnapExt
{


struct snapBsIdxCache;

class snapBs
{
  using CMMData = SnapExt::CMMData;

  public:
    snapBs(double const jmax, const snapCg &cg, double const *coefs, double factor);

    int set_neighbours(double const *rx, double const *ry, double const *rz, const double * rd, double rcut, size_t N_atom);
  void compute_cmm(double rcut,double rfac0,double rmin0); 			// compute cmm and dcmm coefficients
    void compute_bs();		// compute bispectrum components and its derivatives
    void compute_bs0();		// compute bispectrum 0 components for an isolated atom with no neighbors  

    static int n_idx_bs(int m_two_jmax);// number of bispectrum components

    inline double en_val() const { return m_energy; }
    inline double en_zero_val() const { return m_energy_zero; }
  
    inline double3d force_val(int i) const { return m_force[i]; }
    inline int get_n_atom() const { return m_N_atom; }

    inline double get_jmax() const { return m_jmax; }
    inline double get_factor() const { return m_factor; }
    inline double get_coefs_0() const { return m_coefs[0]; }
    inline const snapBsIdxCache* get_idx_cache() const { return m_idx_cache; }

    void generate_constant_tables(std::ofstream& fout);
    void generate_lbs_compute_blocks(std::vector< SnapExt::BSFullBlockWorkItem >&  bs_fblock );
    
  private:
    int m_N_atom;			// Number of atoms in neighbourhood
//    int m_N_neigh;			// Number of atoms in neighbourhood closer than rcut
    int m_nidx;				// Number of bispectrum components
    int m_two_jmax;
    int m_gsh_size;			// number of gsh coefficients

    const double m_jmax;
    const double m_factor;		// input: weight factors

    const double* __restrict__ m_coefs = nullptr;		// input: machine learning coefficients (number: nb of bispectrum components * nb species)
    const double* __restrict__ m_rx = nullptr;
    const double* __restrict__ m_ry = nullptr;
    const double* __restrict__ m_rz = nullptr;
    const double* __restrict__ m_radius = nullptr;

    snapGsh m_gsh;		// general spherical harmonics

    //std::vector< Complexd > m_cmm;
    std::vector< CMMData > m_cmm;
    
    SnapVectorT<double3d> m_force;	// force for each neighbor
    double m_energy = 0.0;		// contribution to energy
    double m_energy_zero = 0.0;	// contribution to energy0 of isolated atom

    const snapCg& m_cg;

    // precomputed stuff
    snapBsIdxCache* m_idx_cache = nullptr;
};

}

