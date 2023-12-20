#pragma once

#include <vector>
#include <complex>
#include <exaStamp/potential/snaplegacy/SnapLegacyCG.h>
#include <exaStamp/potential/snaplegacy/SnapLegacyGSH.h>

/*
 * References:
 * - Bartok (PhD):
 *   Gaussian Approximation Potential: an interatomic potential derived from first principles Quantum Mechanics.
 * - Thompson (JCP):
 *   Spectral neighbor analysis method for automated generation of quantum-accurate interatomic potentials.
 */
//struct complex3d
//{
//	complex<double> x;
//	complex<double> y;
//	complex<double> z;
//};

class SnapLegacyBS
{
  template<typename T> using vector = std::vector<T>;
  template<typename T> using complex = std::complex<T>;

  public:
    SnapLegacyBS(double const jmax /*, double const rcut*/ /*, int myspecy*/, double const *coefs, double const *factor);

    int compute_bs(int specy, double rcut, const SnapLegacyCG &cg);		// compute bispectrum components and its derivatives

    int compute_cmm(double rcut); 			// compute cmm and dcmm coefficients

    int set_neighbours(double const *rx, double const *ry, double const *rz, int const * _species, double rcut, size_t N_atom); //list of neighbours positions, species, and number of neighbours

    static int n_idx_bs(int m_two_jmax);// number of bispectrum components

    complex<double> bs_val(int idx) {	// just in case and for tests
      return m_bs[idx];
    }

    complex3d dbs_val(int idx) {	// just in case and for tests
      return m_dbs[idx];
    }

    complex<double> en_val(int specy);		// SNAP energy of the atom
    complex3d force_val(int i);		// contribution to the force from ith atom of neighbourhood

//    inline double get_rcut() const { return m_rcut; }
//    inline int get_my_specy() const { return m_myspecy; }
    inline int get_n_atom() const { return m_N_atom; }

  private:
    int m_N_atom;			// Number of atoms in neighbourhood
    int m_N_neigh;			// Number of atoms in neighbourhood closer than rcut
    int m_nidx;				// Number of bispectrum components
    double m_jmax;
//    double m_rcut;
    //const int m_myspecy;			// input: atom specy
    const double* m_coefs;		// input: machine learning coefficients (number: nb of bispectrum components * nb species)
    const double* m_factor;		// input: weight factors
    int m_two_jmax;
    int m_gsh_size;			// number of gsh coefficients
    vector<complex<double> > m_bs;	// bispectrum components
    vector<complex3d> m_dbs;		// bispectrum derivatives
    vector<complex<double> > m_cmm;	// cmm coefficients (Thompson)
    vector<complex3d> m_dcmm;		// cmm coefficients derivatives (Thompson)
    vector<double> m_radius;		// relative distance of neighbours
    vector<double3d> m_r_vec;		// position of neighbours
    vector<SnapLegacyGSH> m_vec_gsh;		// general spherical harmonics objects vector
    complex<double> m_energy;		// contribution to energy
    vector<int> m_species;		// input: neighbours species

    inline int idx(int j, int m1, int m2) {
      return j*(j+1)*(2*j+1)/6+(m2+j)/2*(j+1)+(m1+j)/2+2; //sum of j first square integers+offset
    }

    inline int didx(int N, int j, int m1, int m2) {
      return N*m_gsh_size+j*(j+1)*(2*j+1)/6+(m2+j)/2*(j+1)+(m1+j)/2+2; //sum of j first square integers+offset
    }

    vector<double> m_default_coefs;
    vector<double> m_default_factors;  
};
