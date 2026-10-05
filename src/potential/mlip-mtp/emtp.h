/* ----------------------------------------------------------------------
   Standalone MTP (Moment Tensor Potential) engine — decoupled from LAMMPS.
   Ported from LAMMPS ML-MTP (pair_mtp.cpp, mtp_rb_chevbyshev_basis.cpp) by
   Richard Meng et al. — see /home/lafourcadep/CODES/ATOMISTIC/FFs/lammps-mtp-kokkos.
------------------------------------------------------------------------- */

#pragma once

#include <string>
#include <vector>

// Reads a real MLIP-3-format .mtp/.almtp file (topology + radial-basis coefficients +
// linear model coefficients) and evaluates, per central atom, either the trained
// energy/force (peratom_energyforce_soa) or the coefficient-free descriptor vector and its
// per-neighbor position-derivative (peratom_descriptors_soa). Only RBChebyshev radial basis
// is supported (the only type LAMMPS's own reference implementation supports). The MaxVol
// active-learning ("#MVS_v1.1") block is intentionally not parsed — out of scope for v1, and
// LAMMPS's own reader never looks for it either.
class EMTP
{
public:
  // ---- parsed potential data ----
  double scaling = 1.0;
  int species_count = 1;
  double min_dist = 0.0, max_dist = 5.0;
  int radial_basis_size = 0;
  int radial_funcs_count = 0;
  // [ (type1*species_count+type2) * radial_funcs_count*radial_basis_size + mu*radial_basis_size + ri ]
  std::vector<double> radial_basis_coeffs;

  int alpha_moments_count = 0;
  int alpha_index_basic_count = 0;
  std::vector<int> alpha_index_basic;   // flat [alpha_index_basic_count][4]: mu,p,q,s
  int alpha_index_times_count = 0;
  std::vector<int> alpha_index_times;   // flat [alpha_index_times_count][4]: a0,a1,multiplier,a3
  int alpha_scalar_moments = 0;
  std::vector<int> alpha_moment_mapping;  // [alpha_scalar_moments]
  std::vector<double> species_coeffs;     // [species_count]
  std::vector<double> linear_coeffs;      // [alpha_scalar_moments] == file's "moment_coeffs", shared across species

  double rcut = 0.0;   // == max_dist, exposed for mtp_init's rcut_max wiring (mirrors EAPOD::rcut)

  int Njmax = 0;

  // Per-atom scratch/output, sized at construction to Njmax. Public so the compute-op
  // functors can read them straight after a call, same convention as EAPOD's soa_fij/bd/bdd.
  std::vector<double> soa_fij;  // [3*Njmax] -- energy/force path output (central-atom-owns/
                                 // neighbor-negated convention: caller does f_central += soa_fij[jj],
                                 // f_neighbor -= soa_fij[jj], matching pair_mtp.cpp's f[i]+=/f[j]-=)
  std::vector<double> bd;       // [alpha_scalar_moments] -- descriptor path: basis-function values B_k
  std::vector<double> bdd;      // descriptor path: per-neighbor-pair Jacobian of B_k, addressed
                                 // [3*jj + 3*Nj*k] using the CALL's actual Nj as stride (not Njmax) --
                                 // matches EAPOD's own bdd addressing convention exactly.

  EMTP(const std::string& mtp_file, int njmax);

  void read_mtp_file(const std::string& mtp_file);

  // Energy+force path. ti_0indexed/tj_0indexed are exaStamp 0-indexed particle types;
  // type_map converts them to MTP's own 0-indexed species index (purely positional --
  // MTP files carry no species names, see mtp_config.h). Returns per-atom energy.
  double peratom_energyforce_soa(const double* drx, const double* dry, const double* drz,
                                  int ti_0indexed, const int* tj_0indexed, int Nj,
                                  const int* type_map);

  // Coefficient-free descriptor path. Unlike POD's peratombase_descriptors_soa, MTP's
  // descriptor genuinely depends on the central atom's type too (the radial basis is
  // per-(central,neighbor)-species-pair) -- so ti_0indexed is required here.
  void peratom_descriptors_soa(const double* drx, const double* dry, const double* drz,
                                int ti_0indexed, const int* tj_0indexed, int Nj,
                                const int* type_map);

private:
  // Per-thread reusable scratch (mirrors EAPOD's own private-tmpmem strategy).
  std::vector<int> soa_tj_;                   // [Njmax] -- neighbor types mapped to MTP species
  std::vector<double> Q_, Qd_;                // [radial_basis_size] -- Chebyshev basis value/derivative
  std::vector<double> radial_vals_, radial_ders_;  // [radial_funcs_count]
  std::vector<double> dist_powers_;           // [max_rank_+1]
  std::vector<double> coord_powers_;          // [(max_rank_+1)*3]
  std::vector<double> moment_tensor_vals_;    // [alpha_moments_count]
  std::vector<double> moment_jacobian_;       // [Njmax*alpha_index_basic_count*3], stride = actual Nj per call
  std::vector<double> nbh_energy_ders_wrt_moments_;  // [alpha_moments_count]
  std::vector<double> moment_ders_wrt_basis_; // [alpha_moments_count*alpha_scalar_moments] -- descriptor-mode backprop matrix
  int max_rank_ = 0;   // max(p+q+s) over alpha_index_basic (dist_powers_/coord_powers_ need indices 0..max_rank_)

  // Direct, unmodified port of pair_mtp.cpp:110-201 (radial basis -> basic moments -> contraction
  // DAG forward pass), shared by both energy/force and descriptor modes. Fills moment_tensor_vals_
  // and moment_jacobian_ (basic-moment Jacobian only) for one atom's full neighbor slice. ti0/tj0
  // are already-mapped MTP 0-indexed species.
  void build_moments_and_jacobian(const double* drx, const double* dry, const double* drz,
                                   int ti0, const int* tj0, int Nj);
};
