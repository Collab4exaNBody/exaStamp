/* ----------------------------------------------------------------------
   Standalone MTP engine — decoupled from LAMMPS.
   Ported from LAMMPS ML-MTP (pair_mtp.cpp, mtp_rb_chevbyshev_basis.cpp).
------------------------------------------------------------------------- */

#include "emtp.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

// ── File-parsing helpers ──────────────────────────────────────────────────
// MLIP-3 .mtp files are mostly "key = value" lines (tab-indented nested keys are just
// whitespace to a tokenizer), except a handful of fields packed as one long brace-delimited
// line ("key = {{a,b,c,d}, {a,b,c,d}, ...}"). Replacing '=', '{', '}', ',', '-' (only used for
// the "type1-type2" radial_coeffs pair header) with spaces before whitespace-tokenizing handles
// every line uniformly -- no need for a stateful tokenizer class.
static std::vector<std::string> mtp_words(std::string line, const std::string& extra_seps)
{
  for (char& c : line) {
    if (c == '=' || extra_seps.find(c) != std::string::npos) c = ' ';
  }
  std::istringstream iss(line);
  std::vector<std::string> words;
  std::string w;
  while (iss >> w) words.push_back(w);
  return words;
}

EMTP::EMTP(const std::string& mtp_file, int njmax) : Njmax(njmax)
{
  read_mtp_file(mtp_file);

  soa_fij.assign(static_cast<size_t>(Njmax) * 3, 0.0);
  bd.assign(alpha_scalar_moments, 0.0);
  bdd.assign(static_cast<size_t>(Njmax) * alpha_scalar_moments * 3, 0.0);

  soa_tj_.assign(Njmax, 0);
  Q_.assign(radial_basis_size, 0.0);
  Qd_.assign(radial_basis_size, 0.0);
  radial_vals_.assign(radial_funcs_count, 0.0);
  radial_ders_.assign(radial_funcs_count, 0.0);
  dist_powers_.assign(max_rank_ + 1, 0.0);
  coord_powers_.assign(static_cast<size_t>(max_rank_ + 1) * 3, 0.0);
  // index 0 is always 1 regardless of dist/direction -- set once, matches LAMMPS's own
  // one-time init (pair_mtp.cpp:647), not reset per neighbor.
  dist_powers_[0] = 1.0;
  coord_powers_[0] = coord_powers_[1] = coord_powers_[2] = 1.0;

  moment_tensor_vals_.assign(alpha_moments_count, 0.0);
  moment_jacobian_.assign(static_cast<size_t>(Njmax) * alpha_index_basic_count * 3, 0.0);
  nbh_energy_ders_wrt_moments_.assign(alpha_moments_count, 0.0);
  moment_ders_wrt_basis_.assign(static_cast<size_t>(alpha_moments_count) * alpha_scalar_moments, 0.0);
}

void EMTP::read_mtp_file(const std::string& mtp_file)
{
  std::ifstream fp(mtp_file);
  if (!fp.is_open())
    throw std::runtime_error("Cannot open MTP file: " + mtp_file);

  auto next_line = [&]() -> std::string {
    std::string line;
    if (!std::getline(fp, line))
      throw std::runtime_error("MTP file ended unexpectedly while parsing: " + mtp_file);
    return line;
  };
  auto words_of = [&](const std::string& extra_seps) {
    return mtp_words(next_line(), extra_seps);
  };

  std::vector<std::string> words = words_of("");
  if (words.empty() || words[0] != "MTP")
    throw std::runtime_error("Not an MTP potential file (missing 'MTP' header): " + mtp_file);

  words = words_of("");
  if (words.empty() || words[0] != "version")
    throw std::runtime_error("MTP file missing 'version': " + mtp_file);

  words = words_of("");
  std::string keyword = words.empty() ? "" : words[0];
  if (keyword == "potential_name") { words = words_of(""); keyword = words.empty() ? "" : words[0]; }

  if (keyword == "scaling") {
    scaling = std::stod(words.at(1));
    words = words_of(""); keyword = words.empty() ? "" : words[0];
  } else {
    scaling = 1.0;
  }

  if (keyword != "species_count")
    throw std::runtime_error("MTP file: 'species_count' not found: " + mtp_file);
  species_count = std::stoi(words.at(1));

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword == "potential_tag") { words = words_of(""); keyword = words.empty() ? "" : words[0]; }

  if (keyword != "radial_basis_type")
    throw std::runtime_error("MTP file: 'radial_basis_type' not found: " + mtp_file);
  if (words.size() < 2 || words[1] != "RBChebyshev")
    throw std::runtime_error("MTP file: unsupported radial_basis_type (only RBChebyshev is implemented): " + mtp_file);

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  // Nested radial-basis "scaling" (if present) is superseded by the outer one -- LAMMPS itself
  // unconditionally overwrites radial_basis->scaling with the outer scaling right after
  // construction (pair_mtp.cpp:416), so its value is discarded here too, only its line consumed.
  if (keyword == "scaling") { words = words_of(""); keyword = words.empty() ? "" : words[0]; }

  if (keyword != "min_val" && keyword != "min_dist")
    throw std::runtime_error("MTP file: 'min_dist' not found: " + mtp_file);
  min_dist = std::stod(words.at(1));

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword != "max_val" && keyword != "max_dist")
    throw std::runtime_error("MTP file: 'max_dist' not found: " + mtp_file);
  max_dist = std::stod(words.at(1));
  rcut = max_dist;

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword != "radial_basis_size")
    throw std::runtime_error("MTP file: 'radial_basis_size' not found: " + mtp_file);
  radial_basis_size = std::stoi(words.at(1));

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword != "radial_funcs_count")
    throw std::runtime_error("MTP file: 'radial_funcs_count' not found: " + mtp_file);
  radial_funcs_count = std::stoi(words.at(1));

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword != "radial_coeffs")
    throw std::runtime_error("MTP file: 'radial_coeffs' not found (magnetic_basis_type is not supported): " + mtp_file);

  const int pairs_count = species_count * species_count;
  const int coeffs_per_pair = radial_basis_size * radial_funcs_count;
  radial_basis_coeffs.assign(static_cast<size_t>(pairs_count) * coeffs_per_pair, 0.0);

  for (int i = 0; i < pairs_count; i++) {
    std::vector<std::string> hdr = words_of("-");
    if (hdr.size() < 2)
      throw std::runtime_error("MTP file: malformed radial_coeffs pair header: " + mtp_file);
    const int type1 = std::stoi(hdr[0]);
    const int type2 = std::stoi(hdr[1]);
    const int pair_offset = (type1 * species_count + type2) * coeffs_per_pair;
    for (int j = 0; j < radial_funcs_count; j++) {
      std::vector<std::string> row = words_of("{,}");
      for (int k = 0; k < radial_basis_size; k++)
        radial_basis_coeffs[pair_offset + j * radial_basis_size + k] = std::stod(row.at(k));
    }
  }

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword != "alpha_moments_count")
    throw std::runtime_error("MTP file: 'alpha_moments_count' not found: " + mtp_file);
  alpha_moments_count = std::stoi(words.at(1));

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword != "alpha_index_basic_count")
    throw std::runtime_error("MTP file: 'alpha_index_basic_count' not found: " + mtp_file);
  alpha_index_basic_count = std::stoi(words.at(1));

  words = words_of("{},");
  if (words.empty() || words[0] != "alpha_index_basic")
    throw std::runtime_error("MTP file: 'alpha_index_basic' not found: " + mtp_file);
  alpha_index_basic.assign(static_cast<size_t>(alpha_index_basic_count) * 4, 0);
  for (int i = 0; i < alpha_index_basic_count; i++)
    for (int j = 0; j < 4; j++)
      alpha_index_basic[i * 4 + j] = std::stoi(words.at(1 + i * 4 + j));

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword != "alpha_index_times_count")
    throw std::runtime_error("MTP file: 'alpha_index_times_count' not found: " + mtp_file);
  alpha_index_times_count = std::stoi(words.at(1));

  words = words_of("{},");
  if (words.empty() || words[0] != "alpha_index_times")
    throw std::runtime_error("MTP file: 'alpha_index_times' not found: " + mtp_file);
  alpha_index_times.assign(static_cast<size_t>(alpha_index_times_count) * 4, 0);
  for (int i = 0; i < alpha_index_times_count; i++)
    for (int j = 0; j < 4; j++)
      alpha_index_times[i * 4 + j] = std::stoi(words.at(1 + i * 4 + j));

  words = words_of(""); keyword = words.empty() ? "" : words[0];
  if (keyword != "alpha_scalar_moments")
    throw std::runtime_error("MTP file: 'alpha_scalar_moments' not found: " + mtp_file);
  alpha_scalar_moments = std::stoi(words.at(1));

  words = words_of("{},");
  if (words.empty() || words[0] != "alpha_moment_mapping")
    throw std::runtime_error("MTP file: 'alpha_moment_mapping' not found: " + mtp_file);
  alpha_moment_mapping.assign(alpha_scalar_moments, 0);
  for (int i = 0; i < alpha_scalar_moments; i++) alpha_moment_mapping[i] = std::stoi(words.at(1 + i));

  words = words_of("{},");
  if (words.empty() || words[0] != "species_coeffs")
    throw std::runtime_error("MTP file: 'species_coeffs' not found: " + mtp_file);
  species_coeffs.assign(species_count, 0.0);
  for (int i = 0; i < species_count; i++) species_coeffs[i] = std::stod(words.at(1 + i));

  words = words_of("{},");
  if (words.empty() || words[0] != "moment_coeffs")
    throw std::runtime_error("MTP file: 'moment_coeffs' not found: " + mtp_file);
  linear_coeffs.assign(alpha_scalar_moments, 0.0);
  for (int i = 0; i < alpha_scalar_moments; i++) linear_coeffs[i] = std::stod(words.at(1 + i));

  // Stop here -- matches LAMMPS's own reader (pair_mtp.cpp:556-569), which never looks for the
  // trailing "#MVS_v1.1" MaxVol active-learning block either. Out of scope for v1.

  max_rank_ = 0;
  for (int i = 0; i < alpha_index_basic_count; i++) {
    const int rank = alpha_index_basic[i * 4 + 1] + alpha_index_basic[i * 4 + 2] + alpha_index_basic[i * 4 + 3];
    if (rank > max_rank_) max_rank_ = rank;
  }
}

void EMTP::build_moments_and_jacobian(const double* drx, const double* dry, const double* drz,
                                       int ti0, const int* tj0, int Nj)
{
  std::fill(moment_tensor_vals_.begin(), moment_tensor_vals_.end(), 0.0);

  const double inv_span = 2.0 / (max_dist - min_dist);
  const double mid = min_dist + max_dist;

  for (int jj = 0; jj < Nj; jj++) {
    const double r[3] = { drx[jj], dry[jj], drz[jj] };
    const double dist = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

    // Chebyshev radial basis + derivative (mtp_rb_chevbyshev_basis.cpp:29-53).
    const double ksi = (2.0 * dist - mid) * (1.0 / (max_dist - min_dist));
    const double dmc = dist - max_dist;
    Q_[0] = scaling * dmc * dmc;
    Qd_[0] = scaling * 2.0 * dmc;
    if (radial_basis_size > 1) {
      Q_[1] = scaling * ksi * dmc * dmc;
      Qd_[1] = scaling * (inv_span * dmc * dmc + 2.0 * ksi * dmc);
    }
    for (int i = 2; i < radial_basis_size; i++) {
      Q_[i] = 2.0 * ksi * Q_[i - 1] - Q_[i - 2];
      Qd_[i] = 2.0 * (inv_span * Q_[i - 1] + ksi * Qd_[i - 1]) - Qd_[i - 2];
    }

    const int jtype = tj0[jj];
    const int pair_offset = (ti0 * species_count + jtype) * radial_funcs_count * radial_basis_size;
    for (int mu = 0; mu < radial_funcs_count; mu++) {
      double val = 0.0, der = 0.0;
      const int off = pair_offset + mu * radial_basis_size;
      for (int ri = 0; ri < radial_basis_size; ri++) {
        val += radial_basis_coeffs[off + ri] * Q_[ri];
        der += radial_basis_coeffs[off + ri] * Qd_[ri];
      }
      radial_vals_[mu] = val;
      radial_ders_[mu] = der;
    }

    for (int k = 1; k <= max_rank_; k++) {
      dist_powers_[k] = dist_powers_[k - 1] * dist;
      coord_powers_[k * 3 + 0] = coord_powers_[(k - 1) * 3 + 0] * r[0];
      coord_powers_[k * 3 + 1] = coord_powers_[(k - 1) * 3 + 1] * r[1];
      coord_powers_[k * 3 + 2] = coord_powers_[(k - 1) * 3 + 2] * r[2];
    }

    for (int k = 0; k < alpha_index_basic_count; k++) {
      const int mu = alpha_index_basic[k * 4 + 0];
      const int p = alpha_index_basic[k * 4 + 1];
      const int q = alpha_index_basic[k * 4 + 2];
      const int s = alpha_index_basic[k * 4 + 3];
      const int rank = p + q + s;

      double val = radial_vals_[mu];
      double der = radial_ders_[mu];
      const double norm_fac = 1.0 / dist_powers_[rank];
      val *= norm_fac;
      der = der * norm_fac - rank * val / dist;

      const double pow0 = coord_powers_[p * 3 + 0];
      const double pow1 = coord_powers_[q * 3 + 1];
      const double pow2 = coord_powers_[s * 3 + 2];
      double powv = pow0 * pow1 * pow2;
      moment_tensor_vals_[k] += val * powv;

      powv *= der / dist;
      double jx = powv * r[0];
      double jy = powv * r[1];
      double jz = powv * r[2];
      if (p != 0) jx += val * p * coord_powers_[(p - 1) * 3 + 0] * pow1 * pow2;
      if (q != 0) jy += val * q * pow0 * coord_powers_[(q - 1) * 3 + 1] * pow2;
      if (s != 0) jz += val * s * pow0 * pow1 * coord_powers_[(s - 1) * 3 + 2];

      const size_t jbase = (static_cast<size_t>(jj) * alpha_index_basic_count + k) * 3;
      moment_jacobian_[jbase + 0] = jx;
      moment_jacobian_[jbase + 1] = jy;
      moment_jacobian_[jbase + 2] = jz;
    }
  }

  // Contraction DAG forward pass (pair_mtp.cpp:196-201) -- a flat instruction tape read entirely
  // from the file, no hardcoded contraction scheme.
  for (int k = 0; k < alpha_index_times_count; k++) {
    const int a0 = alpha_index_times[k * 4 + 0];
    const int a1 = alpha_index_times[k * 4 + 1];
    const int mult = alpha_index_times[k * 4 + 2];
    const int a3 = alpha_index_times[k * 4 + 3];
    moment_tensor_vals_[a3] += mult * moment_tensor_vals_[a0] * moment_tensor_vals_[a1];
  }
}

double EMTP::peratom_energyforce_soa(const double* drx, const double* dry, const double* drz,
                                      int ti_0indexed, const int* tj_0indexed, int Nj,
                                      const int* type_map)
{
  const int ti0 = type_map[ti_0indexed];
  for (int j = 0; j < Nj; j++) soa_tj_[j] = type_map[tj_0indexed[j]];

  build_moments_and_jacobian(drx, dry, drz, ti0, soa_tj_.data(), Nj);

  double energy = species_coeffs[ti0];
  std::fill(nbh_energy_ders_wrt_moments_.begin(), nbh_energy_ders_wrt_moments_.end(), 0.0);
  for (int k = 0; k < alpha_scalar_moments; k++) {
    const int m = alpha_moment_mapping[k];
    energy += linear_coeffs[k] * moment_tensor_vals_[m];
    nbh_energy_ders_wrt_moments_[m] = linear_coeffs[k];
  }

  // Reverse-mode backprop through the contraction DAG (pair_mtp.cpp:221-233) -- walk
  // alpha_index_times BACKWARD, product rule. Both += lines execute even when a0==a1 (no
  // special-casing needed: val3 is read before either target is written).
  for (int k = alpha_index_times_count - 1; k >= 0; k--) {
    const int a0 = alpha_index_times[k * 4 + 0];
    const int a1 = alpha_index_times[k * 4 + 1];
    const int mult = alpha_index_times[k * 4 + 2];
    const int a3 = alpha_index_times[k * 4 + 3];
    const double val0 = moment_tensor_vals_[a0];
    const double val1 = moment_tensor_vals_[a1];
    const double val3 = nbh_energy_ders_wrt_moments_[a3];
    nbh_energy_ders_wrt_moments_[a1] += val3 * mult * val0;
    nbh_energy_ders_wrt_moments_[a0] += val3 * mult * val1;
  }

  // Per-neighbor jacobian dot (pair_mtp.cpp:236-254) -- basic moments only (indices
  // [0, alpha_index_basic_count) in moment_tensor_vals_, by the file's own convention).
  for (int jj = 0; jj < Nj; jj++) {
    double fx = 0.0, fy = 0.0, fz = 0.0;
    for (int k = 0; k < alpha_index_basic_count; k++) {
      const size_t jbase = (static_cast<size_t>(jj) * alpha_index_basic_count + k) * 3;
      const double d = nbh_energy_ders_wrt_moments_[k];
      fx += d * moment_jacobian_[jbase + 0];
      fy += d * moment_jacobian_[jbase + 1];
      fz += d * moment_jacobian_[jbase + 2];
    }
    soa_fij[jj * 3 + 0] = fx;
    soa_fij[jj * 3 + 1] = fy;
    soa_fij[jj * 3 + 2] = fz;
  }

  return energy;
}

void EMTP::peratom_descriptors_soa(const double* drx, const double* dry, const double* drz,
                                    int ti_0indexed, const int* tj_0indexed, int Nj,
                                    const int* type_map)
{
  const int ti0 = type_map[ti_0indexed];
  for (int j = 0; j < Nj; j++) soa_tj_[j] = type_map[tj_0indexed[j]];

  build_moments_and_jacobian(drx, dry, drz, ti0, soa_tj_.data(), Nj);

  for (int k = 0; k < alpha_scalar_moments; k++) bd[k] = moment_tensor_vals_[alpha_moment_mapping[k]];

  // Descriptor-mode backprop -- no LAMMPS analogue (LAMMPS only ever backprops a single
  // linear-coefficient-weighted scalar). Seed moment_ders_wrt_basis_[q][k] = d(moment[q])/d(B_k)
  // as an identity at each B_k's own moment slot, then walk the same reverse DAG applied to every
  // column k at once (product rule, same aliasing-safe read-before-write shape as above).
  std::fill(moment_ders_wrt_basis_.begin(), moment_ders_wrt_basis_.end(), 0.0);
  for (int k = 0; k < alpha_scalar_moments; k++)
    moment_ders_wrt_basis_[static_cast<size_t>(alpha_moment_mapping[k]) * alpha_scalar_moments + k] = 1.0;

  for (int t = alpha_index_times_count - 1; t >= 0; t--) {
    const int a0 = alpha_index_times[t * 4 + 0];
    const int a1 = alpha_index_times[t * 4 + 1];
    const int mult = alpha_index_times[t * 4 + 2];
    const int a3 = alpha_index_times[t * 4 + 3];
    const double val0 = moment_tensor_vals_[a0];
    const double val1 = moment_tensor_vals_[a1];
    const double* d3 = &moment_ders_wrt_basis_[static_cast<size_t>(a3) * alpha_scalar_moments];
    double* d1 = &moment_ders_wrt_basis_[static_cast<size_t>(a1) * alpha_scalar_moments];
    double* d0 = &moment_ders_wrt_basis_[static_cast<size_t>(a0) * alpha_scalar_moments];
    for (int k = 0; k < alpha_scalar_moments; k++) {
      const double dv3 = d3[k];
      d1[k] += dv3 * mult * val0;
      d0[k] += dv3 * mult * val1;
    }
  }

  // Per-neighbor, per-basis-function jacobian: dot moment_ders_wrt_basis_ against the *basic*-
  // moment Jacobian only (contracted moments' Jacobians are never needed standalone, exactly as
  // in the energy/force backprop). Addressed with the actual call's Nj as stride (matches
  // EAPOD's own bdd addressing convention), NOT the fixed Njmax buffer capacity.
  for (int jj = 0; jj < Nj; jj++) {
    for (int k = 0; k < alpha_scalar_moments; k++) {
      double vx = 0.0, vy = 0.0, vz = 0.0;
      for (int q = 0; q < alpha_index_basic_count; q++) {
        const double w = moment_ders_wrt_basis_[static_cast<size_t>(q) * alpha_scalar_moments + k];
        if (w == 0.0) continue;
        const size_t jbase = (static_cast<size_t>(jj) * alpha_index_basic_count + q) * 3;
        vx += w * moment_jacobian_[jbase + 0];
        vy += w * moment_jacobian_[jbase + 1];
        vz += w * moment_jacobian_[jbase + 2];
      }
      // Negated: vx/vy/vz above is d(B_k)/d(r), r = r_neighbor - r_central (the RELATIVE vector
      // moment_jacobian_ is defined against, same as the force path's temp_force). For a genuine
      // (non-force) derivative there is no extra "F=-dE/dr" sign flip, so d(B_k)/d(r_central_abs)
      // = -d(B_k)/d(r) and d(B_k)/d(r_neighbor_abs) = +d(B_k)/d(r) -- the OPPOSITE relationship
      // from the force case. Storing bdd pre-negated here lets mtp_descriptor_op.h use the same
      // central+=/neighbor-= scatter idiom as everywhere else in this codebase (matches k2b's own
      // documented sign-convention note in k2b_descriptor_op.h; verified against a direct
      // multi-atom finite-difference check on the total system descriptor sum -- the naive
      // unnegated sign failed that check by O(0.1), the negated one passes to ~1e-9).
      const size_t out = 3 * static_cast<size_t>(jj) + 3 * static_cast<size_t>(Nj) * k;
      bdd[out + 0] = -vx;
      bdd[out + 1] = -vy;
      bdd[out + 2] = -vz;
    }
  }
}
