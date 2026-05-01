/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <yaml-cpp/yaml.h>

#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <onika/file_utils.h>
#include <onika/log.h>
#include <onika/memory/allocator.h>
#include <onika/cuda/ro_shallow_copy.h>
#include <exaStamp/potential_factory/pair_potential.h>
#include <exaStamp/unit_system.h>

#include <onika/cuda/cuda.h>

namespace exaStamp
{
  using namespace exanb;

  struct alignas(sizeof(double)*8) SplineCoeffs
  {
    static constexpr size_t N_SPLINE_POINTS_STORAGE = 8;
    alignas(sizeof(double)*8) double coeffs[N_SPLINE_POINTS_STORAGE] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  };
  static_assert( sizeof(SplineCoeffs) == ( SplineCoeffs::N_SPLINE_POINTS_STORAGE * sizeof(double) ) );

  // ===========================================================================
  // AireboParams  — full parameter set, lives on the CPU heap (inside a slot).
  // Plain C arrays store data for CPU-side access and file reading.
  // CudaMMVector mirrors hold the same data in GPU-accessible unified memory;
  // they are filled by finalize_spline_data() after file loading.
  // ===========================================================================

  struct AireboParams
  {
    std::string potential_file; // path to CH.airebo or CH.airebo-m
    double cutlj     = 3.0;    // LJ cutoff multiplier (in units of sigma)
    double sigmin    = 2.16;   // inner LJ switching lower bound (in units of sigma); default = cutlj-0.84
    bool   ljflag    = true;   // include Lennard-Jones term
    bool   torflag   = true;   // include torsion term
    bool   morseflag = false;  // use Morse instead of LJ for non-bonded

    // atom type indices: C=0, H=1 (used for [i][j] pair indexing below)

    // ===== REBO scalar parameters =====
    double smin  = 0.0, Nmin  = 0.0, Nmax  = 0.0;
    double NCmin = 0.0, NCmax = 0.0;

    double rcmin[2][2]   = {};   // minimum REBO cutoff
    double rcmax[2][2]   = {};   // maximum REBO cutoff
    double rcmaxsq[2][2] = {};   // rcmax squared (precomputed)
    double rcmaxp[2][2]  = {};   // rcmax for the pi-bond cutoff

    double Q[2][2]       = {};   // Brenner Q parameter
    double alpha[2][2]   = {};   // exponential decay
    double A[2][2]       = {};   // repulsive pre-factor
    double rho[2][2]     = {};   // Morse-like length scale
    double BIJc[2][2][3] = {};   // attractive pre-factors   [i][j][term 0..2]
    double Beta[2][2][3] = {};   // attractive decay lengths [i][j][term 0..2]

    // ===== LJ / torsion scalar parameters =====
    double rcLJmin[2][2]   = {};
    double rcLJmax[2][2]   = {};
    double rcLJmaxsq[2][2] = {};
    double bLJmin[2][2]    = {};
    double bLJmax[2][2]    = {};
    double epsilon[2][2]   = {};  // LJ well depth
    double sigma[2][2]     = {};  // LJ radius
    double epsilonT[2][2]  = {};  // torsion energy scale (CC-CC, CC-CH, HC-CH)

    // ===== 1D quintic splines (CPU-side plain arrays) =====
    double gCdom[5]  = {};
    double gC1[4][6] = {};
    double gC2[4][6] = {};
    double gHdom[4]  = {};
    double gH[3][6]  = {};

    // ===== 2D bicubic splines (CPU-side) =====
    double pCCdom[2][2]  = {};
    double pCC[4][4][16] = {};
    double pCHdom[2][2]  = {};
    double pCH[4][4][16] = {};

    // ===== 3D tricubic splines (CPU-side) =====
    double piCCdom[3][2]     = {};
    double piCC[4][4][9][64] = {};
    double piCHdom[3][2]     = {};
    double piCH[4][4][9][64] = {};
    double piHHdom[3][2]     = {};
    double piHH[4][4][9][64] = {};
    double Tijdom[3][2]      = {};
    double Tijc[4][4][9][64] = {};

    // ===== GPU-accessible flat mirrors (filled by finalize_spline_data) =====
    onika::memory::CudaMMVector<double> gCdom_data;
    onika::memory::CudaMMVector<double> gC1_data;
    onika::memory::CudaMMVector<double> gC2_data;
    onika::memory::CudaMMVector<double> gHdom_data;
    onika::memory::CudaMMVector<double> gH_data;
    onika::memory::CudaMMVector<double> pCCdom_data;
    onika::memory::CudaMMVector<double> pCC_data;
    onika::memory::CudaMMVector<double> pCHdom_data;
    onika::memory::CudaMMVector<double> pCH_data;
    onika::memory::CudaMMVector<double> piCCdom_data;
    onika::memory::CudaMMVector<double> piCC_data;
    onika::memory::CudaMMVector<double> piCHdom_data;
    onika::memory::CudaMMVector<double> piCH_data;
    onika::memory::CudaMMVector<double> piHHdom_data;
    onika::memory::CudaMMVector<double> piHH_data;
    onika::memory::CudaMMVector<double> Tijdom_data;
    onika::memory::CudaMMVector<double> Tijc_data;

    // Copy all spline arrays into their GPU-accessible CudaMMVector mirrors.
    // Must be called once after the parameter file has been fully read.
    inline void finalize_spline_data()
    {
      auto fill = [](onika::memory::CudaMMVector<double>& dst, const double* src, size_t n) {
        dst.resize(n);
        std::copy(src, src + n, dst.data());
      };
      fill(gCdom_data,  gCdom,            5      );
      fill(gC1_data,    &gC1[0][0],       4*6    );
      fill(gC2_data,    &gC2[0][0],       4*6    );
      fill(gHdom_data,  gHdom,            4      );
      fill(gH_data,     &gH[0][0],        3*6    );
      fill(pCCdom_data, &pCCdom[0][0],    4      );
      fill(pCC_data,    &pCC[0][0][0],    4*4*16 );
      fill(pCHdom_data, &pCHdom[0][0],    4      );
      fill(pCH_data,    &pCH[0][0][0],    4*4*16 );
      fill(piCCdom_data,&piCCdom[0][0],   6      );
      fill(piCC_data,   &piCC[0][0][0][0],4*4*9*64);
      fill(piCHdom_data,&piCHdom[0][0],   6      );
      fill(piCH_data,   &piCH[0][0][0][0],4*4*9*64);
      fill(piHHdom_data,&piHHdom[0][0],   6      );
      fill(piHH_data,   &piHH[0][0][0][0],4*4*9*64);
      fill(Tijdom_data, &Tijdom[0][0],    6      );
      fill(Tijc_data,   &Tijc[0][0][0][0],4*4*9*64);
    }
  };

  // ===========================================================================
  // AireboParamsRO  — lightweight read-only view passed to compute functors.
  // Scalars are copied inline; spline arrays are accessed via raw pointers into
  // the CudaMMVector buffers of AireboParams.  Designed to fit within the
  // 1024-byte functor size limit enforced by onika::parallel::block_parallel_for.
  // ===========================================================================

  struct AireboParamsRO
  {
    // --- configuration flags ---
    bool   ljflag    = true;
    bool   torflag   = true;
    bool   morseflag = false;
    double cutlj     = 3.0;
    double sigmin    = 2.16;

    // --- REBO scalar parameters (copied inline) ---
    double smin=0, Nmin=0, Nmax=0, NCmin=0, NCmax=0;
    double rcmin[2][2]   = {};
    double rcmax[2][2]   = {};
    double rcmaxsq[2][2] = {};
    double rcmaxp[2][2]  = {};
    double Q[2][2]       = {};
    double alpha[2][2]   = {};
    double A[2][2]       = {};
    double rho[2][2]     = {};
    double BIJc[2][2][3] = {};
    double Beta[2][2][3] = {};

    // --- LJ / torsion scalar parameters (copied inline) ---
    double rcLJmin[2][2]   = {};
    double rcLJmax[2][2]   = {};
    double rcLJmaxsq[2][2] = {};
    double bLJmin[2][2]    = {};
    double bLJmax[2][2]    = {};
    double epsilon[2][2]   = {};
    double sigma[2][2]     = {};
    double epsilonT[2][2]  = {};

    // --- spline data: raw pointers into GPU-accessible CudaMMVector buffers ---
    const double* __restrict__ gCdom   = nullptr;   // [5]
    const double* __restrict__ gC1     = nullptr;   // [4*6]
    const double* __restrict__ gC2     = nullptr;   // [4*6]
    const double* __restrict__ gHdom   = nullptr;   // [4]
    const double* __restrict__ gH      = nullptr;   // [3*6]
    const double* __restrict__ pCCdom  = nullptr;   // [2*2]
    const double* __restrict__ pCC     = nullptr;   // [4*4*16]
    const double* __restrict__ pCHdom  = nullptr;   // [2*2]
    const double* __restrict__ pCH     = nullptr;   // [4*4*16]
    const double* __restrict__ piCCdom = nullptr;   // [3*2]
    const double* __restrict__ piCC    = nullptr;   // [4*4*9*64]
    const double* __restrict__ piCHdom = nullptr;   // [3*2]
    const double* __restrict__ piCH    = nullptr;   // [4*4*9*64]
    const double* __restrict__ piHHdom = nullptr;   // [3*2]
    const double* __restrict__ piHH    = nullptr;   // [4*4*9*64]
    const double* __restrict__ Tijdom  = nullptr;   // [3*2]
    const double* __restrict__ Tijc    = nullptr;   // [4*4*9*64]

    AireboParamsRO() = default;
    AireboParamsRO(const AireboParamsRO&) = default;

    inline AireboParamsRO(const AireboParams& p)
      : ljflag(p.ljflag), torflag(p.torflag), morseflag(p.morseflag)
      , cutlj(p.cutlj), sigmin(p.sigmin)
      , smin(p.smin), Nmin(p.Nmin), Nmax(p.Nmax), NCmin(p.NCmin), NCmax(p.NCmax)
      , gCdom  (p.gCdom_data.data())
      , gC1    (p.gC1_data.data())
      , gC2    (p.gC2_data.data())
      , gHdom  (p.gHdom_data.data())
      , gH     (p.gH_data.data())
      , pCCdom (p.pCCdom_data.data())
      , pCC    (p.pCC_data.data())
      , pCHdom (p.pCHdom_data.data())
      , pCH    (p.pCH_data.data())
      , piCCdom(p.piCCdom_data.data())
      , piCC   (p.piCC_data.data())
      , piCHdom(p.piCHdom_data.data())
      , piCH   (p.piCH_data.data())
      , piHHdom(p.piHHdom_data.data())
      , piHH   (p.piHH_data.data())
      , Tijdom (p.Tijdom_data.data())
      , Tijc   (p.Tijc_data.data())
    {
      // copy 2D scalar arrays
#define AIREBO_COPY2(dst, src) \
      dst[0][0]=src[0][0]; dst[0][1]=src[0][1]; dst[1][0]=src[1][0]; dst[1][1]=src[1][1];
      AIREBO_COPY2(rcmin,   p.rcmin)
      AIREBO_COPY2(rcmax,   p.rcmax)
      AIREBO_COPY2(rcmaxsq, p.rcmaxsq)
      AIREBO_COPY2(rcmaxp,  p.rcmaxp)
      AIREBO_COPY2(Q,       p.Q)
      AIREBO_COPY2(alpha,   p.alpha)
      AIREBO_COPY2(A,       p.A)
      AIREBO_COPY2(rho,     p.rho)
      AIREBO_COPY2(rcLJmin,   p.rcLJmin)
      AIREBO_COPY2(rcLJmax,   p.rcLJmax)
      AIREBO_COPY2(rcLJmaxsq, p.rcLJmaxsq)
      AIREBO_COPY2(bLJmin,    p.bLJmin)
      AIREBO_COPY2(bLJmax,    p.bLJmax)
      AIREBO_COPY2(epsilon,   p.epsilon)
      AIREBO_COPY2(sigma,     p.sigma)
      AIREBO_COPY2(epsilonT,  p.epsilonT)
#undef AIREBO_COPY2
      // copy BIJc and Beta [2][2][3]
      for(int i=0;i<2;i++) for(int j=0;j<2;j++) for(int k=0;k<3;k++) {
        BIJc[i][j][k] = p.BIJc[i][j][k];
        Beta[i][j][k] = p.Beta[i][j][k];
      }
    }
  };

} // namespace exaStamp

namespace onika
{
  namespace cuda
  {
    template<> struct ReadOnlyShallowCopyType<exaStamp::AireboParams>
    {
      using type = exaStamp::AireboParamsRO;
    };
  }
}

namespace YAML
{
  using exaStamp::AireboParams;
  using onika::physics::Quantity;

  template<> struct convert<AireboParams>
  {
    static bool decode(const Node& node, AireboParams& v)
    {
      using exanb::lout;
      using exanb::lerr;

      if (!node.IsMap())           { return false; }
      if (!node["potential_file"]) { return false; }

      v.potential_file = node["potential_file"].as<std::string>();
      if (node["cutlj"])     v.cutlj     = node["cutlj"].as<Quantity>().convert();
      v.sigmin = v.cutlj - 0.84;  // default: LAMMPS sigwid = 0.84
      if (node["sigmin"])    v.sigmin    = node["sigmin"].as<Quantity>().convert();
      if (node["ljflag"])    v.ljflag    = node["ljflag"].as<bool>();
      if (node["torflag"])   v.torflag   = node["torflag"].as<bool>();
      if (node["morseflag"]) v.morseflag = node["morseflag"].as<bool>();

      const std::string file_path = onika::data_file_path(v.potential_file);
      std::ifstream f(file_path);
      if (!f.is_open()) {
        lerr << "Cannot open AIREBO potential file: " << file_path << std::endl;
        return false;
      }
      lout << "Reading AIREBO potential file: " << file_path << std::endl;

      // Read the next double from the file, skipping blank lines and comment
      // lines (lines whose first non-whitespace character is '#' or that contain
      // no parseable number).  Scalar-param lines have the format "value  label";
      // only the first token is consumed.
      auto next_double = [&]() -> double {
        std::string line;
        while (std::getline(f, line)) {
          std::istringstream iss(line);
          double val;
          if (iss >> val) return val;
        }
        lerr << "AIREBO file: unexpected end of file while reading " << file_path << std::endl;
        return 0.0;
      };
      auto next_int = [&]() -> int { return static_cast<int>(next_double()); };

      // Read N consecutive doubles into a flat C array
      auto read_n = [&](double* dst, int n) {
        for (int i = 0; i < n; ++i) dst[i] = next_double();
      };

      // ======================================================================
      // Scalar parameters — file order matches the CH.airebo / CH.airebo-m
      // parameter list (same order as LAMMPS PairAIREBO::read_file).
      // ======================================================================

      // REBO cutoffs
      double rcmin_CC  = next_double(), rcmin_CH  = next_double(), rcmin_HH  = next_double();
      double rcmax_CC  = next_double(), rcmax_CH  = next_double(), rcmax_HH  = next_double();
      double rcmaxp_CC = next_double(), rcmaxp_CH = next_double(), rcmaxp_HH = next_double();

      v.smin  = next_double(); v.Nmin  = next_double(); v.Nmax  = next_double();
      v.NCmin = next_double(); v.NCmax = next_double();

      // REBO pair parameters
      double Q_CC    = next_double(), Q_CH    = next_double(), Q_HH    = next_double();
      double alpha_CC = next_double(), alpha_CH = next_double(), alpha_HH = next_double();
      double A_CC    = next_double(), A_CH    = next_double(), A_HH    = next_double();

      double BIJc_CC1 = next_double(), BIJc_CC2 = next_double(), BIJc_CC3 = next_double();
      double BIJc_CH1 = next_double(), BIJc_CH2 = next_double(), BIJc_CH3 = next_double();
      double BIJc_HH1 = next_double(), BIJc_HH2 = next_double(), BIJc_HH3 = next_double();

      double Beta_CC1 = next_double(), Beta_CC2 = next_double(), Beta_CC3 = next_double();
      double Beta_CH1 = next_double(), Beta_CH2 = next_double(), Beta_CH3 = next_double();
      double Beta_HH1 = next_double(), Beta_HH2 = next_double(), Beta_HH3 = next_double();

      double rho_CC = next_double(), rho_CH = next_double(), rho_HH = next_double();

      // LJ parameters
      double rcLJmin_CC = next_double(), rcLJmin_CH = next_double(), rcLJmin_HH = next_double();
      double rcLJmax_CC = next_double(), rcLJmax_CH = next_double(), rcLJmax_HH = next_double();
      double bLJmin_CC  = next_double(), bLJmin_CH  = next_double(), bLJmin_HH  = next_double();
      double bLJmax_CC  = next_double(), bLJmax_CH  = next_double(), bLJmax_HH  = next_double();
      double epsilon_CC = next_double(), epsilon_CH = next_double(), epsilon_HH = next_double();
      double sigma_CC   = next_double(), sigma_CH   = next_double(), sigma_HH   = next_double();

      // Torsion energy scales
      double epsilonT_CCCC = next_double(), epsilonT_CCCH = next_double(), epsilonT_HCCH = next_double();

      // ------------------------------------------------------------------
      // Assign into 2D arrays with C=0, H=1 and [i][j]=[j][i] symmetry
      // ------------------------------------------------------------------

      v.rcmin[0][0]  = rcmin_CC;   v.rcmin[0][1]  = rcmin_CH;
      v.rcmin[1][0]  = rcmin_CH;   v.rcmin[1][1]  = rcmin_HH;

      v.rcmax[0][0]  = rcmax_CC;   v.rcmax[0][1]  = rcmax_CH;
      v.rcmax[1][0]  = rcmax_CH;   v.rcmax[1][1]  = rcmax_HH;

      v.rcmaxsq[0][0] = rcmax_CC*rcmax_CC;  v.rcmaxsq[0][1] = rcmax_CH*rcmax_CH;
      v.rcmaxsq[1][0] = rcmax_CH*rcmax_CH;  v.rcmaxsq[1][1] = rcmax_HH*rcmax_HH;

      v.rcmaxp[0][0] = rcmaxp_CC;  v.rcmaxp[0][1] = rcmaxp_CH;
      v.rcmaxp[1][0] = rcmaxp_CH;  v.rcmaxp[1][1] = rcmaxp_HH;

      v.Q[0][0]     = Q_CC;      v.Q[0][1]     = Q_CH;
      v.Q[1][0]     = Q_CH;      v.Q[1][1]     = Q_HH;

      v.alpha[0][0] = alpha_CC;  v.alpha[0][1] = alpha_CH;
      v.alpha[1][0] = alpha_CH;  v.alpha[1][1] = alpha_HH;

      v.A[0][0]     = A_CC;      v.A[0][1]     = A_CH;
      v.A[1][0]     = A_CH;      v.A[1][1]     = A_HH;

      v.rho[0][0]   = rho_CC;    v.rho[0][1]   = rho_CH;
      v.rho[1][0]   = rho_CH;    v.rho[1][1]   = rho_HH;

      v.BIJc[0][0][0] = BIJc_CC1; v.BIJc[0][0][1] = BIJc_CC2; v.BIJc[0][0][2] = BIJc_CC3;
      v.BIJc[0][1][0] = BIJc_CH1; v.BIJc[0][1][1] = BIJc_CH2; v.BIJc[0][1][2] = BIJc_CH3;
      v.BIJc[1][0][0] = BIJc_CH1; v.BIJc[1][0][1] = BIJc_CH2; v.BIJc[1][0][2] = BIJc_CH3;
      v.BIJc[1][1][0] = BIJc_HH1; v.BIJc[1][1][1] = BIJc_HH2; v.BIJc[1][1][2] = BIJc_HH3;

      v.Beta[0][0][0] = Beta_CC1; v.Beta[0][0][1] = Beta_CC2; v.Beta[0][0][2] = Beta_CC3;
      v.Beta[0][1][0] = Beta_CH1; v.Beta[0][1][1] = Beta_CH2; v.Beta[0][1][2] = Beta_CH3;
      v.Beta[1][0][0] = Beta_CH1; v.Beta[1][0][1] = Beta_CH2; v.Beta[1][0][2] = Beta_CH3;
      v.Beta[1][1][0] = Beta_HH1; v.Beta[1][1][1] = Beta_HH2; v.Beta[1][1][2] = Beta_HH3;

      v.rcLJmin[0][0] = rcLJmin_CC;  v.rcLJmin[0][1] = rcLJmin_CH;
      v.rcLJmin[1][0] = rcLJmin_CH;  v.rcLJmin[1][1] = rcLJmin_HH;

      v.rcLJmax[0][0] = rcLJmax_CC;  v.rcLJmax[0][1] = rcLJmax_CH;
      v.rcLJmax[1][0] = rcLJmax_CH;  v.rcLJmax[1][1] = rcLJmax_HH;

      v.rcLJmaxsq[0][0] = rcLJmax_CC*rcLJmax_CC;  v.rcLJmaxsq[0][1] = rcLJmax_CH*rcLJmax_CH;
      v.rcLJmaxsq[1][0] = rcLJmax_CH*rcLJmax_CH;  v.rcLJmaxsq[1][1] = rcLJmax_HH*rcLJmax_HH;

      v.bLJmin[0][0] = bLJmin_CC;  v.bLJmin[0][1] = bLJmin_CH;
      v.bLJmin[1][0] = bLJmin_CH;  v.bLJmin[1][1] = bLJmin_HH;

      v.bLJmax[0][0] = bLJmax_CC;  v.bLJmax[0][1] = bLJmax_CH;
      v.bLJmax[1][0] = bLJmax_CH;  v.bLJmax[1][1] = bLJmax_HH;

      v.epsilon[0][0] = epsilon_CC;  v.epsilon[0][1] = epsilon_CH;
      v.epsilon[1][0] = epsilon_CH;  v.epsilon[1][1] = epsilon_HH;

      v.sigma[0][0] = sigma_CC;  v.sigma[0][1] = sigma_CH;
      v.sigma[1][0] = sigma_CH;  v.sigma[1][1] = sigma_HH;

      v.epsilonT[0][0] = epsilonT_CCCC;  v.epsilonT[0][1] = epsilonT_CCCH;
      v.epsilonT[1][0] = epsilonT_CCCH;  v.epsilonT[1][1] = epsilonT_HCCH;

      // ======================================================================
      // Spline blocks — each block is preceded by a "# name" comment line
      // (consumed transparently by next_double).
      // ======================================================================

      // --- gC : carbon angular function, two variants (low/high N_C) ---
      // limit = number of knot positions (5 → 4 intervals)
      {
        int limit = next_int();
        read_n(v.gCdom, limit);
        for (int i = 0; i < limit-1; i++) read_n(v.gC1[i], 6);
        for (int i = 0; i < limit-1; i++) read_n(v.gC2[i], 6);
      }

      // --- gH : hydrogen angular function ---
      // limit = number of knot positions (4 → 3 intervals)
      {
        int limit = next_int();
        read_n(v.gHdom, limit);
        for (int i = 0; i < limit-1; i++) read_n(v.gH[i], 6);
      }

      // --- pCC : 2D bicubic bond-order correction for C-C ---
      // limit/2 × limit/2 domain values, then pCCdom[0][1] × pCCdom[1][1] patches of 16
      {
        int limit = next_int();
        for (int i = 0; i < limit/2; i++) read_n(v.pCCdom[i], limit/2);
        for (int i = 0; i < (int)v.pCCdom[0][1]; i++)
          for (int j = 0; j < (int)v.pCCdom[1][1]; j++)
            read_n(v.pCC[i][j], 16);
      }

      // --- pCH : 2D bicubic bond-order correction for C-H ---
      {
        int limit = next_int();
        for (int i = 0; i < limit/2; i++)
          for (int j = 0; j < limit/2; j++)
            v.pCHdom[i][j] = next_double();
        for (int i = 0; i < (int)v.pCHdom[0][1]; i++)
          for (int j = 0; j < (int)v.pCHdom[1][1]; j++)
            read_n(v.pCH[i][j], 16);
      }

      // --- piCC : 3D tricubic pi-bond correction for C-C ---
      // limit/2 × limit/3 domain values, then piCCdom[0][1] × [1][1] × [2][1] patches of 64
      {
        int limit = next_int();
        for (int i = 0; i < limit/2; i++)
          for (int j = 0; j < limit/3; j++)
            v.piCCdom[i][j] = next_double();
        for (int i = 0; i < (int)v.piCCdom[0][1]; i++)
          for (int j = 0; j < (int)v.piCCdom[1][1]; j++)
            for (int k = 0; k < (int)v.piCCdom[2][1]; k++)
              read_n(v.piCC[i][j][k], 64);
      }

      // --- piCH : 3D tricubic pi-bond correction for C-H ---
      {
        int limit = next_int();
        for (int i = 0; i < limit/2; i++)
          for (int j = 0; j < limit/3; j++)
            v.piCHdom[i][j] = next_double();
        for (int i = 0; i < (int)v.piCHdom[0][1]; i++)
          for (int j = 0; j < (int)v.piCHdom[1][1]; j++)
            for (int k = 0; k < (int)v.piCHdom[2][1]; k++)
              read_n(v.piCH[i][j][k], 64);
      }

      // --- piHH : 3D tricubic pi-bond correction for H-H ---
      {
        int limit = next_int();
        for (int i = 0; i < limit/2; i++)
          for (int j = 0; j < limit/3; j++)
            v.piHHdom[i][j] = next_double();
        for (int i = 0; i < (int)v.piHHdom[0][1]; i++)
          for (int j = 0; j < (int)v.piHHdom[1][1]; j++)
            for (int k = 0; k < (int)v.piHHdom[2][1]; k++)
              read_n(v.piHH[i][j][k], 64);
      }

      // --- Tij : 3D tricubic torsion correction ---
      {
        int limit = next_int();
        for (int i = 0; i < limit/2; i++)
          for (int j = 0; j < limit/3; j++)
            v.Tijdom[i][j] = next_double();
        for (int i = 0; i < (int)v.Tijdom[0][1]; i++)
          for (int j = 0; j < (int)v.Tijdom[1][1]; j++)
            for (int k = 0; k < (int)v.Tijdom[2][1]; k++)
              read_n(v.Tijc[i][j][k], 64);
      }

      // Populate GPU-accessible flat mirrors
      v.finalize_spline_data();

      lout << "AIREBO potential file read successfully." << std::endl;
      return true;
    }
  };
}
