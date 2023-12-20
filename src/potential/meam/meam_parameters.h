#pragma once

#include <onika/cuda/cuda.h>

#ifndef XSTAMP_MEAM_MAX_NEIGHBORS
#define XSTAMP_MEAM_MAX_NEIGHBORS 32
#endif

/// @brief Structure to handle parameters of a MEAM potential
struct MeamParameters
{
  static inline constexpr size_t MAX_PARTICLE_NEIGHBORS = XSTAMP_MEAM_MAX_NEIGHBORS;

	/// @brief Default constructors
  MeamParameters() = default;
  MeamParameters(const MeamParameters&) = default;
  MeamParameters(MeamParameters &&) = default;

  MeamParameters& operator = (const MeamParameters&) = default;
  MeamParameters& operator = (MeamParameters &&) = default;

  /// @brief Constructor from values
  /// @param [in] rmax_ Maximal distance for density contribution
  /// @param [in] rmin_ Minimal distance for density contribution
  /// @param [in] Ecoh_ Cohesive energy
  /// @param [in] E0_ Energy
  /// @param [in] A_ Empirical parameter A
  /// @param [in] r0_ Empirical parameter r0
  /// @param [in] alpha_ Parameter alpha
  /// @param [in] delta_ Parameter delta
  /// @param [in] beta0_ Empirical parameter r0
  /// @param [in] beta1_ Attenuation rate
  /// @param [in] beta2_ Attenuation rate
  /// @param [in] beta3_ Attenuation rate
  /// @param [in] t0_ Unnamed coefficient t0
  /// @param [in] t1_ Unnamed coefficient t1
  /// @param [in] t2_ Unnamed coefficient t2
  /// @param [in] t3_ Unnamed coefficient t3
  /// @param [in] s0_ Unnamed coefficient s0
  /// @param [in] s1_ Unnamed coefficient s1
  /// @param [in] s2_ Unnamed coefficient s2
  /// @param [in] s3_ Unnamed coefficient s3
  /// @param [in] Cmin_ Boundary inf used by screening function
  /// @param [in] Cmax_ Boundary sup used by screening function
  /// @param [in] Z_ Number of neighbors in the reference structure
  /// @param [in] rc_ Distance f(rc)=0
  /// @param [in] rp_ Interval [rc-rp,rc] where we apply a polynome such that P(rc)=0
  ONIKA_HOST_DEVICE_FUNC inline 
  MeamParameters(double rmax_, double rmin_, double Ecoh_, double E0_, double A_, double r0_, double alpha_, double delta_,
		  double beta0_ , double beta1_, double beta2_, double beta3_,
		  double t0_ , double t1_ , double t2_ , double t3_ ,
		  double s0_ , double s1_ , double s2_ , double s3_ ,
		  double Cmin_ , double Cmax_ , double Z_, double rc_, double rp_)
    : rmax(rmax_)
    , rmin(rmin_)
    , Ecoh(Ecoh_)
    , E0(E0_)
    , A(A_)
    , r0(1./r0_)
    , alpha(alpha_)
    , delta(delta_) 
    , beta0(beta0_)
    , beta1(beta1_)
    , beta2(beta2_)
    , beta3(beta3_)
    , t0(t0_)
    , t1(t1_)
    , t2(t2_)
    , t3(t3_)
    , s0(s0_)
    , s1(s1_)
    , s2(s2_)
    , s3(s3_)
    , Cmin(Cmin_)
    , Cmax(Cmax_)
    , Z(Z_)
    , rc(rc_)
    , rp(rp_) {}

  double rmax=0.0;  ///< Maximal distance for density contribution
  double rmin=0.0;  ///< Minimal distance for density contribution
  double Ecoh=0.0;  ///< Cohesive energy
  double E0=0.0;    ///< Energy
  double A=0.0;     ///< Empirical parameter A
  double r0=0.0;    ///< Empirical parameter r0
  double alpha=0.0; ///< Parameter alpha
  double delta=0.0; ///< Parameter delta
  double beta0=0.0; ///< Attenuation rate
  double beta1=0.0; ///< Attenuation rate
  double beta2=0.0; ///< Attenuation rate
  double beta3=0.0; ///< Attenuation rate
  double t0=0.0; 		///< Unnamed coefficient t0
  double t1=0.0; 		///< Unnamed coefficient t1
  double t2=0.0; 		///< Unnamed coefficient t2
  double t3=0.0; 		///< Unnamed coefficient t3
  double s0=0.0; 		///< Unnamed coefficient s0
  double s1=0.0; 		///< Unnamed coefficient s1
  double s2=0.0; 		///< Unnamed coefficient s2
  double s3=0.0; 		///< Unnamed coefficient s3
  double Cmin=0.0;  ///< Boundary inf used by screening function
  double Cmax=0.0;  ///< Boundary sup used by screening function
  double Z=0.0;     ///< Number of neighbors in the reference structure
  double rc=0.0;    ///< Distance f(rc)=0
  double rp=0.0;    ///< Interval [rc-rp,rc] where we apply a polynome such that P(rc)=0
  // -------------
  double rd=0.0;
//  double m[5][4];   ///< Filled with polynomial coefficients of Emu,rho_i between [rc-rp,rc]
  double Rcut=0.0;
  double PolynomeMeam[5*4] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  
};

