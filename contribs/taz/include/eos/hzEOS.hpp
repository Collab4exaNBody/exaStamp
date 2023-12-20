///@file
///@ brief HZ equation of state

#ifndef __HZEOS_HPP_INCLUDED
#define __HZEOS_HPP_INCLUDED

#include "eos/eos.hpp"
#include "simd/hzEOS.hpp"

/// @class HZEOS
/// @brief Class for HZ equation of state
class HZEOS : public IEOS{
public:
  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that equation of state
  /// @brief Shorcut for the simd version of the equation of state
  typedef simd::kernels::HZ<double> simd_opt_t;

  /// @brief Constructor
  /// @param [in] size_ Mesoparticle size
  /// @param [in] mass_ Mesoparticle mass
  /// @param [in] gamma0_ Gruneisen parameter
  /// @param [in] rho0_ Reference density
  /// @param [in] c0_ Parameter c_0
  /// @param [in] cv_ Heat capacity
  /// @param [in] s_ Parameter s
  HZEOS(const uint& size_, const double& mass_, const double& gamma0_, const double& rho0_, const double& c0_, const double& cv_, const double& s_) :
    size(size_), mass(mass_), gamma0(gamma0_), rho0(rho0_), c02(c0_*c0_), cv(cv_), icv(1./cv_), s(s_) {
    
    deltaT0 = convert(298.13, SI_Units_base::kelvin, Stamp_Units::temperature) - convert(1e5, SI_Units_base::joule/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::energy/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length)/gamma0/rho0*icv;
    
  }
  
  /// @brief Destructor
  virtual ~HZEOS() {}

  /// @brief Get EOS type
  /// @return EOS type
  virtual Type getType() const {
    return Type::HZ;
  }

  /// @brief Get EOS name
  /// @return EOS name
  virtual std::string getName() const {
    return "HZ";
  }

  /// @brief Get mesoparticle mass
  /// @return Mass
  inline double getMass() const {
    return mass;
  }

  /// @brief Get mesoparticle size
  /// @return Mesoparticle size
  inline uint getSize() const {
    return size;
  }

  virtual double initSample(double rho_, double T_);

  /// @brief Get energy from density and entropy
  /// @param [in] rho Density
  /// @param [in] S entropy
  /// @return Energy
  inline double getEnergy(const double& rho, const double& S) {
    double y = rho0/rho;
    return (cv*(exp(S/mass*icv-gamma0*y) - deltaT0*exp(gamma0*(1-y))) + eref(rho))*mass;
  }

  /// @brief Get temperature from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Temperature
  inline double getTemperature(const double& rho, const double& e) {
    return (e/mass-eref(rho))*icv + deltaT0*exp(gamma0*(1-rho0/rho));
  }

  /// @brief Get Inverse temperature from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Inverse temperature \f[ \beta \f]
  inline double getBeta(const double& rho, const double& e) {
    return 1./getTemperature(rho,e)/Stamp_Constant::boltzmann;
  }

  /// @brief Get entropy, temperature and pressure from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Entropy, temperature and pressure
  inline vec3<double> getSTP(const double& rho, const double& e) {
    double egx = exp(gamma0*(1-rho0/rho));
    double dem = e/mass-eref(rho);
    
    double S = cv*mass*(log( dem*icv + deltaT0*egx ) + gamma0*rho0/rho);
    double T = dem*icv + deltaT0*egx;
    double P = pref(rho) + rho0*gamma0*dem;
    return vec3<double>(S,T,P);
  }

  /// @brief Get entropy, temperature and pressure from density and energy, along with their derivatives
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @param [out] drt Derivative of temperature with respect to density
  /// @param [out] det Derivative of temperature with respect to energy
  /// @param [out] drp Derivative of pressure with respect to density
  /// @param [out] dep Derivative of pressure with respect to energy
  /// @return Entropy, temperature and pressure
  vec3<double> getSTPAndDerivatives(const double& rho, const double& e, double& drt, double& det, double& drp, double& dep) {
    double egx = exp(gamma0*(1-rho0/rho));
    double dem = e/mass-eref(rho);
    
    double S = cv*mass*(log( dem*icv + deltaT0*egx ) + gamma0*rho0/rho);
    double T = dem*icv + deltaT0*egx;
    double P = pref(rho) + rho0*gamma0*dem;

    drt = (-pref(rho0)*icv + deltaT0*gamma0*rho0*egx)/rho/rho;
    det = icv/mass;
    drp = dpref(rho) - rho0*gamma0*pref(rho)/rho/rho;
    dep = rho0*gamma0/mass;
    
    return vec3<double>(S,T,P);
  }

  /// @brief Get energy, temperature and pressure from density and entropy
  /// @param [in] rho Density
  /// @param [in] S Entropy
  /// @return Energy, temperature and pressure
  inline vec3<double> getETP(const double& rho, const double& S) {
    double egx = exp(gamma0*(1-rho0/rho));
    double dem = cv*(exp(S/mass*icv-gamma0*rho0/rho) - deltaT0*egx);

    double E = mass*(dem+eref(rho));
    double T = dem*icv + deltaT0*egx;
    double P = pref(rho) + rho0*gamma0*dem;
    return vec3<double>(E,T,P);
  }

  /// @brief Get energy, temperature and pressure from density and entropy, along with their derivatives
  /// @param [in] rho Density
  /// @param [in] S Entropy
  /// @param [out] drt Derivative of temperature with respect to density
  /// @param [out] dst Derivative of temperature with respect to entropy
  /// @param [out] drp Derivative of pressure with respect to density
  /// @param [out] dsp Derivative of pressure with respect to entropy
  /// @return Energy, temperature and pressure
  vec3<double> getETPAndDerivatives(const double& rho, const double& S, double& drt, double& dst, double& drp, double& dsp) {
    double egx = exp(gamma0*(1-rho0/rho));
    double dem = cv*(exp(S/mass*icv-gamma0*rho0/rho) - deltaT0*egx);

    double E = mass*(dem+eref(rho));
    double T = dem*icv + deltaT0*egx;
    double P = pref(rho) + rho0*gamma0*dem;

    drt = gamma0*rho0/rho/rho*exp(S/mass*icv-gamma0*rho0/rho);
    dst = icv/mass*exp(S/mass*icv-gamma0*rho0/rho);
    drp = dpref(rho) + rho0*rho0*gamma0*gamma0*dem/rho/rho;
    dsp = rho0*gamma0/mass*exp(S/mass*icv-gamma0*rho0/rho);

    return vec3<double>(E,T,P);
  }

  /// @brief Set the parameters for the simd version of the equation of state
  /// @param [in,out] opt Simd version of the equation of state
  void setParameters(simd::kernels::HZ<double>& opt) {
    opt.setParameters(mass,gamma0,rho0,c02,cv,s,deltaT0);
  }

private:
  double eref(const double& rho);
  double pref(const double& rho);
  double dpref(const double& rho);

  uint size; ///< Size of the mesoparticle
  double mass; ///< Mass of the mesoparticle
  
  double gamma0; ///< Gruneisen coefficient
  double rho0; ///< Reference density
  double c02; ///< Squared parameter \f[ c_0 \f]
  double cv; ///< Heat capacity
  double icv; ///< Inverse heat capacity
  double s; ///< Parameter \f[ s \f]
  double deltaT0; ///< Temperature shift (\f[ T_0-T_{00} \f])

};



/// @brief Initialize internal energy
/// @param [in] rho_ Density
/// @param [in] T_ Temperature
/// @return Energy
inline double HZEOS::initSample(double rho_, double T_) {

  return mass*(eref(rho_) + cv*(T_-deltaT0*exp(gamma0*(1-rho0/rho_)))); 

}


/// @brief Reference energy in the HZ EOS
/// @param [in] rho Density
inline double HZEOS::eref(const double& rho) {

  double x = 1-rho0/rho;
  double x2 = x*x;
  if (x > 0)
    return 0.5*c02*x2/(1-s*x)*(1+s*x/3-s*(gamma0-s)*x2/6);
  else
    return 0.5*c02*x2/(1-s*x);
  
}


/// @brief Reference pressure in the HZ EOS
/// @param [in] rho Density
inline double HZEOS::pref(const double& rho) {

  double x = 1-rho0/rho;
  double sx = s*x;
  if (x > 0) {
    double x2 = x*x;
    return rho0*c02*x/(1-sx)/(1-sx)*(1-s*gamma0*x2/3+s*s*(gamma0-s)*x2*x/4);
  }
  else
    return rho0*c02*x/(1-sx)/(1-sx)*(1-0.5*sx);
  
}

/// @brief Derivative of the reference pressure in the HZ EOS
/// @param [in] rho Density
inline double HZEOS::dpref(const double& rho) {

  double x = 1-rho0/rho;
  double sx = s*x;
  if (x > 0) {
    double sx2 = sx*x;
    return rho0*rho0*c02/(1-sx)/(1-sx)/(1-sx)*(1+sx-sx2*gamma0+sx2*sx*(4./3*gamma0-s)-0.5*sx2*sx*sx*(gamma0-s))/rho/rho;
  }
  else
    return rho0*rho0*c02*x/(1-sx)/(1-sx)/(1-sx)/rho/rho;
  
}

#endif // __HZEOS_HPP_INCLUDED
