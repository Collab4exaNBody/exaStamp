///@file
///@ brief JWL equation of state

#ifndef __JWLEOS_HPP_INCLUDED
#define __JWLEOS_HPP_INCLUDED

#include "eos/eos.hpp"
#include "simd/jwlEOS.hpp"

/// @class JWLEOS
/// @brief Jones-Wilkins-Lee equation of state
class JWLEOS : public IEOS{
public:
  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that potential
  /// @brief Shorcut for the simd version of the potential
  typedef simd::kernels::JWL<double> simd_opt_t;

  /// @brief Constructor
  /// @param [in] size_ Mesoparticle size
  /// @param [in] mass_ Mesoparticle mass
  /// @param [in] gamma0_ Gruneisen parameter
  /// @param [in] rho0_ Reference density
  /// @param [in] E0_ Reference energy
  /// @param [in] Dcj_ Detonation velocity
  /// @param [in] Pcj_ Pressure at the CJ point
  /// @param [in] Tcj_ Temperature at the CJ point
  /// @param [in] cv_ Heat capacity
  /// @param [in] a_ Parameter \f$ a \f$
  /// @param [in] b_ Parameter \f$ b \f$
  /// @param [in] R1_ Parameter \f$ R_1 \f$
  /// @param [in] R2_ Parameter \f$ R_2 \f$
  JWLEOS(const uint& size_, const double& mass_, const double& gamma0_, const double& rho0_, const double& E0_, const double& Dcj_, const double& Pcj_, const double& Tcj_, const double& cv_, const double& a_, const double& b_, const double& R1_, const double& R2_) :
    size(size_), mass(mass_), gamma0(gamma0_), rho0(rho0_), cv(cv_), icv(1./cv_), a(a_), b(b_), R1(R1_), R2(R2_) {

    double rhocj =  1./(1./rho0 - Pcj_/(rho0*rho0*Dcj_*Dcj_));
    double Ecj   = 0.5*Pcj_*(1./rho0-1./rhocj) + E0_;
    double Pk1cj = a*exp(-R1*rho0/rhocj) + b*exp(-R2*rho0/rhocj);
    cek   = -a/rho0/R1*exp(-R1*rho0/rhocj) - b/rho0/R2*exp(-R2*rho0/rhocj) - (Pcj_-Pk1cj)/rhocj/gamma0 + Ecj;
    K     = (Pcj_ - Pk1cj - cv*gamma0*Tcj_*rhocj) * pow(rho0/rhocj,gamma0+1); 

  }
  
  /// @brief Destructor
  virtual ~JWLEOS() {}

  /// @brief Get EOS type
  /// @return EOS type
  virtual Type getType() const {
    return Type::JWL;
  }

  /// @brief Get EOS name
  /// @return EOS name
  virtual std::string getName() const {
    return "Jones-Wilkins-Lee";
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
    return mass*(ek(rho) + cv*exp(S/mass*icv+gamma0*log(rho)));
  }

  /// @brief Get temperature from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Temperature
  inline double getTemperature(const double& rho, const double& e) {
    return (e/mass-ek(rho))*icv;
  }

  /// @brief Get Inverse temperature from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Inverse temperature \f$ \beta \f$
  inline double getBeta(const double& rho, const double& e) {
    return 1./getTemperature(rho,e)/Stamp_Constant::boltzmann;
  }

  /// @brief Get entropy, temperature and pressure from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Entropy, temperature and pressure
  inline vec3<double> getSTP(const double& rho, const double& e) {
    double dem = e/mass-ek(rho);
    
    double S = cv*mass*(log(dem*icv) - gamma0*log(rho));
    double T = dem*icv;
    double P = pk(rho) + rho*gamma0*dem;
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
    double dem = e/mass-ek(rho);
    
    double S = cv*mass*(log(dem*icv) - gamma0*log(rho));
    double T = dem*icv;
    double P = pk(rho) + rho*gamma0*dem;

    drt = -pk(rho)*icv/rho/rho;
    det = icv/mass;
    drp = dpk(rho) + gamma0*dem - gamma0*pk(rho)/rho;
    dep = rho*gamma0/mass;

    return vec3<double>(S,T,P);
  }

  /// @brief Get energy, temperature and pressure from density and entropy
  /// @param [in] rho Density
  /// @param [in] S Entropy
  /// @return Energy, temperature and pressure
  inline vec3<double> getETP(const double& rho, const double& S) {
    double dem = cv*exp(S/mass*icv+gamma0*log(rho));

    double E = mass*(ek(rho)+dem);
    double T = dem*icv;
    double P = pk(rho) + rho*gamma0*dem;
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
    double dem = cv*exp(S/mass*icv+gamma0*log(rho));

    double E = mass*(ek(rho)+dem);
    double T = dem*icv;
    double P = pk(rho) + rho*gamma0*dem;

    drt = gamma0/rho*T;
    dst = icv/mass*exp(S/mass*icv+gamma0*log(rho));
    drp = dpk(rho) + gamma0*dem + gamma0*gamma0*dem;
    dsp = rho*gamma0/mass*exp(S/mass*icv+gamma0*log(rho));

    return vec3<double>(E,T,P);
  }


  /// @brief Set the parameters for the simd version of the equation of state
  /// @param [in,out] opt Simd version of the equation of state
  void setParameters(simd::kernels::JWL<double>& opt) {
    opt.setParameters(mass,gamma0,rho0,cv,a,b,R1,R2,cek,K);
  }

private:
  double ek(const double& rho);
  double pk(const double& rho);
  double dpk(const double& rho);

  uint size; ///< Size of the mesoparticle
  double mass; ///< Mass of the mesoparticle
  
  double gamma0; ///< Gruneisen coefficient
  double rho0; ///< Reference density
  
  double cv; ///< Heat capacity
  double icv; ///< Inverse heat capacity
  double a; ///< Parameter \f$ a \f$
  double b; ///< Parameter \f$ B \f$
  double R1; ///< Parameter \f$ R_1 \f$
  double R2; ///< Parameter \f$ R_2 \f$

  double cek; ///< Constant \f$ C_{ek} \f$
  double K; ///< Constant \f$ K \f$
  
};



/// @brief Initialize internal energy
/// @param [in] rho_ Density
/// @param [in] T_ Temperature
/// @return Energy
inline double JWLEOS::initSample(double rho_, double T_) {

  return mass*(ek(rho_) + cv*T_);

}

/// @brief Reference energy in the JWL EOS
/// @param [in] rho Density
inline double JWLEOS::ek(const double& rho) {

  double y = rho0/rho;
  
  return a/rho0/R1*exp(-R1*y) + b/rho0/R2*exp(-R2*y) + K/rho0/gamma0*pow(y,-gamma0) + cek;
  
}

/// @brief Reference pressure in the JWL EOS
/// @param [in] rho Density
inline double JWLEOS::pk(const double& rho) {

  double y = rho0/rho;
  
  return a*exp(-R1*y) + b*exp(-R2*y) + K*pow(y,-gamma0-1);
  
}

/// @brief Derivative of the reference pressure in the JWL EOS
/// @param [in] rho Density
inline double JWLEOS::dpk(const double& rho) {

  double y = rho0/rho;
  
  return (a*R1*exp(-R1*y) + b*R2*exp(-R2*y)+ K*(gamma0+1)*pow(y,-gamma0-2))*rho0/rho/rho ; 
  
}

#endif // __JWLEOS_HPP_INCLUDED
