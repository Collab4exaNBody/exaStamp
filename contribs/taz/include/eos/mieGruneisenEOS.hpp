/// @file
/// @brief Definition of MieGruneisenEOS class

#ifndef __MIE_GRUNEISEN_EOS_HPP_INCLUDED
#define __MIE_GRUNEISEN_EOS_HPP_INCLUDED


#include "eos/eos.hpp"

#include "simd/mgEOS.hpp"

#include "utils/stampUnits.hpp"


/// @brief Class for Mie Gruneisen equation of state
/// @see http://www.sciencedirect.com/science/article/pii/S1631072112002082
class MieGruneisenEOS : public IEOS{
public:
  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that potential
  /// @brief Shorcut for the simd version of the potential
  typedef simd::kernels::MieGruneisen<double> simd_opt_t;

  /// @brief Constructor
  /// @param [in] size_ Size of the mesoparticles
  /// @param [in] mass_ Mass of the mesoparticles
  /// @param [in] g0 \f[ \Gamma_0 \f] 
  /// @param [in] ginf \f[ \Gamma_{\infty} \f]  
  /// @param [in] t0 \f[ \theta_0 \f] 
  /// @param [in] q_ \f[ q \f] 
  /// @param [in] d0 \f[ \rho_0 \f]
  /// @param [in] ks_ \f[ K_S \f]
  /// @param [in] ns_ \f[ N_S \f] 
  /// @param [in] ds \f[ \rho_S \f] 
  /// @param [in] ur_ \f[ u_r \f]
  /// @param [in] cvr_ \f[ Cv_r \f]
  MieGruneisenEOS(const uint& size_, const double& mass_, const double& g0, const double& ginf, const double& t0, const double& q_,
		  const double& d0, const double& ks_, const double& ns_, const double& ds, const double& ur_, const double& cvr_) :
    size(size_), mass(mass_), gamma0(g0), gammaInf(ginf), theta0(t0), q(q_), rho0(d0), ks(ks_), npu(ns_+1.), rhos(ds), ur(ur_), cvr(cvr_), icvr(1./cvr_) {}

  /// @brief Destructor
  virtual ~MieGruneisenEOS() {}

  /// @brief Get EOS type
  /// @return EOS type
  virtual Type getType() const {
    return Type::MIE_GRUNEISEN;
  }

  /// @brief Get EOS name
  /// @return EOS name
  virtual std::string getName() const {
    return "Mie-Gruneisen";
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
    return (cvr*ur*theta(rho)*(exp(S/mass*icvr)-1) + pote(rho))*mass;
  }

  /// @brief Get temperature from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Temperature
  inline double getTemperature(const double& rho, const double& e) {
    return redu(rho,e/mass)*theta(rho);
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
    double S = mass*cvr*log(redu(rho,e/mass)/ur);
    double T = redu(rho,e/mass)*theta(rho);
    double P = potp(rho) + rho*gamma(rho)*(e/mass-pote(rho));
    return vec3<double>(S,T,P);
  }

  /// @brief Get energy, temperature and pressure from density and entropy
  /// @param [in] rho Density
  /// @param [in] S Entropy
  /// @return Energy, temperature and pressure
  inline vec3<double> getETP(const double& rho, const double& S) {
    double E = cvr*ur*theta(rho)*(exp(S/mass*icvr)-1) + pote(rho);
    double T = redu(rho,E)*theta(rho);;
    double P = potp(rho) + rho*gamma(rho)*(E-pote(rho));
    E *= mass;
    return vec3<double>(E,T,P);
  }

  /// @brief Set the parameters for the simd version of the equation of state
  /// @param [in,out] opt Simd version of the equation of state
  void setParameters(simd::kernels::MieGruneisen<double>& opt) {
    opt.setParameters(mass,gamma0,gammaInf,theta0,q,rho0,ks,npu,rhos,ur,cvr);
  }

private:
  /// @brief \f[ \Gamma \f] function in the Mie-Gruneisen EoS
  /// \f[ \Gamma(\rho) = \Gamma_{\infty} + (\Gamma_0-\Gamma_{\infty})\left(\frac{\rho_0}{\rho}\right)^q \f]
  /// @param [in] rho Density
  double gamma(const double& rho);
  
  /// @brief \f[ \theta \f] function in the Mie-Gruneisen EoS
  /// \f[ \theta(\rho) = \theta_0\left(\frac{\rho_0}{\rho}\right)^{-\Gamma_{\infty}}\exp\left[\frac{\Gamma_0-\Gamma_{\infty}}{q}\left(1-\left(\frac{\rho_0}{\rho}\right)^q\right)\right] \f]
  /// @param [in] rho Density
  double theta(const double& rho);

  /// @brief Potential function in the Mie-Gruneisen EoS
  /// \f[ E_s(\rho) = \frac{K_s}{\rho_s}\left(\frac{\exp((N_s+1)x)}{(N_s+1)^2} - \frac{x}{N_s+1}\right) \f]
  /// with \f[ x = 1-\frac{\rho_s}{\rho} \f]
  /// @param [in] rho Density
  double pote(const double& rho);

  /// @brief Pressure potential function in the Mie-Gruneisen EoS
  /// \f[ P_s(\rho) = \frac{K_s}{N_s+1}\left(\exp((N_s+1)x) -1 \right) \f]
  /// with \f[ x = 1-\frac{\rho_s}{\rho} \f]
  /// @param [in] rho Density
  double potp(const double& rho);

  /// @brief Reduced energy function in the Mie-Gruneisen EoS
  /// \f[ u(e,\rho) = \frac{e-E_s(\rho)}{Cv_r\theta(\rho)}+u_r \f]
  /// @param [in] rho Density
  /// @param [in] e Energy per mass
  double redu(const double& rho, const double& e);
  
  uint size; ///< Size of the mesoparticle
  double mass; ///< Mass of the mesoparticle

  double gamma0; ///< Parameter \f[ \Gamma_0 \f] of the Mie-Gruneisen EoS
  double gammaInf; ///< Parameter \f[ \Gamma_{\infty} \f] of the Mie-Gruneisen EoS
  double theta0; ///< Parameter \f[ \theta_0 \f] of the Mie-Gruneisen EoS
  double q; ///< Parameter \f[ q \f] of the Mie-Gruneisen EoS
  double rho0; ///< Parameter \f[ \rho_0 \f] of the Mie-Gruneisen EoS

  double ks; ///< Parameter \f[ K_s \f] of the Mie-Gruneisen EoS
  double npu; ///< Parameter \f[ N_s+1 \f] of the Mie-Gruneisen EoS
  double rhos; ///< Parameter \f[ \rho_s \f] of the Mie-Gruneisen EoS

  double ur; ///< Parameter \f[ u_r \f] of the Mie-Gruneisen EoS
  double cvr; ///< Parameter \f[ Cv_r \f] of the Mie-Gruneisen EoS
  double icvr; ///< Parameter \f[ \frac1{Cv_r} \f] of the Mie-Gruneisen EoS
};


/// @brief Initialize internal energy
/// @param [in] rho_ Density
/// @param [in] T_ Temperature
/// @return Energy
inline double MieGruneisenEOS::initSample(double rho_, double T_) {

  return (cvr*(T_-ur*theta(rho_)) + pote(rho_))*mass;

}

inline double MieGruneisenEOS::gamma(const double& rho) {

  return gammaInf + (gamma0-gammaInf)*exp(q*log(rho0/rho));
  
}

inline double MieGruneisenEOS::theta(const double& rho) {

  double lir = log(rho0/rho);
  return theta0 * exp(-gammaInf*lir + (gamma0-gammaInf)/q*(1-exp(q*lir)));
  
}

inline double MieGruneisenEOS::pote(const double& rho) {

  double x = 1-rhos/rho;
  return ks/rhos*(exp(npu*x)/npu - x)/npu; 
  
}

inline double MieGruneisenEOS::potp(const double& rho) {

  return ks/npu*(exp(npu*(1-rhos/rho)) - 1); 
  
}

inline double MieGruneisenEOS::redu(const double& rho, const double& e) {

  return (e-pote(rho))*icvr/theta(rho) + ur;
  
}


#endif // __MIE_GRUNEISEN_EOS_HPP_INCLUDED
