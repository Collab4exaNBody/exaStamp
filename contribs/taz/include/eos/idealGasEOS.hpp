/// @file
/// @brief Definition of IdealGasEOS class

#ifndef __IDEAL_GAS_EOS_HPP_INCLUDED
#define __IDEAL_GAS_EOS_HPP_INCLUDED


#include "eos/eos.hpp"

#include "simd/igEOS.hpp"

#include "utils/stampUnits.hpp"


/// @brief Class for Ideal Gas equation of state
class IdealGasEOS : public IEOS{
public:
  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that potential
  /// @brief Shorcut for the simd version of the potential
  typedef simd::kernels::IdealGas<double> simd_opt_t;

  /// @brief Constructor
  /// @param [in] size_ Mesoparticle size
  /// @param [in] mass_ Mesoparticle mass
  IdealGasEOS(const uint& size_, const double& mass_) : size(size_), mass(mass_), cv(1.5*(size_-1)*Stamp_Constant::boltzmann), icv(1./cv), cd(cv*2./3.), cp(2./3./mass_) {}

  /// @brief Destructor
  virtual ~IdealGasEOS() {}

  /// @brief Get EOS type
  /// @return EOS type
  virtual Type getType() const {
    return Type::IDEAL_GAS;
  }


  /// @brief Get EOS name
  /// @return EOS name
  virtual std::string getName() const {
    return "Ideal gas";
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

  /// @brief Get energy form density and entropy
  /// @param [in] rho Density
  /// @param [in] S entropy
  /// @return Energy
  inline double getEnergy(const double& rho, const double& S) {
    return exp(S/cv)*cbrt(rho*rho);
  }

  /// @brief Get energy from temperature
  /// @param [in] T Temperature
  /// @return Energy
  inline double getEnergy(const double& T) {
    return cv*T;
  }

  /// @brief Get temperature from energy
  /// @param [in] e Energy
  /// @return Temperature
  inline double getTemperature(const double& e) {
    return e*icv;
  }

  /// @brief Get Inverse temperature from energy
  /// @param [in] e Energy
  /// @return Inverse temperature \f[ \beta \f]
  inline double getBeta(const double& e) {
    return cv/e/Stamp_Constant::boltzmann;
  }

  /// @brief Get entropy, temperature and pressure from density and energy
  /// @param [in] rho Density
  /// @param [in] e Energy
  /// @return Entropy, temperature and pressure
  inline vec3<double> getSTP(const double& rho, const double& e) {
    double S = cv*log(e) - cd*log(rho);
    double T = e*icv;
    double P = rho*e*cp;
    return vec3<double>(S,T,P);
  }

  /// @brief Get energy, temperature and pressure from density and entropy
  /// @param [in] rho Density
  /// @param [in] S Entropy
  /// @return Energy, temperature and pressure
  inline vec3<double> getETP(const double& rho, const double& S) {
    double E = exp(S*icv)*cbrt(rho*rho);
    double T = E*icv;
    double P = rho*E*cp;
    return vec3<double>(E,T,P);
  }

  /// @brief Set the parameters for the simd version of the equation of state
  /// @param [in,out] opt Simd version of the equation of state
  void setParameters(simd::kernels::IdealGas<double>& opt) {
    opt.setParameters(size,mass);
  }

private:
  uint size; ///< Mesoparticle size
  double mass; ///< Mesoparticle mass
  double cv; ///< Heat capacity
  double icv; ///< Inverse heat capacity
  double cd; ///< Coefficient of the density contribution to entropy
  double cp; ///< Constant for pressure

};


/// @brief Initialize internal energy
/// @param [in] rho_ Density
/// @param [in] T_ Temperature
/// @return Energy
inline double IdealGasEOS::initSample(double rho_, double T_) {

  return getEnergy(T_);

}

#endif // __IDEAL_GAS_EOS_HPP_INCLUDED

