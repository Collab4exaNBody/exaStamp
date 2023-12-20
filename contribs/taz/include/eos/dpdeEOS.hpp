/// @file
/// @brief Definition of DPDEEOS class

#ifndef __DPDEEOS_HPP_INCLUDED
#define __DPDEEOS_HPP_INCLUDED


#include "eos/eos.hpp"

#include "utils/stampUnits.hpp"


/// @brief Class for DPDE equation of state
class DPDEEOS : public IEOS{
public:
  static constexpr bool has_simd_opt = false; ///< Indicates that there is a simd version of that potential
  /// @brief Constructor
  /// @param [in] cv_ Heat capacity
  DPDEEOS(const double& cv_) : cv(cv_), icv(1./cv) {}

  /// @brief Destructor
  virtual ~DPDEEOS() {}

  /// @brief Get EOS type
  /// @return EOS type
  virtual Type getType() const {
    return Type::DPDE_ET;
  }

  /// @brief Get EOS name
  /// @return EOS name
  virtual std::string getName() const {
    return "DPDE (constant)";
  }

  virtual double initSample(double rho_, double T_);

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

private:
  double cv; ///< Heat capacity
  double icv; ///< Inverse heat capacity

};


/// @brief Initialize internal energy
/// @param [in] rho_ Density
/// @param [in] T_ Temperature
/// @return Energy
inline double DPDEEOS::initSample(double rho_, double T_) {

  return getEnergy(T_);

}

#endif // __DPDEEOS_HPP_INCLUDED
