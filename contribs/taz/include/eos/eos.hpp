///@file
///@brief Interface for EOS object

#ifndef __EOS_HPP_INCLUDED
#define __EOS_HPP_INCLUDED


#include <string>


/// @brief Interface for a @c EOS
class IEOS {
  
public:
  /// @brief EOS type
  enum Type {
    DPDE_ET, ///< Simple EOS for DPDE
    IDEAL_GAS, ///< Ideal gas EOS
    MIE_GRUNEISEN, ///< Mie-Gruneisen EOS
    HZ, ///< HZ EOS
    JWL, ///< Jones-Wilkins-Lee EOS
    REACTIVE, ///< Reactive EOS
  };


  /// @brief Destructor
  virtual ~IEOS() {}

  /// @brief Get EOS type
  virtual Type getType() const = 0;

  /// @brief Get EOS name
  virtual std::string getName() const = 0;

  /// @brief Initialize energy
  /// @param [in] rho_ Density
  /// @param [in] T_ Temperature
  /// @return Internal energy
  virtual double initSample(double rho_, double T_) = 0;  
};



#endif // __EOS_HPP_INCLUDED
