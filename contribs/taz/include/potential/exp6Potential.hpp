/// @file 
/// @brief Definition of exponential-6 potential

#ifndef __EXP6_POTENTIAL_HPP_INCLUDED
#define __EXP6_POTENTIAL_HPP_INCLUDED


#include "potential/pairPotential.hpp"

#include "utils/auxMath.hpp"


/// @brief Structure to handle parameters of an exponential-6 potential
struct Exp6Parameters {
public:

  /// @brief Default constructor
	///
	///
  Exp6Parameters() : A(0.), B(1.), C(0.), D(0.) {}

  /// @brief Arguments constructor
  /// @param [in] a Parameter A
  /// @param [in] b Parameter B, must not be zero
  /// @param [in] c Parameter C
  /// @param [in] d Parameter D
  Exp6Parameters(double a, double b, double c, double d) : A(a), B(b), C(c), D(d) {}

  /// @brief Destructor (nothing to do)
  ~Exp6Parameters() {}

  double A; ///< Parameters A
  double B; ///< Parameters B, must not be zero
  double C; ///< Parameters C
  double D; ///< Parameters D
};


/// @brief Exponential-6 potential
/// @warning This potential is not correctly initialized and should not be used until further developments
/// @note This potential is called Exp6v1 in the original Stamp
/// @note Is this different from Buckingham potential ?
class Exp6Potential : public PairPotential {

public:

  static constexpr bool has_simd_opt = false; ///< Indicates that there is no simd version of that potential

  /// @brief Get the type of the potential
  /// @return Type
  Type getType() const { 
    return Type::EXP6;
  }

  /// @brief Get the name of the potential
  /// @return Name
  virtual std::string getName() const {
    return "Exp-6";
  }

  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const {
    return 1.5;
  }
  

  /// @brief Constructor with arguments
  /// @param [in] rcut Cutoff radius
  /// @param [in] exp6Param Parameters
  Exp6Potential(const double& rcut, const Exp6Parameters& exp6Param) 
    : PairPotential(rcut), p(exp6Param) {
    this->initializeCutoffVariables();
  }

  /// @brief Destructor (nothing to do)
  ~Exp6Potential() {}
  
  /// @brief Function call operator : calculate the force and energy for the potential
  ///
  /// The energy is computed according to the following formula : 
  /// \f[ V(r) \, = A \exp\left(-Br\right) \,-\, \frac{C}{r^6} \,+\, D\left(\frac{12}{Br}\right)^{12} \f]
  /// @param [in] r Interatomic distance
  /// @param [out] e Energy
  /// @param [out] de Energy derivative with respect to the interatomic distance
  /// @warning Energy derivative is not calculated here (set to 0)
  virtual void operator () (double r, double& e, double& de) {
    
    double t1, t2, t3;
    
    t1 = p.A * auxExp(-p.B*r);
    
    double r6=r*r*r;
    r6*=r6;
    t2 = p.C/r6;
    
    double ratio = 12/(p.B*r);
    ratio*= ratio*ratio;
    ratio*=ratio;
    ratio*=ratio;
    t3 = p.D*ratio;
    
    e = t1 - t2 + t3;
    de = 0.;
  }
  
private:

  Exp6Parameters p; ///< Parameters

};

#endif // __EXP6_POTENTIAL_HPP_INCLUDED
