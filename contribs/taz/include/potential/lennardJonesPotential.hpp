/// @file 
/// @brief Definition of Lennard-Jones potential

#ifndef __LENNARD_JONES_POTENTIAL_HPP_INCLUDED
#define __LENNARD_JONES_POTENTIAL_HPP_INCLUDED


#include "potential/pairPotential.hpp"

#include "simd/lj.hpp"


/// @brief Structure to handle parameters of a Lennard-Jones potential
struct LennardJonesParameters {

  /// @brief Default constructor
  ///
  ///
  LennardJonesParameters() : epsilon(0.), sigma(0.) {}

  /// @brief Argument constructor
  /// @param [in] e Epsilon
  /// @param [in] s Sigma
  LennardJonesParameters(double e, double s) : epsilon(e), sigma(s) {}

  /// @brief Destructor (nothing to do)
  ~LennardJonesParameters() {}

  double epsilon;  ///< Depth of the potential well
  double sigma;    ///< Finite distance at which the inter-particle potential is zero
};


/// @brief Lennard-Jones potential
/// @see http://en.wikipedia.org/wiki/Lennard-Jones_potential
class LennardJonesPotential : public PairPotential {

public:

  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that potential
  /// @brief Shorcut for the simd version of the potential
  typedef simd::kernels::LennardJones<double> simd_opt_t;

  /// @brief Get the type of the potential
  /// @return Type
  virtual Type getType() const {
    return Type::LJ;
  }

  /// @brief Get the name of the potential
  /// @return Name
  virtual std::string getName() const {
    return "Lennard-Jones";
  }

  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const { return 1; }

  /// @brief Constructor with arguments
  /// @param [in] rcut Cutoff radius
  /// @param [in] ljParam Parameters
  LennardJonesPotential(const double& rcut, const LennardJonesParameters& ljParam) 
    : PairPotential(rcut), p(ljParam) {
    this->initializeCutoffVariables();
  }

  /// @brief Destructor
  ~LennardJonesPotential() {}

  /// @brief Function call operator : calculate the force and energy for the potential
  /// 
  /// The energy is computed according to the following formula : 
  /// \f[ V(r) \, = \, 4\epsilon \, \left[ \left( \frac{\sigma}{r}\right)^{12} - \left( \frac{\sigma}{r}\right)^6 \right] \f]
  /// @param [in] r Interatomic distance
  /// @param [out] e Energy
  /// @param [out] de Energy derivative with respect to the interatomic distance
  virtual inline void operator () (double r, double& e, double& de) {

    double ir = 1./r;
    double ratio = p.sigma*ir;
    double ratio6 = ratio*ratio*ratio; ratio6*=ratio6;
    double ratio12 = ratio6*ratio6;

    e = 4*p.epsilon * (ratio12-ratio6) - eCut;
    de = 24*p.epsilon * (2.*ratio12-ratio6) * ir;
  }

  /// @brief Acessor the the parameters
  const LennardJonesParameters& getParameters() { return p; }

  /// @brief Set the parameters for the simd version of the potential
  /// @param [in,out] opt Simd version of the potential
  void setParameters(simd_opt_t& opt) {
    opt.setParameters(p.sigma, p.epsilon, eCut);
  }

private:

  LennardJonesParameters p; ///< Parameters

};

#endif // __LENNARD_JONES_POTENTIAL_HPP_INCLUDED
