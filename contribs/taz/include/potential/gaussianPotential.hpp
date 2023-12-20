/// @file 
/// @brief Definition of gaussian potential

#ifndef __GAUSSIAN_POTENTIAL_HPP_INCLUDED
#define __GAUSSIAN_POTENTIAL_HPP_INCLUDED


#include "potential/pairPotential.hpp"

#include "simd/gpot.hpp"


/// @brief Structure to handle parameters of a Gaussian potential
/// 
struct GaussianParameters {

  /// @brief Default constructor
  ///
  ///
  GaussianParameters() : r2Attractive(1.), r2Repulsive(1.), ratio(1.), epsilon(0.) {}

  /// @brief Argument constructor
  /// @param [in] ra Attractive radius
  /// @param [in] rb Repulsive radius
  /// @param [in] A Attractive/repulsive ratio
  /// @param [in] e Epsilon
  GaussianParameters(double ra, double rb, double A, double e) : r2Attractive(ra*ra), r2Repulsive(rb*rb), ratio(A), epsilon(e) {}

  /// @brief Destructor (nothing to do)
  ~GaussianParameters() {}

  double r2Attractive;   ///< Square radius of the attractive gaussian
  double r2Repulsive;    ///< Square radius of the repulsive gaussian
  double ratio;          ///< Attractive/repulsive ratio
  double epsilon;        ///< Force magnitude
};


/// @brief Gaussian potential
/// @see Lei, Mundy, Schenter, Voulgarakis in J. Chem. Phys., 2015
class GaussianPotential : public PairPotential {

public:

  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that potential
  /// @brief Shorcut for the simd version of the potential
  typedef simd::kernels::Gaussian<double> simd_opt_t;

  /// @brief Get the type of the potential
  /// @return Type
  virtual Type getType() const {
    return Type::GAUSSIAN;
  }

  /// @brief Get the name of the potential
  /// @return Name
  virtual std::string getName() const {
    return "Gaussian";
  }

  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const { return 1; }

  /// @brief Constructor with arguments
  /// @param [in] rcut Cutoff radius
  /// @param [in] gParam Parameters
  GaussianPotential(const double& rcut, const GaussianParameters& gParam) 
    : PairPotential(rcut), p(gParam) {
    this->initializeCutoffVariables();
  }

  /// @brief Destructor
  ~GaussianPotential() {}

  /// @brief Function call operator : calculate the force and energy for the potential
  /// 
  /// The energy is computed according to the following formula : 
  /// \f[ V(r) \, = \, \epsilon \, \left[ Ar_a^2 \exp \left( -\frac12 \frac{r^2}{r_a^2}\right) - r_b^2\exp\left( -\frac12 \frac{r^2}{r_b^2} \right) \right] \f]
  /// @param [in] r Interatomic distance
  /// @param [out] e Energy
  /// @param [out] de Energy derivative with respect to the interatomic distance
  virtual inline void operator () (double r, double& e, double& de) {

    double expAtt = exp(-0.5*r*r/p.r2Attractive);
    double expRep = exp(-0.5*r*r/p.r2Repulsive);

    e = p.epsilon * (p.ratio*p.r2Attractive*expAtt-p.r2Repulsive*expRep) - eCut;
    de = p.epsilon * (p.ratio*expAtt-expRep);
  }

  /// @brief Acessor the the parameters
  const GaussianParameters& getParameters() { return p; }

  /// @brief Set the parameters for the simd version of the potential
  /// @param [in,out] opt Simd version of the potential
  void setParameters(simd_opt_t& opt) {
    opt.setParameters(p.r2Attractive, p.r2Repulsive, p.ratio, p.epsilon, eCut);
  }

private:

  GaussianParameters p; ///< Parameters

};

#endif // __GAUSSIAN_POTENTIAL_HPP_INCLUDED
