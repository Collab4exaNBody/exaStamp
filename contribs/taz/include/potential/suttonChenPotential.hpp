/// @file 
/// @brief Definition of Sutton-Chen potential

#ifndef __SUTTON_CHEN_POTENTIAL_HPP_INCLUDED
#define __SUTTON_CHEN_POTENTIAL_HPP_INCLUDED


#include "potential/eamPotential.hpp"

#include "simd/suttonChen.hpp"


/// @brief Structure to handle parameters of a Sutton-Chen potential
struct SuttonChenParameters {

	/// @brief Default constructor
	///
	///
  SuttonChenParameters()
    : c(0.), epsilon(0.), a0(0.), n(0), m(0) {}

  /// @brief Argument contructor
  /// @param [in] c_ Parameter c
  /// @param [in] e_ Epsilon
  /// @param [in] a_ a0
  /// @param [in] n_ Parameter n
  /// @param [in] m_ Parameter m
  SuttonChenParameters(double c_, double e_, double a_, double n_, double m_)
    : c(c_), epsilon(e_), a0(a_), n(n_), m(m_) {}

  /// @brief Destructor (nothing to do)
  ~SuttonChenParameters() {}

  double c; ///< Shape parameter c
  double epsilon; ///< Energy scale
  double a0; ///< Length scale
  double n; ///< Shape parameter n
  double m; ///< Shape parameter m

};


/// @brief Sutton-Chen potential
///
/// A Sutton-Chen potential is an EAM potential where :
/// \f[ \rho_i = \sum_{j \in N(i)}{\left(\frac{a_0}{r_{ij}}\right)^m} \f]
/// \f[ \phi(r) = \varepsilon \left(\frac{a_0}{r}\right)^n \f]
/// \f[ F(\sum_{j \in N(i)}{\rho_j}) = -c\varepsilon \sqrt{\sum_{j \in N(i)}{\rho_j}} \f]
class SuttonChenPotential : public EAMPotential {

public:

  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that potential
  /// @brief Shorcut for the simd version of the potential
  typedef simd::kernels::SuttonChen<double> simd_opt_t;

  /// @brief Get the type of the potential
  /// @return Type
  virtual Type getType() const {
    return Type::SUTTON_CHEN;
  }

  /// @brief Get the name of the potential
  /// @return Name
  virtual std::string getName() const {
    return "Sutton-Chen";
  }

  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const {
    return 3.;
  }

  /// @brief Constructor
	/// @param [in] rcut Cutoff radius
  /// @param [in] scParam Parameters
  SuttonChenPotential(const double& rcut, const SuttonChenParameters& scParam) 
    : EAMPotential(rcut), p(scParam) {
    this->initializeCutoffVariables();
  }

  /// @brief Destructor (nothing to do)
  virtual ~SuttonChenPotential() {}

  /// @brief Calculate the \f$ \phi(r) \f$ term of the energy and its derivative
  ///
	/// Where \f$ \phi(r) = \varepsilon \left(\frac{a_0}{r}\right)^n \f$
  /// @param [in] r Interatomic distance
  /// @param [in] phi \f$ \phi(r) \f$ term
  /// @param [in] dphi Derivative of \f$ \phi(r) \f$ with respect to the interatomic distance
  inline void phi(double r, double& phi, double& dphi) {
    double ratio = p.a0/r;
    phi  = p.epsilon * auxPow(ratio, p.n);
    dphi = -1 * p.n * phi / r;
    phi -= phiCut;
  }

  /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$) and its derivative
  ///
	/// Where \f$ \rho_i = \sum_{j \in N(i)}{\left(\frac{a_0}{r_{ij}}\right)^m} \f$
  /// @param [in] r Interatomic distance
  /// @param [in] rho Density contribution
  /// @param [in] drho Derivative of the density contribution with respect to the interatomic distance
  inline void rho(double r, double& rho, double& drho) {
    double ratio = p.a0/r;
    rho  = auxPow(ratio, p.m);
    drho = -1 * p.m * rho / r;
    rho -= rhoCut;
  }

  /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and its derivative
  ///
	/// Where \f$ F(\sum_{j \in N(i)}{\rho_j}) = -c\varepsilon \sqrt{\sum_{j \in N(i)}{\rho_j}} \f$
  /// @param [in] rho Sum of the density contributions on the neighbors \f$ \sum_{j \in N(i)}{\rho_j} \f$
  /// @param [in] f F(rho)
  /// @param [in] df Derivative of F(rho) with respect to rho
  inline void fEmbed(double rho, double& f, double& df) {
    double sqrtRho = auxSqrt(rho);
    f  = -1. * p.c * p.epsilon * sqrtRho;
    df = 0.5 * f / (rho>0? rho : 0);
    f -= fCut;
  }

  /// @brief Acessor to the parameters
  const SuttonChenParameters& getParameters() { return p; }

  /// @brief Set the parameters for the simd version of the potential
  /// @param [in,out] opt Simd version of the potential
  void setParameters(simd_opt_t& opt) {
    opt.setParameters(p.epsilon, p.a0, p.m, p.n, rhoCut, fCut, phiCut);
  }

public:

  SuttonChenParameters p; ///< Parameters

};



#endif // __SUTTON_CHEN_POTENTIAL_HPP_INCLUDED
