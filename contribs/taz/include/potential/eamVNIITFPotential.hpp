/// @file 
/// @brief Definition of EAM VNIITF potential

#ifndef __EAM_VNIITF_POTENTIAL_HPP_INCLUDED
#define __EAM_VNIITF_POTENTIAL_HPP_INCLUDED


#include "potential/eamPotential.hpp"

#include "simd/eamVNIITF.hpp"

#include "utils/spline.hpp"


/// @brief Structure to handle parameters of a EAM VNIITF potential
struct EamVniitfParameters {

	/// @brief Default constructor
  EamVniitfParameters() {}

  /// @brief Constructor
	/// @param [in] rmax_ Parameter rmax
  /// @param [in] rmin_ Parameter rmin
  /// @param [in] rt0_ Parameter rt0
  /// @param [in] Ecoh_ Parameter Ecoh
  /// @param [in] E0_ Parameter E0
  /// @param [in] beta_ Parameter beta
  /// @param [in] A_ Parameter A
  /// @param [in] Z_ Parameter Z
  /// @param [in] n_ Parameter n
  /// @param [in] alpha_ Parameter alpha
  /// @param [in] D_ Parameter D
  /// @param [in] mu_ Parameter mu
  /// @param [in] eta_ Parameter eta
  EamVniitfParameters(double rmax_, double rmin_, double rt0_, double Ecoh_, double E0_, double beta_, 
		      double A_, double Z_, double n_, double alpha_, double D_, double eta_, double mu_)
    : rmax(rmax_), rmin(rmin_), rt0(rt0_), Ecoh(Ecoh_), E0(E0_), beta(beta_), A(A_), Z(Z_), n(n_), alpha(alpha_), D(D_), eta(eta_), mu(mu_) {}

  /// @brief Destructor (nothing to do)
  ~EamVniitfParameters() {}

  double rmax; ///< Maximal distance for density contribution
  double rmin; ///< Minimal distance for density contribution
  double rt0; ///< Characteristic radius
  double Ecoh; ///< Cohesive energy
  double E0; ///< Energy parameter
  double beta; ///< Attenuation rate
  double A; ///< Empirical parameter A
  double Z; ///< Number of neighbors in the reference structure
  double n; ///< Empirical parameter n
  double alpha; ///< Parameter alpha
  double D; ///< Parameter D
  double eta; ///< Parameter eta
  double mu; ///< Parameter mu

};


/// @brief EAM VNIITF potential 
///
/// This is a EAM potential developed to calculate the properties of tin near melting curve
/// @see http://iopscience.iop.org/article/10.1088/1742-6596/500/3/032017/pdf
class EamVniitfPotential : public EAMPotential {

public:

  static constexpr bool has_simd_opt = true; ///< Indicates that there is a simd version of that potential
  /// @brief Shorcut for the simd version of the potential
  typedef simd::kernels::EamVniitf<double> simd_opt_t;

  /// @brief Get the type of the potential
  /// @return Type
  virtual Type getType() const {
    return Type::EAM_VNIITF;
  }

  /// @brief Get the name of the potential
  /// @return Name
  virtual std::string getName() const {
    return "EAM-Vniitf";
  }

  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const {
    return 2.;
  }

  /// @brief Constructor
	/// @param [in] rcut Cutoff radius
  /// @param [in] parameters Parameters
  EamVniitfPotential(const double& rcut, const EamVniitfParameters& parameters) 
    : EAMPotential(rcut), p(parameters) {
    this->initializeCutoffVariables();
  }

  /// @brief Destructor (nothing to do)
  virtual ~EamVniitfPotential() {}

  /// @brief Calculate the \f$ \phi(r) \f$ term of the energy and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] phi \f$ \phi(r) \f$ term
  /// @param [in] dphi Derivative of \f$ \phi(r) \f$ with respect to the interatomic distance
  inline virtual void phi(double r, double& phi, double& dphi) override final {

    double ir   = 1/r;
    double irt0 = 1/p.rt0;
    double dr   = r*irt0 - 1.0;
    double dr2  = dr*dr;
    double a    = -2*p.Ecoh/p.Z;
    double b    = p.alpha*p.alpha*p.alpha*p.D*p.rt0;

    double f1  = a*( 1 + p.alpha*dr + p.eta*dr2 + (p.mu+b*ir)*dr2*dr );
    double df1 = a*( p.alpha*irt0 + 2*p.eta*irt0*dr + 3*p.mu*irt0*dr2 + b*(3*irt0/p.alpha-dr*ir)*dr2*ir );

    double f2  = auxExp(-p.alpha*dr);
    double df2 = -p.alpha*irt0*f2;

    double drS = (p.rmax-r)/(p.rmax-p.rmin);
    double S  = Spline::S<3>(drS);
    double dS = Spline::dS<3>(drS) / (p.rmin-p.rmax);

    // Maybe this could be made before
    if (drS<0) {
      S  = 0.;
      dS = 0.;
    }
    else if (drS>1) {
      S  = 1.;
      dS = 0.;
    }

    phi  = (p.E0 + f1*f2) *  S;
    dphi = (p.E0 + f1*f2) * dS + (f1*df2 + f2*df1) * S;

  }

  /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$) and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] rho Density contribution
  /// @param [in] drho Derivative of the density contribution with respect to the interatomic distance
  inline virtual void rho(double r, double& rho, double& drho) override final {

    double irt0 = 1/p.rt0;

    double F  = auxExp(-p.beta * (r*irt0 - 1.0) ) / p.Z;
    double dF = -p.beta*F*irt0;

    double drS = (p.rmax-r)/(p.rmax-p.rmin);
    double S  = Spline::S<3>(drS);
    double dS = Spline::dS<3>(drS) / (p.rmin-p.rmax);

    // Maybe this could be made before
    if (drS<0) {
      S  = 0.;
      dS = 0.;
    }
    else if (drS>1) {
      S  = 1.;
      dS = 0.;
    }

    rho  = F*S;
    drho = F*dS + S*dF;

  }

  /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and its derivative
  /// @param [in] rho Sum of the density contributions on the neighbors \f$ \sum_{j \in N(i)}{\rho_j} \f$
  /// @param [in] f F(rho)
  /// @param [in] df Derivative of F(rho) with respect to rho
  inline virtual void fEmbed(double rho, double& f, double& df)  override final{

    if (rho<=0.) {

      f  = 0.;
      df = 0.;

    }
    else {

      double a = auxPow(rho, p.n);
      double b = p.A*p.Ecoh * a;
      double c = auxLog(a);
      
      f  = b * (c-1);
      df = p.n*b*c / rho; 

    }

  }

  /// @brief Acessor to the parameters
  const EamVniitfParameters& getParameters() { return p; }

  /// @brief Set the parameters for the simd version of the potential
  /// @param [in,out] opt Simd version of the potential
  void setParameters(simd_opt_t& opt) {
    opt.setParameters(p.rmax, p.rmin, p.rt0, p.Ecoh, p.E0, p.beta, p.Z, p.alpha, p.D, p.eta, p.mu, rhoCut, phiCut);
  }

private:

  EamVniitfParameters p; ///< Parameters

};



#endif // __EAM_VNIITF_POTENTIAL_HPP_INCLUDED
