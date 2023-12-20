/// @file 
/// @brief Definition of the EAMPotential subclass

#ifndef __EAM_POTENTIAL_HPP_INCLUDED
#define __EAM_POTENTIAL_HPP_INCLUDED


#include "potential/shortRangePotential.hpp"


/// @brief Interface for an EAM potential
///
/// An EAM potential is a potential where the interaction energy between two atoms
/// depends on the atomic density around these atoms.
/// The energy on an atom i whose neighbors list is \f$ N(i) \f$ can be written in the general form :
/// \f[ \frac{1}{2} \sum_{j \in N(i)}{\phi(r_{ij})} + F(\sum_{j \in N(i)}{\rho_j})\f]
/// where the density contribution \f$ \rho_j \f$ depend on the neighbors of atom j
class EAMPotential : public ShortRangePotential {

public:

  /// @brief Get the subclass of the potential
  /// @return Subclass
  virtual Traversal getTraversal() const {
    return Traversal::EAM;
  }

	/// @brief Get the ghost thickness to use with this potential
  ///
	///
  virtual uint getGhostThickness() const {
    return 1;
  }

  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const = 0;

  /// @brief Constructor
  /// @param [in] rCut Cutoff Radius
  EAMPotential(double rCut)
    : ShortRangePotential(rCut),
      phiCut(0.), rhoCut(0.), fCut(0.) {}

  /// @brief Destructor (nothing to do)
  virtual ~EAMPotential() {}

  /// @brief Calculate the \f$ \phi(r) \f$ term of the energy and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] phi \f$ \phi(r) \f$ term
  /// @param [in] dphi Derivative of \f$ \phi(r) \f$ with respect to the interatomic distance
  virtual void phi(double r, double& phi, double& dphi) = 0;
  /// @brief Calculate the density contribution of a neighbor atom (\f$ \rho(r) \f$) and its derivative
  /// @param [in] r Interatomic distance
  /// @param [in] rho Density contribution
  /// @param [in] drho Derivative of the density contribution with respect to the interatomic distance
  virtual void rho(double r, double& rho, double& drho) = 0;
  /// @brief Calculate the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term of the energy and its derivative
  /// @param [in] rho Sum of the density contributions on the neighbors \f$ \sum_{j \in N(i)}{\rho_j} \f$
  /// @param [in] f F(rho)
  /// @param [in] df Derivative of F(rho) with respect to rho
  virtual void fEmbed(double rho, double& f, double& df) = 0;

  /// @brief Accessor to the \f$ \phi(r) \f$ term at cutoff radius
  inline double getPhiCut() { return phiCut; }
  /// @brief Accessor to the density contribution at cutoff radius
  inline double getRhoCut() { return rhoCut; }
  /// @brief Accessor to the \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term at cutoff radius
  inline double getFCut()   { return fCut; }

public:

  /// @brief Calculate shifts used to render the potential continuous at cutoff distance :
  /// calculate \f$ \phi(r) \f$ term, density contribution and \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term at cutoff radius
  virtual void initializeCutoffVariables() {
    double phiAtRCut=0.0, dPhi=0.0;
    double rhoAtRCut=0.0, dRho=0.0;
    phiCut = 0.0;
    rhoCut = 0.0;
    this->phi(cutoffRadius, rhoAtRCut, dRho);
    this->rho(cutoffRadius, phiAtRCut, dPhi);
    rhoCut = rhoAtRCut;
    phiCut = phiAtRCut;
  }

  double phiCut; ///< \f$ \phi(r) \f$ term at cutoff radius
  double rhoCut; ///< Density contribution at cutoff radius
  double fCut; ///< \f$ F(\sum_{j \in N(i)}{\rho_j}) \f$ term at cutoff radius

};



#endif // __EAM_POTENTIAL_HPP_INCLUDED
