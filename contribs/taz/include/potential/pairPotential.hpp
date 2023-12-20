/// @file 
/// @brief Definition of the PairPotential subclass

#ifndef __PAIR_POTENTIAL_HPP_INCLUDED
#define __PAIR_POTENTIAL_HPP_INCLUDED


#include "potential/shortRangePotential.hpp"


/// @brief Base class for pair potentials
///
/// A pair potential is a potential where the interaction energy between two atoms
/// is a function of the positions of these atoms
class PairPotential : public ShortRangePotential {

public:

  /// @brief Get the subclass of the potential
  /// @return Subclass
  virtual Traversal getTraversal() const {
    return Traversal::PAIR;
  }

	/// @brief Get the ghost thickness to use with this potential
  ///
	/// The ghost thickness is 1 since only the position of the neighbors is used
  /// @return Ghost thickness
  virtual uint getGhostThickness() const {
    return 1;
  }

  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const = 0;

  /// @brief Constructor from a cutoff radius
  /// @param [in] rcut Cutoff radius
  PairPotential(double rcut) : ShortRangePotential(rcut), eCut(0.) {}

  /// @brief Destructor (nothing to do)
  virtual ~PairPotential() {}

  /// @brief Function call operator : calculate the force and energy for the potential
  /// @param [in] r Interatomic distance
  /// @param [out] e Energy
  /// @param [out] de Energy derivative with respect to the interatomic distance
  virtual void operator () (double r, double &e, double &de) = 0;

  /// @brief Accessor to the energy at cutoff radius
  double getEcut() {
    return eCut;
  }


protected:

  /// @brief Calculate shifts used to render the potential continuous at cutoff distance :
  /// calculate energy at cutoff radius
  virtual void initializeCutoffVariables() {
    double tmp;
    (*this)(cutoffRadius, eCut, tmp);
  }

  double eCut; ///< Energy at the cutoff radius

};

#endif // __PAIR_POTENTIAL_HPP_INCLUDED
