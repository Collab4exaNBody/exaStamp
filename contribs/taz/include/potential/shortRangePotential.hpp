/// @file 
/// @brief Definition of ShortRangePotential class

#ifndef __SHORT_RANGE_POTENTIAL_HPP_INCLUDED
#define __SHORT_RANGE_POTENTIAL_HPP_INCLUDED


#include "potential/potential.hpp"


/// @brief Base class for short range potentials
///
/// Short range potentials are potentials which limit the interactions between remote particles
/// by mean of a cutoff radius to speed up the calculation
class ShortRangePotential : public IPotential {

public:

	/// @brief Get the ghost thickness to use with this potential
	/// @return Ghost thickness
  virtual uint getGhostThickness() const = 0;
  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const = 0;

  /// @brief Destructor
  virtual ~ShortRangePotential() {}

  /// @brief Accessor to the cutoff radius
  double getCutoffRadius() const {
    return cutoffRadius;
  }

protected:

  /// @brief Constructor from a cutoff radius
  /// @param [in] rcut Cutoff radius
  ShortRangePotential(double rcut) : cutoffRadius(rcut) {}

  /// @brief Calculate shifts used to render the potential continuous at cutoff distance
  virtual void initializeCutoffVariables() = 0;

  double cutoffRadius;  ///< Cut-off radius

};



#endif // __SHORT_RANGE_POTENTIAL_HPP_INCLUDED
