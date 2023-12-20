/// @file 
/// @brief Interface for a potential object

#ifndef __POTENTIAL_HPP_INCLUDED
#define __POTENTIAL_HPP_INCLUDED


#include <string>


//     A potential is a model which means to describe how to particles
// interact. The simplest ones ( @c PairPotential) contain only one function
// which enables to compute the energy blablabla to be continued ...

/// @brief Base class to define a potential
class IPotential {

public:

  // Potential type :
  //     As you can see in the inheritance diagram, there are a lot of
  // different potentials. The outer nodes of this tree are all implementations
  // of a physical potentials, and their parents define potential types. To
  // type their is a specific set of function to compute the interaction.

	/// @brief Enumeration of the subclasses of potential
	///
	/// In a subclass of potential you will find potentials using the same scheme to compute the forces
	/// but differing by their mathematical formulation
  enum Traversal {
  	IDEAL_GAS, ///< Ideal gas interaction
  	PAIR, ///< Pair potential
  	EAM, ///< EAM potential
  };

  /// @brief Enumeration of the types of potential
  ///
	/// Each type of potential defines a specific mathematical formulation
  enum Type { 
    IDEAL_GAS_, ///< Ideal gas interaction
    EXP6, ///< Exponential 6 potential
    LJ, ///< Lennard-Jones potential
    SUTTON_CHEN, ///< Sutton-Chen potential
    EAM_VNIITF, ///< EAM-Vniitf potential
    MEAM_, ///< MEAM potential
    GAUSSIAN ///< Gaussian potential
  };

  /// @brief Destructor (nothing to do)
  virtual ~IPotential() {}

  /// @brief Get the subclass of the potential
  /// @return Subclass
  virtual Traversal getTraversal() const = 0;
  /// @brief Get the type of the potential
  /// @return Type
  virtual Type getType() const = 0;
  /// @brief Get the name of the potential
  /// @return Name
  virtual std::string getName() const = 0;

  /// @brief Get the cost to apply the potential
  /// @return Cost
  virtual double cost() const = 0;

};


/// @brief Not meant to be used, it is just to protect the name
class Potential {};

#endif // __POTENTIAL_HPP_INCLUDED
