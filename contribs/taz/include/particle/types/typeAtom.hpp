/// @file 
/// @brief Class TypeAtom

#ifndef __TYPE_ATOM_HPP_INCLUDED
#define __TYPE_ATOM_HPP_INCLUDED


#include "particle/types/typeParticle.hpp"


/// @brief Class to store the type dependent properties of an atom
class TypeAtom : public TypeParticle {

public:

  /// @brief Constructor
  /// @param [in] name_ Name
  /// @param [in] mass_ Mass
  /// @param [in] atomicNumber_ Atomic number
  TypeAtom(const std::string& name_, double mass_, int atomicNumber_) :
  	TypeParticle(name_, mass_),
  	atomicNumber(atomicNumber_)
  	{}

  /// @brief Destructor (nothing to do)
  virtual ~TypeAtom() {}

  /// @brief Accessor to the atomic number
  inline uint getAtomicNumber() const {
    return atomicNumber;
  }
  
protected:

  uint atomicNumber;  ///< Atomic number

// For now charge is not used there and therefore commented.
// If you considered using this charge please first check that
// none of the following solutions is adapted to your needs :
//  - store the charges as parameter of the potential (see the
// example of the Buckingham-Wolf potential);
//  - include the charge in the particle, like in the
// MPI__AtomCharged class (use only if there is a big amount of
// different charges).
//
//  double charge;      ///< Charge

};

#endif // __TYPE_ATOM_HPP_INCLUDED
