/// @file 
/// @brief Class TypeParticle

#ifndef __TYPE_PARTICLE_HPP_INCLUDED
#define __TYPE_PARTICLE_HPP_INCLUDED


#include <string>

#include "utils/vec3/vec3.hpp"


/// @brief Class to store the type dependent properties of a particle
class TypeParticle {

public:

  /// @brief Constructor
  /// @param [in] name_ Name
  /// @param [in] mass_ Mass
  TypeParticle(const std::string& name_, double mass_)
    : name(name_), mass(mass_), oneOverMass(1./mass_) {}

  /// @brief The destructor
  virtual ~TypeParticle() {}

  /// @brief Accessor to the name
  inline const std::string& getName() const {
    return name;
  }

  ///@brief Accessor to the mass
  inline double getMass() const { 
    return mass; 
  }

  /// @brief Accessor to the mass inverse
  inline double getInvMass() const { 
    return oneOverMass; 
  }

  /// @brief Get size of particle
  /// @return 1 (not a mesoparticle)
  virtual inline double getSize() const {
    return 1;
  }

  /// @brief Get fricition coefficient in the parallel direction
  virtual inline double getGammaParallel() const {
    return 0.;
  }

  /// @brief Get fricition coefficient in the orthogonal directions
  virtual inline double getGammaOrthogonal() const {
    return 0.;
  }

  /// @brief Get smoothing length
  virtual inline double getSmoothingLength() const {
    return 0.;
  }

  /// @brief Get velocity
  virtual inline vec3<double> getVelocity() const {
    return 0.;
  }

protected:

  std::string name; ///< Name

  double mass; ///< Mass
  double oneOverMass; ///< Mass inverse

};

#endif // __TYPE_PARTICLE_HPP_INCLUDED
