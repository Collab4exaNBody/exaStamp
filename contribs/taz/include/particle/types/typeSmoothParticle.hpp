/// @file 
/// @brief class TypeSmoothParticle

#ifndef __TYPE_SMOOTH_PARTICLE_HPP_INCLUDED
#define __TYPE_SMOOTH_PARTICLE_HPP_INCLUDED


#include "particle/types/typeMesoparticle.hpp"


/// @brief Particle type for smooth particles
class TypeSmoothParticle : public TypeMesoparticle {

public:

  /// @brief Constructor for real smooth particles
  /// @param [in] name_ Name
  /// @param [in] mass_ Mass
  /// @param [in] unitMass_ Mass of a microscopic particle
  /// @param [in] bulkViscosity_ Bulk viscosity
  /// @param [in] shearViscosity_ Shear viscosity
  /// @param [in] smoothingLength_ Smoothing length
  TypeSmoothParticle(const std::string& name_, double mass_, double unitMass_, double bulkViscosity_, double shearViscosity_, double smoothingLength_)
    : TypeMesoparticle(name_, mass_, unitMass_,(10*shearViscosity_/3.+4*bulkViscosity_)*smoothingLength_/mass_,(5*shearViscosity_/3.-bulkViscosity_)*smoothingLength_/mass_), smoothingLength(smoothingLength_) {}

 
  /// @brief Destructor
  virtual ~TypeSmoothParticle() {}

  /// @brief Accessor to the smoothing length
  /// @return Smoothing length
  virtual inline double getSmoothingLength() const {
    return smoothingLength;
  }

  /// @brief Accessor to the particle type velocity
  /// Should be 0 for non wall types
  /// @return Velocity
  virtual inline vec3<double> getVelocity() const {
    return 0.;
  }

protected:

  double smoothingLength; ///< Smoothing length (for SDPD)
};





#endif // __TYPE_SMOOTH_PARTICLE_HPP_INCLUDED
