/// @file 
/// @brief class TypeSmoothWallParticle

#ifndef __TYPE_SMOOTH_WALL_PARTICLE_HPP_INCLUDED
#define __TYPE_SMOOTH_WALL_PARTICLE_HPP_INCLUDED


#include "particle/types/typeSmoothParticle.hpp"


/// @brief Particle type for wall particles in SDPD
class TypeSmoothWallParticle : public TypeSmoothParticle {

public:

  /// @brief Constructor for virtual smooth particles
  /// The inverse mass is set to 0 to prevent any movement
  /// @param [in] name_ Name
  /// @param [in] mass_ Mass
  /// @param [in] unitMass_ Mass of a microscopic particle
  /// @param [in] smoothingLength_ Smoothing length
  /// @param [in] velocity_ Wall velocity
  TypeSmoothWallParticle(const std::string& name_, double mass_, double unitMass_, double smoothingLength_, vec3<double> velocity_)
    : TypeSmoothParticle(name_, mass_, unitMass_,0.,0.,smoothingLength_), velocity(velocity_) {
    oneOverMass = 0.;
  }
  
  /// @brief Destructor
  virtual ~TypeSmoothWallParticle() {}

  /// @brief Accessor to the smoothing length
  /// @return Smoothing length
  virtual inline double getSmoothingLength() const {
    return smoothingLength;
  }

  /// @brief Accessor to the wall velocity
  /// @return Wall velocity
  virtual inline vec3<double> getVelocity() const {
    return velocity;
  }

protected:

  vec3<double> velocity; ///< Constant velocity at which the wall moves
};


#endif // __TYPE_SMOOTH_WALL_PARTICLE_HPP_INCLUDED
