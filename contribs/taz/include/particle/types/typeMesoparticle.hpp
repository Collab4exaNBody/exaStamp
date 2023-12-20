/// @file 
/// @brief class TypeMesoparticle

#ifndef __TYPE_MESOPARTICLE_HPP_INCLUDED
#define __TYPE_MESOPARTICLE_HPP_INCLUDED


#include "particle/types/typeParticle.hpp"


/// @brief Particle type for mesoparticle
class TypeMesoparticle : public TypeParticle {

public:

  /// @brief Constructor
  /// @param [in] name_ Name
  /// @param [in] mass_ Mass
  /// @param [in] unitMass_ Mass of a microscopic particle
  /// @param [in] gammaParallel_ Friction parameter in the parallel direction
  /// @param [in] gammaOrthogonal_ Friction parameter in the orthogonal direction
  TypeMesoparticle(const std::string& name_, double mass_, int unitMass_, double gammaParallel_, double gammaOrthogonal_)
    : TypeParticle(name_, mass_), unitMass(unitMass_), gammaParallel(gammaParallel_), gammaOrthogonal(gammaOrthogonal_) {}

  /// @brief Destructor
  virtual ~TypeMesoparticle() {}

  /// @brief Accessor to the unit mass
  /// @return Mass of one component
  virtual inline double getUnitMass() const {
    return unitMass;
  }

  /// @brief Get size of particle
  /// @return 1 (not a mesoaprticle)
  virtual inline double getSize() const {
    return this->mass/unitMass;
  }

  /// @brief Accessor to the parallel friction parameter
  /// @return Friction parameter 
  virtual inline double getGammaParallel() const {
    return gammaParallel;
  }

  /// @brief Accessor to the orthogonal friction parameter
  /// @return Friction parameter 
  virtual inline double getGammaOrthogonal() const {
    return gammaOrthogonal;
  }

  /// @brief Accessor to the smoothing length (dummy function for mesoparticle)
  /// Should be 0 for non smooth particle
  /// @ return Smoothing length
  virtual inline double getSmoothingLength() const {
    return 0;
  }

  /// @brief Accessor to the particle type velocity
  /// Should be 0 for non wall types
  /// @return Velocity
  virtual inline vec3<double> getVelocity() const {
    return 0.;
  }

protected:

  uint unitMass;  ///< Mass of one consistuent small particle
  double gammaParallel; ///< Friction coefficient in the parallel direction
  double gammaOrthogonal; ///< Friction coefficient in the orthogonal directions

};





#endif // __TYPE_MESOPARTICLE_HPP_INCLUDED
