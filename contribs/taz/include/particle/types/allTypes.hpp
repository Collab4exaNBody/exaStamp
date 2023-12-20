/// @file 
/// @brief Gathering of all implementations based on TypeParticle

#ifndef __ALL_TYPES_HPP_INCLUDED
#define __ALL_TYPES_HPP_INCLUDED


#include <string>


class Input;
template <class T> class Configuration;
class TypeParticle;


/// @brief Temporary structure gathering the data to initialize all types of particles
template <> class Configuration<TypeParticle> {
public:

  /// @brief Default constructor
  Configuration() {}

  /// @brief Destructor (nothing to do)
  ~Configuration() {}

  Configuration(const Input& input);

  Array<std::string> atomNames;         ///< Atom names
  Array<int>         atomAtomicNumbers; ///< Atom atomic numbers
  Array<double>      atomMasses;        ///< Atom masses
  Array<double>      atomCharges;       ///< Atom charges (atom gather atom and ion)

  Array<std::string> mesoNames;         ///< mesoparticle names
  Array<double>      mesoMasses;        ///< mesoparticle masses
  Array<double>      mesoUnitMasses;    ///< mesoparticle unit masses
  Array<double>      mesoGammaPara;     ///< mesoparticle parallel friction coefficient
  Array<double>      mesoGammaOrtho;    ///< mesoparticle orthogonal friction coefficient

  Array<std::string> smoothNames;             ///< smoothparticle names
  Array<double>      smoothMasses;            ///< smoothparticle masses
  Array<double>      smoothUnitMasses;        ///< smoothparticle unit masses
  Array<double>      smoothBulkViscosity;     ///< smoothparticle bulk viscosity
  Array<double>      smoothShearViscosity;    ///< smoothparticle shear viscosity
  Array<double>      smoothSmoothingLength;   ///< smoothparticle smoothing lengh
  Array<std::string> smoothKernel;            ///< smoothparticle kernel function

  Array<std::string>    wallNames;             ///< wall particle names
  Array<double>         wallMasses;            ///< wall particle masses
  Array<double>         wallUnitMasses;        ///< wall particle unit masses
  Array<double>         wallSmoothingLength;   ///< wall particle smoothing lengh
  Array<std::string>    wallKernel;            ///< wall particle kernel function
  Array< vec3<double> > wallVelocity;          ///< wall velocity

};

#endif // __ALL_TYPES_HPP_INCLUDED
