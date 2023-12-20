/// @file
/// @brief Gathering of all implementations based on TypeParticle


#include "io/input.hpp"

#include "particle/types/allTypes.hpp"


 /// @brief Constructor from an input structure
 /// @param [in] input Input structure
Configuration<TypeParticle>::Configuration(const Input& input)
  : atomNames(input.atomNames.size()),
    atomAtomicNumbers(input.atomAtomicNumbers.size()),
    atomMasses(input.atomMasses.size()),
    atomCharges(input.atomCharges.size()),
    mesoNames(input.mesoNames.size()),
    mesoMasses(input.mesoMasses.size()),
    mesoUnitMasses(input.mesoUnitMasses.size()),
    mesoGammaPara(input.mesoGammaPara.size()),
    mesoGammaOrtho(input.mesoGammaOrtho.size()),
    smoothNames(input.smoothNames.size()),
    smoothMasses(input.smoothMasses.size()),
    smoothUnitMasses(input.smoothUnitMasses.size()),
    smoothBulkViscosity(input.smoothBulkViscosity.size()),
    smoothShearViscosity(input.smoothShearViscosity.size()),
    smoothSmoothingLength(input.smoothSmoothingLength.size()),
    smoothKernel(input.smoothKernel.size()),
    wallNames(input.wallNames.size()),
    wallMasses(input.wallMasses.size()),
    wallUnitMasses(input.wallUnitMasses.size()),
    wallSmoothingLength(input.wallSmoothingLength.size()),
    wallKernel(input.wallKernel.size()),
    wallVelocity(input.wallVelocityX.size()) {

  
  int numAtoms = input.atomNames.size();
  for (int i=0; i<numAtoms; ++i) {
    atomNames[i] = input.atomNames[i];
    atomAtomicNumbers[i] = input.atomAtomicNumbers[i];
    atomMasses[i] = input.atomMasses[i];
    atomCharges[i] = input.atomCharges[i];
  }

  int numMesos = input.mesoNames.size();
  for (int i=0; i<numMesos; ++i) {
    mesoNames[i] = input.mesoNames[i];
    mesoMasses[i] = input.mesoMasses[i];
    mesoUnitMasses[i] = input.mesoUnitMasses[i];
    mesoGammaPara[i] = input.mesoGammaPara[i];
    mesoGammaOrtho[i] = input.mesoGammaOrtho[i];
  }

  int numSmooths = input.smoothNames.size();
  for (int i=0; i<numSmooths; ++i) {
    smoothNames[i] = input.smoothNames[i];
    smoothMasses[i] = input.smoothMasses[i];
    smoothUnitMasses[i] = input.smoothUnitMasses[i];
    smoothBulkViscosity[i] = input.smoothBulkViscosity[i];
    smoothShearViscosity[i] = input.smoothShearViscosity[i];
    smoothSmoothingLength[i] = input.smoothSmoothingLength[i];
    smoothKernel[i] = input.smoothKernel[i];
  }

  int numWalls = input.wallNames.size();
  for (int i=0; i<numWalls; ++i) {
    wallNames[i] = input.wallNames[i];
    wallMasses[i] = input.wallMasses[i];
    wallUnitMasses[i] = input.wallUnitMasses[i];
    wallSmoothingLength[i] = input.wallSmoothingLength[i];
    wallKernel[i] = input.wallKernel[i];
    wallVelocity[i] = vec3<double>(input.wallVelocityX[i],input.wallVelocityY[i],input.wallVelocityZ[i]);
  }

}
