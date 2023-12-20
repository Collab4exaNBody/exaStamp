/// @file
/// @brief Declaration of functions to read from a Stamp Dump file

#ifndef __PARTICLE_INPUT_HPP_INCLUDED
#define __PARTICLE_INPUT_HPP_INCLUDED


#include <string>

#include "io/outputManager.hpp"


class LegacyHeaderIOStruct;
class MPI__Particle;
class MPI__Mesoparticle;


void readLegacyHeader(LegacyHeaderIOStruct* header, const std::string& filename, int init_step, InputOutputManager* ioManager);


void readLegacyParticles(InputOutputManager::FileId& fileId, InputOutputManager::Offset& offset, uint64_t size, MPI__Particle* particles);

void readLegacyDPDEParticles(InputOutputManager::FileId& fileId, InputOutputManager::Offset& offset, uint64_t size, MPI__Mesoparticle* particles);

#endif // __PARTICLE_INPUT_HPP_INCLUDED

