/// @file 
/// @brief Definition of a tool to print the particle output in DAT format

#ifndef __PARTICLE_OUTPUT_DAT_HPP_INCLUDED
#define __PARTICLE_OUTPUT_DAT_HPP_INCLUDED


#include "io/particleOutput.hpp"


/// @brief Tool to write the particles in DAT format
struct ParticleWriterDAT : public ParticleWriter {

  /// @brief Default constructor
  ParticleWriterDAT() {}
  /// @brief Destructor (nothing to do)
  virtual ~ParticleWriterDAT() {}


  /// @brief Constructor
  /// @param [in] buff The particles to be written
  /// @param [in] writeT Indicates if the types are to be written
  /// @param [in] writeV Indicate if the velocities are to be written
  /// @param [in] writeE Indicate if the velocities are to be written
  /// @param [in] writeP Indicate if the progress variables are to be written
  /// @param [in] name Root name of the file where the particles will be written
  ParticleWriterDAT(ParticleOutput* buff, bool writeT, bool writeV, bool writeE, bool writeP, const std::string& name)
    : ParticleWriter(buff, writeT, writeV, writeE, writeP, name, ".dat") {
  }

  void write(uint step);

};

#endif // __PARTICLE_OUTPUT_HPP_INCLUDED
