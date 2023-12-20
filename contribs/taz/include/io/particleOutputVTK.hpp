/// @file 
/// @brief Definition of a tool to print the particle output in VTK format

#ifndef __PARTICLE_OUTPUT_VTK_HPP_INCLUDED
#define __PARTICLE_OUTPUT_VTK_HPP_INCLUDED


#include <fstream>
#include <string>

#include "io/particleOutput.hpp"


/// @brief Tool to write the particles in VTK format
struct ParticleWriterVTK : public ParticleWriter {

  /// @brief Default constructor
  ParticleWriterVTK() {}
  /// @brief Destructor (nothing to do)
  virtual ~ParticleWriterVTK() {}



  /// @brief Constructor
  /// @param [in] buff The particles to be written
  /// @param [in] writeT Indicates if the types are to be written (not used)
  /// @param [in] writeV Indicate if the velocities are to be written (not used)
  /// @param [in] writeE Indicate if the internal energies are to be written (not used)
  /// @param [in] writeP Indicate if the progress variables are to be written (not used)
  /// @param [in] name Root name of the file where the particles will be written
  ParticleWriterVTK(ParticleOutput* buff,  bool writeT, bool writeV, bool writeE, bool writeP, const std::string& name) 
    : ParticleWriter(buff, writeT, writeV, writeE, writeP, name, ".vtk"), writeCorners(false) {
  }

  void write(uint step);

  void writeHeader(std::ofstream& flux);
  void writePositions(std::ofstream& flux);
  void writeIndexes(std::ofstream& flux);
  void writeTypes(std::ofstream& flux);
  void writeNbZ(std::ofstream& flux);
  void writeVelocities(std::ofstream& flux);
  void writeEints(std::ofstream& flux);
  void writeProgresses(std::ofstream& flux);

  /// @brief Get the number of object to print do to corners writing
  /// @return Number of corners or 0 if they are not written
  uint extraSize() {
    return writeCorners ? 8 : 0 ;
  }

  bool writeCorners; ///< Indicates if the corners of the simulation box must be written (not used)

};

#endif // __PARTICLE_OUTPUT_HPP_INCLUDED
