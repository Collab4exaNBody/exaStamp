/// @file 
/// @brief Definition of tools to print the particle output in DAT and VTK formats

#ifndef __PARTICLE_OUTPUT_HPP_INCLUDED
#define __PARTICLE_OUTPUT_HPP_INCLUDED


#include <string>


#include "utils/array/array.hpp"
#include "utils/vec3/vec3.hpp"


/// @brief Define the maximum number of digit in the number of steps
#define LOG_MAX_STEPS 9


class CommManager;


/// @brief Store particles data used in output writing
class ParticleOutput {
public:
  /// @brief Default constructor
  ParticleOutput():writeType(false),writeVelocity(false),writeEint(false),writeProgress(false) {}
  /// @brief Destructor (nothing to do)
  virtual ~ParticleOutput() {}

  /// @brief Constructor when you don't know the number of particles
  /// @param [in] writeT Indicate if the types are to be written
  /// @param [in] writeV Indicate if the velocities are to be written
  /// @param [in] writeE Indicate if the internal energies are to be written
  /// @param [in] writeP Indicate if the progress variables are to be written
  // When you don't know the size, it's to gather all on a node
  ParticleOutput(bool writeT, bool writeV, bool writeE, bool writeP) 
    : writeType(writeT), writeVelocity(writeV), writeEint(writeE), writeProgress(writeP) {}

  /// @brief Constructor when you do know the number of particles
  /// @param [in] n Number of particles
  /// @param [in] writeT Indicates if the types are to be written
  /// @param [in] writeV Indicates if the velocities are to be written
  /// @param [in] writeE Indicate if the internal energies are to be written
  /// @param [in] writeP Indicate if the progress variables are to be written
  // When you know it, it's the final on 0
  ParticleOutput(uint n, bool writeT, bool writeV, bool writeE, bool writeP) 
    : writeType(writeT), writeVelocity(writeV), writeEint(writeE), writeProgress(writeP),
      id(n), type(writeType ? n : 0), r(n), v(writeVelocity ? n : 0), ei(writeEint ? n : 0), progress(writeProgress ? n : 0) {}

  void gather(CommManager* comm, ParticleOutput* sendBuffer, Array<int>& counts, Array<int>& disps);


  bool writeType; ///< Indicates if the types must be written
  bool writeVelocity; ///< Indicates if the velocities must be written
  bool writeEint; ///< Indicates if the internal energies must be written
  bool writeProgress; ///< Indicates if the progress variables must be written

  Array<uint> id; ///< Indexes of the particles
  Array<uint8_t> type; ///< Types of the particles

  Array< vec3<double> > r; ///< Positions of the particles
  Array< vec3<double> > v; ///< Velocities of the particles
  Array<double> ei; ///< Internal energies of the particles
  Array<double> progress; ///< Progress variables of the particles

};


/// @brief Base tool to write the particles in DAT or VTK format
class ParticleWriter {
public:
  /// @brief Default constructor
  ParticleWriter() : particles(nullptr) {}

  /// @brief Destructor (nothing to do)
  virtual ~ParticleWriter() {}


  /// @brief Constructor
  /// @param [in] buff The particles to be written
  /// @param [in] writeT Indicates if the types are to be written
  /// @param [in] writeV Indicate if the velocities are to be written
  /// @param [in] writeE Indicate if the velocities are to be written
  /// @param [in] writeP Indicate if the progress variables are to be written
  /// @param [in] name Root name of the file where the particles will be written
  /// @param [in] ext Extension for the file where the particles will be written
  ParticleWriter(ParticleOutput* buff, bool writeT, bool writeV, bool writeE, bool writeP, const std::string& name, const std::string& ext)
    : rootName(name), extension(ext), particles(buff) {
  }

  void setFilename(uint step, std::string& filename);

  std::string rootName; ///< Root name of the file where the particles will be written
  std::string extension; ///< Extension for the file where the particles will be written

  ParticleOutput* particles; ///< Storage for the particles

};

#endif // __PARTICLE_OUTPUT_HPP_INCLUDED
