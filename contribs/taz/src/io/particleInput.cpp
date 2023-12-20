/// @file 
/// @brief Implementation of functions to read from a Stamp Dump file


#include <fstream>
#include <string>

#include "io/outputManager.hpp"
#include "io/StampV3LegacyIOStructures.hpp"

#include "parallel/types/MPI_particle.hpp"
#include "parallel/types/MPI_mesoparticle.hpp"

#include "utils/stampUnits.hpp"


/// @brief Read the header of the legacy file
/// @param [out] header Header to fill
/// @param [in] _filename Name of the legacy file
/// @param [in] init_step Initial step
/// @param [in] ioManager Input/Output manager
void readLegacyHeader(LegacyHeaderIOStruct* header, const std::string& _filename, int init_step, InputOutputManager* ioManager) {

  // Get the file name or create it from the initial step number
  std::string filename=_filename;
  
  if (filename == ""){
    std::string tmp = std::to_string(init_step);
    while (tmp.size() < 9) tmp.insert(0,"0");
    filename="StampV3prot_"+tmp+".MpiIO";
  }
  // Test file
  std::ifstream infile(filename); 
  if (!infile.good()) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": "
	     << "in function 'readLegacyHeader(LegacyHeaderIOStruct*, const std::string&, CommManager&)' : unable to open file " << filename << ". STOP."
	     << std::endl;
    exit(-1);
  }

  // Read the begining of the file file as an array of headers
  InputOutputManager::FileId fileId;
  InputOutputManager::Offset offset = 0;
  Array<LegacyHeaderIOStruct> array(1);

  ioManager->openIOFile(filename, "r", fileId);
  InputOutputManager::read(array, fileId, offset);
  ioManager->closeIOFile(fileId);

  // Keep the first and only element as the header
  *header = array[0];

}


/// @brief Read the header of the Hercule file (not used)
/// @param [out] header Header to fill
/// @param [in] filename Name of the legacy file
/// @param [in] commManager Communication manager (not used)
/// @param [in] inputOutputManager Input/Output manager
void readHerculeHeader(LegacyHeaderIOStruct* header, const std::string& filename, CommManager& commManager,InputOutputManager* inputOutputManager) {

  // Test file
  std::ifstream infile(filename);
  if (!infile.good()) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": "
	     << "in function 'readLegacyHeader(readHerculeHeader*, const std::string&, CommManager&)' : unable to open file " << filename << ". STOP."
	     << std::endl;
    exit(-1);
  }

  // Read the begining of the file file as an array of headers
  InputOutputManager::FileId fileId;
  InputOutputManager::Offset offset = 0;
  Array<LegacyHeaderIOStruct> array(1);

  inputOutputManager->openIOFile(filename, "r", fileId);
  InputOutputManager::read(array, fileId, offset);
  inputOutputManager->closeIOFile(fileId);

  // Keep the first and only element as the header
  *header = array[0];

}


/// @brief Read particles from the legacy file
/// @param [in] fileId Legacy file
/// @param [in] offset Position to read at
/// @param [in] size Number of particles to read
/// @param [out] particles Particles array to fill
void readLegacyParticles(InputOutputManager::FileId& fileId, InputOutputManager::Offset& offset, uint64_t size, MPI__Particle* particles) {

  const auto meterPerSecond = SI_Units_base::meter/SI_Units_base::second;

  // Get an array of legacy-structured particles
  Array<LegacyParticleIOStruct> array(size);
  InputOutputManager::read(array, fileId, offset);

  // Convert the legacy structured particles into MPI_Particles
  for (uint64_t i=0; i<size; ++i) {

    MPI__Particle& tmp = particles[i];
    LegacyParticleIOStruct& tmpLegacy = array[i];

    tmp.id  = tmpLegacy.particleID - 1;
    tmp.ti  = tmpLegacy.particleType;
    tmp.r.x = convert(tmpLegacy.coordinates[0], SI_Units_base::meter, Stamp_Units::length);
    tmp.r.y = convert(tmpLegacy.coordinates[1], SI_Units_base::meter, Stamp_Units::length);
    tmp.r.z = convert(tmpLegacy.coordinates[2], SI_Units_base::meter, Stamp_Units::length);
    tmp.v.x = convert(tmpLegacy.velocity[0], meterPerSecond, Stamp_Units::speed);
    tmp.v.y = convert(tmpLegacy.velocity[1], meterPerSecond, Stamp_Units::speed);
    tmp.v.z = convert(tmpLegacy.velocity[2], meterPerSecond, Stamp_Units::speed);
    
  }

}


/// @brief Read particles from the legacy file
/// @param [in] fileId Legacy file
/// @param [in] offset Position to read at
/// @param [in] size Number of particles to read
/// @param [out] particles Particles array to fill
void readLegacyDPDEParticles(InputOutputManager::FileId& fileId, InputOutputManager::Offset& offset, uint64_t size, MPI__Mesoparticle* particles) {

  const auto meterPerSecond = SI_Units_base::meter/SI_Units_base::second;

  // Get an array of legacy-structured particles
  Array<LegacyDPDEParticleIOStruct> array(size);
  InputOutputManager::read(array, fileId, offset);

  // Convert the legacy structured particles into MPI_Particles
  for (uint64_t i=0; i<size; ++i) {

    MPI__Mesoparticle& tmp = particles[i];
    LegacyDPDEParticleIOStruct& tmpLegacy = array[i];

    tmp.id  = tmpLegacy.particleID - 1;
    tmp.ti  = tmpLegacy.particleType;
    tmp.r.x = convert(tmpLegacy.coordinates[0], SI_Units_base::meter, Stamp_Units::length);
    tmp.r.y = convert(tmpLegacy.coordinates[1], SI_Units_base::meter, Stamp_Units::length);
    tmp.r.z = convert(tmpLegacy.coordinates[2], SI_Units_base::meter, Stamp_Units::length);
    tmp.v.x = convert(tmpLegacy.velocity[0], meterPerSecond, Stamp_Units::speed);
    tmp.v.y = convert(tmpLegacy.velocity[1], meterPerSecond, Stamp_Units::speed);
    tmp.v.z = convert(tmpLegacy.velocity[2], meterPerSecond, Stamp_Units::speed);

    tmp.ei  = convert(tmpLegacy.internalEnergy, SI_Units_base::joule, Stamp_Units::energy);
    tmp.progress = tmpLegacy.progress;
    
  }

}
