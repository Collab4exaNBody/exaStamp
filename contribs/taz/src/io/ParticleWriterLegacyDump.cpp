/// @file
/// @brief Implementations for the particleWriterLegacyDump


#include "io/ParticleWriterLegacyDump.hpp"


/// @brief Get the name of the legacy dump file
/// @param [in] step Step where the dump file is written
void  ParticleWriterLegacyDump::setFileName(uint step){
  std::string tmp = std::to_string(step);
  while (tmp.size() < 9) tmp.insert(0,"0");
  legacyDumpName="StampV3prot_"+tmp+".MpiIO";
}


/// @brief Open the dump file
void  ParticleWriterLegacyDump::open(){
  // Collective open
  m_inputOutputManager->openIOFile(legacyDumpName,"w",mpiIOFileId);
}


/// @brief Write the header
/// @param [in] header Header to write
void ParticleWriterLegacyDump:: writeHeader(Array<LegacyHeaderIOStruct>& header){
  // Only the master write the header
  localOffset=0;
  if (comm->getRank()==Global::masterNode)  // idem isMaster (class Node)
    InputOutputManager::write(header,mpiIOFileId,localOffset);

  // Set offset to write after the header
  localOffset=sizeof(LegacyHeaderIOStruct);
}


/// @brief Write the particles
/// @param [in] particlesArray Particles to write
void  ParticleWriterLegacyDump::writeParticles(Array<LegacyParticleIOStruct>& particlesArray){
  // Share the number of particles on each node to set the offsets
  Array<int> localNodeParticlesNumbers(1, particlesArray.size());
  Array<int> arrayOfAllLocalNodeParticlesNumbers(comm->getNumberOfNodes(),0);
  comm->gather(localNodeParticlesNumbers,arrayOfAllLocalNodeParticlesNumbers);
  comm->broadcast(arrayOfAllLocalNodeParticlesNumbers);

  for (int i=0; i<comm->getRank(); ++i)
    localOffset += arrayOfAllLocalNodeParticlesNumbers[i]*sizeof(LegacyParticleIOStruct);
  // Write particles in the legacy Dump File
  InputOutputManager::write(particlesArray,mpiIOFileId,localOffset);
}


/// @brief Write the particles
/// @param [in] particlesArray Particles to write
void  ParticleWriterLegacyDump::writeParticles(Array<LegacyDPDEParticleIOStruct>& particlesArray){
  // Share the number of particles on each node to set the offsets
  Array<int> localNodeParticlesNumbers(1, particlesArray.size());
  Array<int> arrayOfAllLocalNodeParticlesNumbers(comm->getNumberOfNodes(),0);
  comm->gather(localNodeParticlesNumbers,arrayOfAllLocalNodeParticlesNumbers);
  comm->broadcast(arrayOfAllLocalNodeParticlesNumbers);

  for (int i=0; i<comm->getRank(); ++i)
    localOffset += arrayOfAllLocalNodeParticlesNumbers[i]*sizeof(LegacyDPDEParticleIOStruct);
  // Write particles in the legacy Dump File
  InputOutputManager::write(particlesArray,mpiIOFileId,localOffset);
}


/// @brief Close the dump file
void  ParticleWriterLegacyDump::close(){
  // Collective Close
  m_inputOutputManager->closeIOFile(mpiIOFileId);
}


/// @brief Write the dump file for the current step
/// @param [in] step Step where the dump file is written
/// @param [in] header Header to write
/// @param [in] particlesArray Particles to write
void ParticleWriterLegacyDump::write(uint step,Array<LegacyHeaderIOStruct>& header,Array<LegacyParticleIOStruct>& particlesArray) {
  // Set file name
  this->setFileName(step);
  // Open file
  this->open();
  // Write header
  this->writeHeader(header);
  // Write particles
  this->writeParticles(particlesArray);
  // Close file
  this->close();

}


/// @brief Write the dump file for the current step
/// @param [in] step Step where the dump file is written
/// @param [in] header Header to write
/// @param [in] particlesArray Particles to write
void ParticleWriterLegacyDump::write(uint step,Array<LegacyHeaderIOStruct>& header,Array<LegacyDPDEParticleIOStruct>& particlesArray) {
  // Set file name
  this->setFileName(step);
  // Open file
  this->open();
  // Write header
  this->writeHeader(header);
  // Write particles
  this->writeParticles(particlesArray);
  // Close file
  this->close();

}
