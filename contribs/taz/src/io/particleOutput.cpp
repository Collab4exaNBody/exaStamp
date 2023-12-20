/// @file
/// @brief Implementation of tools to print the particle output in DAT and VTK formats


#include "io/particleOutput.hpp"
#include "parallel/commManager.hpp"


/// @brief Gather particle output from all nodes on one node
/// @param [in] comm Communication manager
/// @param [in] sendBuffer Particles on each node
/// @param [in] counts Number of particle to get on each node
/// @param [in] disps Position where to put the data from each node
void ParticleOutput::gather(CommManager* comm, ParticleOutput* sendBuffer, Array<int>& counts, Array<int>& disps) {

  // Gather indexes
  comm->gatherV(sendBuffer->id, this->id, counts, disps);
  
  // Gather types if required
  if (writeType)
    comm->gatherV(sendBuffer->type, this->type, counts, disps);
  
  // Gather positions
  comm->gatherV(sendBuffer->r, this->r, counts, disps);
    
  // Gather velocities if required
  if (writeVelocity)
    comm->gatherV(sendBuffer->v, this->v, counts, disps);

  if (writeEint)
    comm->gatherV(sendBuffer->ei, this->ei, counts, disps);

  if (writeProgress)
    comm->gatherV(sendBuffer->progress, this->progress, counts, disps);
}


/// @brief Create the name of the output file
/// @param [in] step Step where the output is written
/// @param [out] filename Output file name
void ParticleWriter::setFilename(uint step, std::string& filename) {

  std::string tmp = std::to_string(step);

  while (tmp.size()<LOG_MAX_STEPS) tmp.insert(0, "0");

  filename = rootName + "_" + tmp + extension;

}
