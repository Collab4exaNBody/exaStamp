/// @file
/// @brief Implementations for the particleWriterVTK


#include <fstream>
#include <string>

#include "globals.hpp"
#include "referenceMap.hpp"

#include "domain/domainInfo.hpp"

#include "io/particleOutputVTK.hpp"

#include "utils/array/array.hpp"
#include "utils/vec3/vec3.hpp"


/// @brief Write the particles
/// @param [in] step Step where the output is written
void ParticleWriterVTK::write(uint step) {

	// Get the file name and open it
  std::string tmp;
  this->setFilename(step, tmp);

  std::ofstream out(tmp.c_str(), std::ios::out);

  // Write the header
  writeHeader(out);
  // Write the positions
  writePositions(out);

  // Write the number of printed objects
  out<< "POINT_DATA " << particles->id.size()+extraSize() << std::endl;

  // Write the indexes
  writeIndexes(out);

  // Write the types if required
  if (particles->writeType)     writeTypes(out);
  // If types required, get the type index/particle name pairings and write the atomic numbers
  if (particles->writeType)     writeNbZ(out);
  // Write the velocities if required
  if (particles->writeVelocity) writeVelocities(out);
  // Write the internal energies if required
  if (particles->writeEint)     writeEints(out);
  // Write the progress varaibles if required
  if (particles->writeProgress) writeProgresses(out);

  // Close the file
  out.close();

}


/// @brief Write the header in the specified file
/// @param [in] flux Output file
void ParticleWriterVTK::writeHeader(std::ofstream& flux) {

  flux<< "# vtk DataFile Version 2.0" << std::endl 
      << "A useless comment" << std::endl 
      << "ASCII"<< std::endl << std::endl;

}


/// @brief Write the particles positions in the specified file
/// @param [in] flux Output file
void ParticleWriterVTK::writePositions(std::ofstream& flux) {

  Array< vec3<double> >& r = particles->r;
  uint size = r.size();

  flux<< "DATASET POLYDATA" << std::endl << "POINTS " << size+extraSize() << " float"  << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << (r[p]*10.0) << std::endl;

  // If required, write the positions of the corners
  if (writeCorners) {

    vec3<double> limits[2];
    limits[0] = Global::domainInfo.getMinBounds();
    limits[1] = Global::domainInfo.getMaxBounds();

    for (int i=0; i<2; ++i) {
      for (int j=0; j<2; ++j) {
      	for (int k=0; k<2; ++k) {
      		vec3<double> R(limits[i].x*10.0, limits[j].y*10.0, limits[k].z*10.0);
      		flux << R << std::endl;
      	}
      }
    }

  }

}


/// @brief Write the particles indexes in the specified file
/// @param [in] flux Output file
void ParticleWriterVTK::writeIndexes(std::ofstream& flux) {

  Array<uint>& id = particles->id;
  uint size = id.size();

  flux<< "SCALARS index float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << id[p] << std::endl;

  // If corners required, write the -1 as index
  if (writeCorners) {
    for (uint p=0; p<extraSize(); ++p) 
      flux << -1 << std::endl;
  }

}


/// @brief Write the particles types in the specified file
/// @param [in] flux Output file
void ParticleWriterVTK::writeTypes(std::ofstream& flux) {

  Array<uint8_t>& type = particles->type;
  uint size = type.size();

  flux<< "SCALARS type float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << static_cast<int>(type[p]) << std::endl;

  // If corners required, write the -1 as type
  if (writeCorners) {
    for (uint p=0; p<extraSize(); ++p) 
      flux << -1 << std::endl;
  }

}


/// @brief Write the atomic numbers in the specified file
/// @param [in] flux Output file
void ParticleWriterVTK::writeNbZ(std::ofstream& flux) {

  Array<uint8_t>& type = particles->type;
  uint size = type.size();

  flux<< "SCALARS NbZ integer 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;
  // Get the number of particle types
  uint8_t nTypes = Global::reference.getNumberOfTypes();
  std::vector<uint> nbZ(nTypes,0);
  // For each type write the index, name and mass, plus atomic number and charge if it's an atom type in the standard error flux
  // and get the atomic number
  for(int i=0;i<nTypes;i++){
     TypeParticle * tp =Global::reference.find(i);
     TypeAtom *ta = dynamic_cast<TypeAtom*>(tp);
     if (ta){
        std::cerr<< " Atom Type " << i << " name=" << tp->getName() << " mass=" << tp->getMass() << " AtomicNumber=" << ta->getAtomicNumber() <<std::endl;
         nbZ[i]=ta->getAtomicNumber();
     }else{
        std::cerr<< " Type " << i << " name=" << tp->getName() << " mass=" << tp->getMass() <<std::endl;
     }
  }

  // Write the atomic number for each particle
  for (uint p=0; p<size; ++p)
    flux << nbZ[(static_cast<int>(type[p]))] << std::endl;

  // If corners required, write the -1 as atomic number
  if (writeCorners) {
    for (uint p=0; p<extraSize(); ++p)
      flux << -1 << std::endl;
  }

}


/// @brief Write the particles velocities in the specified file
/// @param [in] flux Output file
void ParticleWriterVTK::writeVelocities(std::ofstream& flux) {

  Array< vec3<double> >& v = particles->v;
  uint size = v.size();

  flux<< "SCALARS velocity float 3" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << v[p] << std::endl;

  // If corners required, write the 0. as velocity
  if (writeCorners) {
    for (uint p=0; p<extraSize(); ++p) 
      flux << vec3<double>(0.) << std::endl;
  }

}


/// @brief Write the particles internal energies in the specified file
/// @param [in] flux Output file
void ParticleWriterVTK::writeEints(std::ofstream& flux) {

  Array< double >& ei = particles->ei;
  uint size = ei.size();

  flux<< "SCALARS eint float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << ei[p] << std::endl;

  if (writeCorners) {
    for (uint p=0; p<extraSize(); ++p) 
      flux << 0. << std::endl;
  }

}


/// @brief Write the particles progres variables in the specified file
/// @param [in] flux Output file
void ParticleWriterVTK::writeProgresses(std::ofstream& flux) {

  Array< double >& progress = particles->progress;
  uint size = progress.size();

  flux<< "SCALARS progress float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << progress[p] << std::endl;

  if (writeCorners) {
    for (uint p=0; p<extraSize(); ++p) 
      flux << 0. << std::endl;
  }

}
