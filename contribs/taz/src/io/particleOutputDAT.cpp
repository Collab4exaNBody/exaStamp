/// @file
/// @brief Implementations for the particleWriterDAT


#include <fstream>
#include <iomanip>
#include <ostream>
#include <string>

#include "io/particleOutputDAT.hpp"


/// @brief Write the particles
/// @param [in] step Step where the output is written
void ParticleWriterDAT::write(uint step) {

  // Get the name of the file
  std::string tmp;
  this->setFilename(step, tmp);

  // Open the file
  std::ofstream out(tmp.c_str(), std::ios::out);

  // Header
  out<< "# INDEX (1)         ";
  if (particles->writeType)
    out << "TYPE (1)                              ";
  out << "POSITION (3)                              ";
  if (particles->writeVelocity)
    out << "VELOCITY (3)                              ";
  if (particles->writeEint)
    out << "INT. ENERGY (1)                      ";
  if (particles->writeProgress)
    out << "PROGRESS (1) ";

  out << std::endl;

  // Print all the particles
  for (uint p=0; p<particles->id.size(); ++p) {

      out << std::setw(10) << (int) particles->id[p] << "      ";

      if (particles->writeType)
          out << std::setw(3) << static_cast<int>(particles->type[p]) << "    ";
      
      out << std::scientific << std::setprecision(16) << std::setw(24) << particles->r[p].x << " "
          << std::scientific << std::setprecision(16) << std::setw(24) << particles->r[p].y << " "
          << std::scientific << std::setprecision(16) << std::setw(24) << particles->r[p].z << "    ";

      if (particles->writeVelocity) { 
	out << std::scientific << std::setprecision(16) << std::setw(24) << particles->v[p].x << " "
	    << std::scientific << std::setprecision(16) << std::setw(24) << particles->v[p].y << " "
	    << std::scientific << std::setprecision(16) << std::setw(24) << particles->v[p].z << " ";
      }

      if (particles->writeEint)
	out << std::scientific << std::setprecision(16) << std::setw(24) << particles->ei[p] << " ";

      if (particles->writeProgress)
	out << std::scientific << std::setprecision(16) << std::setw(24) << particles->progress[p] << " ";

      out << std::endl;
      
  }

  // Close the file
  out.close();
  
}
