#pragma once

#include <cstdint>
#include <vector>
#include <cstdlib>

#include <onika/math/basic_types.h>

namespace exaStamp
{

  using namespace exanb;

  struct CellParticleLocalMetrics
  {
    // The size of vectors below equals the number of particles in cell
    // All variables are per-atom quantities

    // Variables contained in grid
    std::vector< double > pe;         // Potential eneregy
    std::vector< Vec3d >  f;          // Force
    std::vector< Vec3d >  v;          // Velocity    
    std::vector< double > q;          // Charge        
    std::vector< Mat3d >  virial;     // Virial
    std::vector< size_t > id;         // Atom id
    std::vector< size_t > idmol;      // Molecule id
    std::vector< size_t > cmol;       // Molecule c(???)
    std::vector< size_t > type;       // Atom type
    std::vector< size_t > typemol;    // Molecule type                 
    std::vector< double > orient;     // Molecule orientation (defined as a double in the grid??)
    std::vector< double > angmom;     // Molecule angular momentum (defined as a double in the grid??)
    std::vector< double > couple;     // Molecule couple (defined as a double in the grid??)

    // Per-atom variables that can be computed (might need to define the computational functor)
    std::vector< double > ke;         // Kinetic energy
    
  };

  using GridParticleLocalMetrics = std::vector< CellParticleLocalMetrics >;

}


