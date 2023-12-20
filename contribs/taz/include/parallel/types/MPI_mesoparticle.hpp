/// @file 
/// @brief Definition of MPI_Mesoparticle struct

#ifndef __MPI_MESOPARTICLE_HPP_INCLUDED
#define __MPI_MESOPARTICLE_HPP_INCLUDED


#include "parallel/types/MPI_particle.hpp"

/// @struct MPI__Mesoparticle
/// @brief MPI type for mesoparticles
struct MPI__Mesoparticle : public MPI__Particle {

  /// @brief Default constructor
  MPI__Mesoparticle() : MPI__Particle(), ei(), progress(0) {}

  /// @brief Destructor (nothing to do)
  ~MPI__Mesoparticle() {}

  /// @brief Copy constructor
  /// @param [in] particle Mesoparticle to copy
  MPI__Mesoparticle(const MPI__Mesoparticle& particle)
    : MPI__Particle(particle), ei(particle.ei), progress(particle.progress) {}

  /// @brief Assignment operator
  /// @param [in] particle Mesoparticle to copy
  /// @return Copy of particle
  MPI__Mesoparticle& operator = (const MPI__Mesoparticle& particle) {
    id = particle.id;
    ti = particle.ti;
    r  = particle.r;
    v  = particle.v;
    ei = particle.ei;
    progress = particle.progress;
    return *this;
  }

  /// @brief Define of the base to this particle ghost
  typedef MPI__Mesoparticle Base;

  double ei; ///< Internal energy
  double progress; ///< Progress variable for chemical reactions
  
};

/// @brief Short name for MPI__Mesoparticle
typedef MPI__Mesoparticle Mesoparticle;


#endif // __MPI_MESOPARTICLE_HPP_INCLUDED
