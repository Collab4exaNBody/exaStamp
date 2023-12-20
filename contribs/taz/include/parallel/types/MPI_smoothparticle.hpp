/// @file 
/// @brief Definition of MPI_SmoothParticle struct

#ifndef __MPI_SMOOTH_PARTICLE_HPP_INCLUDED
#define __MPI_SMOOTH_PARTICLE_HPP_INCLUDED


#include "parallel/types/MPI_mesoparticle.hpp"

/// @struct MPI__SmoothParticle
/// @brief MPI type for smoothparticles
struct MPI__SmoothParticle : public MPI__Mesoparticle {

  /// @brief Default constructor
  MPI__SmoothParticle() : MPI__Mesoparticle(), rho(0) {}

  /// @brief Destructor (nothing to do)
  ~MPI__SmoothParticle() {}

  /// @brief Copy constructor
  /// @param [in] particle Smoothparticle to copy
  MPI__SmoothParticle(const MPI__SmoothParticle& particle)
    : MPI__Mesoparticle(particle), rho(particle.rho) {}

  /// @brief Assignment operator
  /// @param [in] particle Smoothparticle to copy
  /// @return Copy of the particle
  MPI__SmoothParticle& operator = (const MPI__SmoothParticle& particle) {
    id = particle.id;
    ti = particle.ti;
    r  = particle.r;
    v  = particle.v;
    ei = particle.ei;
    progress = particle.progress;
    rho = particle.rho;
    return *this;
  }

  /// @brief Define of the base to this particle ghost
  typedef MPI__SmoothParticle Base;

  double rho; ///< Density

};

/// @brief short name for MPI_SmoothParticle
typedef MPI__SmoothParticle SmoothParticle;


#endif // __MPI_SMOOTH_PARTICLE_HPP_INCLUDED
