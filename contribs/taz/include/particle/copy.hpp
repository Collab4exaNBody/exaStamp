/// @file
/// @brief Tools to copy exchanged particles onto particles (obsolete but could be used someday)

#ifndef __COPY_HPP_INCLUDED
#define __COPY_HPP_INCLUDED


#include "parallel/types/MPI_particle.hpp"


/// @brief Copy an MPI particle onto a particle (case of and MPI_ParticleBase)
/// @tparam N Not used
/// @param [in] from MPI particle to copy
/// @param [out] to Recipient for the copy
template <const uint N>
inline void copy(const MPI__ParticleBase& from, Particle& to) {
  
  to.id() = from.id;
  to.ti() = from.ti;

  to.rx() = from.r.x;
  to.ry() = from.r.y;
  to.rz() = from.r.z;

}


/// @brief Copy an MPI particle onto a particle (case of and MPI_Particle)
/// @tparam N Not used
/// @param [in] from MPI particle to copy
/// @param [out] to Recipient for the copy
template <const uint N>
inline void copy(const MPI__Particle& from, Particle& to) {

  copy(static_cast<const MPI__ParticleBase&>(from), to);

  to.vx() = from.v.x;
  to.vy() = from.v.y;
  to.vz() = from.v.z;

}


/// @brief Copy an MPI particle onto a particle (case of and MPI_Ghost)
/// @tparam P Particle class in the ghost
/// @tparam N Not used
/// @param [in] from MPI particle to copy
/// @param [out] to Recipient for the copy
template <class P, const uint N>
inline void copy(const MPI__Ghost<P>& from, Particle& to) {
  copy(from.base, to);
}

#endif //  __COPY_HPP_INCLUDED
