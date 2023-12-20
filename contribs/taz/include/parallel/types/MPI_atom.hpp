/// @file 
/// @brief Definition of types to exchange atoms

#ifndef __MPI_ATOM_HPP_INCLUDED
#define __MPI_ATOM_HPP_INCLUDED


#include "parallel/types/MPI_particle.hpp"
#include "parallel/types/MPI_ghost.hpp"


/// @brief Type to exchange an atom with its velocity (same as for a particle)
typedef MPI__Particle      MPI__Atom;
/// @brief Type to communicate the atoms of the ghost (same as for a particle)
typedef MPI__Ghost<MPI__ParticleBase> MPI__AtomGhost;
/// @brief Short name for MPI__Atom
typedef MPI__Atom Atom;

#endif // __MPI_ATOM_HPP_INCLUDED
