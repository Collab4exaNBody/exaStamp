#pragma once

inline void Octree::fillVerlet()
{
	// The position of atoms is stored at the iteration of the update of the verlet lists.
	// This position is used to test if an atom is moved by more than 1/2 of Rverlet.
	this->neighborList.check(this->size);
	this->neighborList.resize_store(this->size);  
	this->neighborList.addCoord( this->size, rx, ry, rz);
}

/// @brief Compute the Verlet List condition
/// @param [in] rVerlet 1/2 of raduis of Verlet 
/// @return Indicates if verlet list should be updated
inline bool Octree::checkVerlet(double rVerlet)
{
  return this->neighborList.checkVerlet(getPositionX(), getPositionY(), getPositionZ(), rVerlet, this->size);
} 
