#pragma once

#include "cellList/AMR/CellListParticleAMR_Base.hpp"
#include <utils/morton.hpp>





inline void leafCell::checkLeafLocals() {
	CellListBase::checkLocals(size);
}



/// @brief Copy the particles of a local cell into a local ghost cell
/// @tparam resize_velocity Indicates if velocities must be resized
/// (= will be calculated)
/// @tparam resize_fe Indicates if force and potential energy must be resized
/// (= will be calculated)
/// @param [out] to Ghost cell to fill
template <bool resize_velocity, bool resize_fe, class CList>
inline void leafCell::ghostCopyAMR(CList& to, const Correcter& correcter){

  assert(disjoint(to));

  int lastSize = to.size;

  to.resizeInfos(lastSize+size);
  to.resizePositions(lastSize+size);
  to.resizeVelocities(lastSize+size);
  to.resizeForces(lastSize+size);

  auxMemCpy(to.getId()+lastSize, getId(), size);
  auxMemCpy(to.getType()+lastSize, getType(), size);
  auxMemCpy(to.rx.data()+lastSize, getPositionX(), size);
  auxMemCpy(to.ry.data()+lastSize, getPositionY(), size);
  auxMemCpy(to.rz.data()+lastSize, getPositionZ(), size);
  
  to.size +=size;

  assert(to.atomDoublon());

}

/// @brief Create a ghost version of the particles of the cell for the recipient
/// and add them to the leaving ghost array
/// @tparam Q MPI type for exchanged ghost particles
/// @param [in] destDomain Recipient domain
/// @param [in] destCell Recipient cell
/// @param [in,out] leaving Leaving ghost array
template <class Q>
inline void leafCell::ghostCopyAMR(const int destDomain, const vec3<int>& destCell, std::vector< std::tuple<int, Q> >& leaving) const {

	// Associate the recipient domain index with a ghost particle
  auto tmp = std::make_tuple(destDomain, Q());
  // Set the recipient cell in the ghost particle
  std::get<1>(tmp).cell = destCell;

  auto& tmp2 = std::get<1>(tmp).base;

  // For each particle
  for (uint i=0; i< size; ++i) {
  
  	// Set the properties of the ghost particle
    tmp2.id  = getId(i);
    tmp2.ti  = getType(i);
    tmp2.r.x = getPositionX(i);
    tmp2.r.y = getPositionY(i);
    tmp2.r.z = getPositionZ(i);

    // Add the ghost particle to the array
    leaving.push_back(tmp);

  }

}




inline void leafCell::checkMtx() {
  this->mtx.check(size);
}

inline void leafCell::avoidDoublon(vec3<int> & numberOfCellsPerDim, int TotalNumberOfCells) {

	for(uint i = 0; i<neighborCells.size();i++)
		for(uint j=i+1; j<neighborCells.size();j++)
		if(neighborCells[i] == neighborCells[j])
		{
			neighborCells.erase(neighborCells.begin()+j);
			j--;
		}
}


inline void leafCell::resetPtrDaughterCell (size_t i) { 
  assert(i<8);
  daughter[i]=nullptr; 
}


inline void leafCell::clearNeighborCells()
{ 
  neighborCells.clear();
}


template<class Cell> inline void leafCell::addNeighborCell (Cell * cell) 
{
	assert(cell->isleaf == true);
	neighborCells.push_back(cell);
}


/// @brief Fill a buffer with all the particles of the cell, case of a legacy Stamp dump buffer
/// @tparam DumpStruct Dump structure (legacy or legacy for DPDE)
/// @param [in,out] buffer Buffer to fill
/// @param [in] start Where to write this cell in the buffer
// Needed to fill the legacy Stamp V3 dump
template <class DumpStruct> void leafCell::fillBufferAMR(Array<DumpStruct>& buffer, const uint start) {
  
  if(isleaf == false)
  {
    int compteur = start;
    for(int i = 0 ; i < 8 ; i++)
    {
      daughter[i]->fillBufferAMR(buffer, compteur);
      compteur += daughter[i]->size;
    }  
  }
  else
  {
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // crade mais c'est bien le MPI Comm World qui est utilisé
  
	  // Get the speed unit
    const auto meterPerSecond = SI_Units_base::meter/SI_Units_base::second;
    
    const auto& ptr_id = getId();
    const auto& ptr_ti = getType();
    const auto& ptr_rx = getPositionX();
    const auto& ptr_ry = getPositionY();
    const auto& ptr_rz = getPositionZ();
    const auto& ptr_vx = getVelocityX();
    const auto& ptr_vy = getVelocityY();
    const auto& ptr_vz = getVelocityZ();

    // For each particule of the cell
    for (size_t i=0; i<size; ++i) {

    	// Get where to write
    	uint j=i+start;

    	// Write the index and type
      buffer[j].particleID = static_cast<int>(ptr_id[i]); // Id
      buffer[j].particleType = static_cast<int>(ptr_ti[i]); // Type

      // Convert and write the position
      buffer[j].coordinates[0] = convert(ptr_rx[i], Stamp_Units::length, SI_Units_base::meter); // coordinate x
      buffer[j].coordinates[1] = convert(ptr_ry[i], Stamp_Units::length, SI_Units_base::meter); // coordinate y
      buffer[j].coordinates[2] = convert(ptr_rz[i], Stamp_Units::length, SI_Units_base::meter); // coordinate z

      // Convert and write the velocity
      buffer[j].velocity[0] = convert(ptr_vx[i],Stamp_Units::speed, meterPerSecond) ;  // velocity x
      buffer[j].velocity[1] = convert(ptr_vy[i],Stamp_Units::speed, meterPerSecond) ;  // velocity y
      buffer[j].velocity[2] = convert(ptr_vz[i],Stamp_Units::speed, meterPerSecond) ;  // velocity z
       
      buffer[j].johnDoe[0] = granny->neighborList.numberOfNeighbors(i+shift); // number of neighbors
      buffer[j].johnDoe[1] = level; // level
      buffer[j].johnDoe[2] = rank; // mpi process
    }
  }

}
