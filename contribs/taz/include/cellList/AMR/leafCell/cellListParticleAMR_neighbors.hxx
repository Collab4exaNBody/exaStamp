#pragma once

#include "simdAMR/neighbors.hpp"

#define SIZEVECTOR 16



/// @brief Sort neighbor list of a given particle
///
/// This function sorts neighbor list of particle index, so that particles
/// of the same type are contiguous in memory
inline void leafCell::sortNeighbors() {

  this->neighborList.sort( [&] (const NeighborList_base::nbr_id& a,
				const NeighborList_base::nbr_id& b) -> bool {

			     return std::get<NeighborList_base::TYPE>(a) < std::get<NeighborList_base::TYPE>(b);

			   });

}

inline void leafCell::neighborListAddNeighbor(const uint16_t index, const nbr_id& nbr)
{
  granny->neighborList.addNeighbor(index+shift,nbr);
}


inline void leafCell::neighborListGetNeighbors(const uint16_t index, const uint8_t type, nbr_id*& start, uint& n)
{
  granny->neighborList.getNeighbors(index+shift, type, start, n);
}  


inline leafCell** leafCell::neighborListGetLeafCellNeighbors(const uint16_t index)
{
	if(isleaf==1)
	{
		return neighborCells.data();
	}
	else
	{
		for(size_t i=0; i<8 ; i++)
		{
			if(index >= daughter[i]->shift && index < daughter[i]->shift+daughter[i]->size)
				return daughter[i]->neighborListGetLeafCellNeighbors(index);
		}
	}

	std::cout << " on ne devrait jamais passer ici " << std::endl;
	std::abort();
	return nullptr;
}  

/// @brief Build neighbor list between the particles of the cell and those of a neighbor cell
/// @param [in] cellIndex Index of the neighbor cell
/// @param [in] cell Reference to the neighbor cell
/// @param [in] mask Projection mask (distance between a particle and a cell, used to prevent useless neighbor computation)
/// @param [in] symmetrize In case of symmetrized force computation, indicates that the neighbor search must be done for these cells
void leafCell::makeNeighborsAMR(const size_t cellIndex, self_type& cell, const bool* mask, const bool symmetrize, const bool perBlock) {

	int sizeOfBlock=1-1;
	if(perBlock) sizeOfBlock = SIZEVECTOR-1;

	if (symmetrize) 
	{
		// Create the identification to put in neighbor list
		NeighborList_base::nbr_id nbr;
		std::get<NeighborList_base::CELL>(nbr) = cellIndex;

		uint8_t * ptr_ti = getType();
		double * ptr_rx = getPositionX();
		double * ptr_ry = getPositionY();
		double * ptr_rz = getPositionZ();


		uint8_t * cell_ptr_ti = cell.getType();
		double * cell_ptr_rx = cell.getPositionX();
		double * cell_ptr_ry = cell.getPositionY();
		double * cell_ptr_rz = cell.getPositionZ();


		//assert(disjoint(cell));

		// Lambda to perform neighboring test
		auto isNeighbor = [&] (const size_t i, const size_t j) -> bool {

			assert(ptr_rx[i] != cell_ptr_rx[j] || 
				ptr_ry[i] != cell_ptr_ry[j] || 
				ptr_rz[i] != cell_ptr_rz[j] 
			);

			return (ptr_rx[i]-cell_ptr_rx[j])*(ptr_rx[i]-cell_ptr_rx[j]) + 
			(ptr_ry[i]-cell_ptr_ry[j])*(ptr_ry[i]-cell_ptr_ry[j]) + 
			(ptr_rz[i]-cell_ptr_rz[j])*(ptr_rz[i]-cell_ptr_rz[j]) < 
			Global::reference.getrVerlet2(ptr_ti[i], cell_ptr_ti[j]);
		};


		// Loop on the particles of the cell
		for (uint i=0; i<size; ++i) {
			// Check mask
			if (mask[i]) continue; 

			// Get the cutoff radii between the current particle and those of the neighbor cell
			for (size_t j=0; j<cell.size; ++j) {
				if (isNeighbor(i,j)) {
					std::get<NeighborList_base::INDX>(nbr) = j;
					std::get<NeighborList_base::TYPE>(nbr) = cell_ptr_ti[j];
					neighborListAddNeighbor(i, nbr);
					j+=sizeOfBlock;
				}
			}
		}
	}
}



/// @brief Build neighbor list between the particles on the cell
/// @param [in] myIndex Index of the cell
/// @param [in] symmetrize Indicates if the force computation is symmetrized
inline void leafCell::makeNeighborsAMR_self(const size_t myIndex, const bool symmetrize, const bool perBlock)
{
	int sizeOfBlock=1-1;
	if(perBlock) sizeOfBlock = SIZEVECTOR-1;

	// Create the identification to put in neighbor list
	NeighborList_base::nbr_id nbr;
	std::get<NeighborList_base::CELL>(nbr) = myIndex;

	// Pointers for fast access
	auto ptr_ti = getType();
	auto ptr_rx = getPositionX();
	auto ptr_ry = getPositionY();
	auto ptr_rz = getPositionZ();

	// Lambda to perform neighboring test
	auto isNeighbor = [&] (const size_t i, const size_t j) -> bool {

		assert(i!=j);

		assert(ptr_rx[i] != ptr_rx[j] || ptr_ry[i] != ptr_ry[j] || ptr_rz[i] != ptr_rz[j]);

		return (ptr_rx[i]-ptr_rx[j])*(ptr_rx[i]-ptr_rx[j]) + 
		(ptr_ry[i]-ptr_ry[j])*(ptr_ry[i]-ptr_ry[j]) + 
		(ptr_rz[i]-ptr_rz[j])*(ptr_rz[i]-ptr_rz[j]) < 
		Global::reference.getrVerlet2(ptr_ti[i], ptr_ti[j]);
	};

	// If symmetrized force computation
	if (symmetrize) 
	{
		for (uint i=0; i<size; ++i) 
		{
			// Build neighbor list
			for (uint j=i+1; j<size; ++j) {
				if (isNeighbor(i,j)) {
					std::get<NeighborList_base::INDX>(nbr) = j;
					std::get<NeighborList_base::TYPE>(nbr) = ptr_ti[j];
					neighborListAddNeighbor(i, nbr);
					j+=sizeOfBlock;
				}
			}
		}
	}
	// If not symmetrized force computation
	else 
	{
		// For each particle
		for (uint i=0; i<size; ++i) 
		{
			// For each particle
			for (uint j=0; j<size; ++j) 
			{
				// Test if neighbor
				if ( i!=j && isNeighbor(i, j) ) 
				{
					// Add to list if it is
					std::get<NeighborList_base::INDX>(nbr) = j;
					std::get<NeighborList_base::TYPE>(nbr) = ptr_ti[j];
					neighborListAddNeighbor(i, nbr);
					j+=sizeOfBlock;
				}
			}
		}
	}
}


/// @brief Compute distance between the particles of the cell and a neighbor cell to prevent useless neighbor computation
/// @tparam elem_chunk Base chunk for the arrays
/// @param [in] coords Coordinates of the cell
/// @param [in] dcell Coordinates of the neighbor cell
/// @param [out] mask Boolean result for each particle (is or isn't near enough of the neighbor cell to consider neighbor computation)
/// @param [in] correct Indicates if the neighbor cell is across the system edge and distances need to be corrected
/// @param [in] correcter Tool to correct distances
inline void leafCell::getProjAMR(const vec3<int>& coords, const vec3<int>& dcell, bool mask[]) {

  // Global variables used several times
  static const int m = 1 << sizeof(int);
  static const double maxRcut2 = auxSq(Global::reference.getMaxRcut());
  const vec3<double> cLength = Global::domainInfo.getCellLength()/(1<<level);
  static const vec3<double>& minBounds  = Global::domainInfo.getMinBounds();

  double prx[size];
  double pry[size];
  double prz[size];

  const auto& ptr_rx = getPositionX();
  const auto& ptr_ry = getPositionY();
  const auto& ptr_rz = getPositionZ();
  
  int c[3] = { coords.x, coords.y, coords.z };

  // Loop on particles of the cell
  for (uint i=0; i<size; ++i) {

    prx[i] = 0.;
    pry[i] = 0.;
    prz[i] = 0.;

    // Distance between the particle and the neighbor cell on the x axis
    if (dcell.x & 0x01 /* true if dcell.x != 0 */) {
      prx[i] = minBounds.x + c[0] * cLength.x - ptr_rx[i];
      if (dcell.x & m /* true if dcell.x==-1 */) 
      	prx[i] -= cLength.x;
    }
    
    // Distance between the particle and the neighbor cell on the y axis
    if (dcell.y & 0x01 /* true if dcell.y != 0 */) {
      pry[i] = minBounds.y + c[1] * cLength.y - ptr_ry[i];
      if (dcell.y & m /* true if dcell.y==-1 */) 
      	pry[i] -= cLength.y;
    }
    
    // Distance between the particle and the neighbor cell on the z axis
    if (dcell.z & 0x01 /* true if dcell.z != 0 */) {
      prz[i] = minBounds.z + c[2] * cLength.z - ptr_rz[i];
      if (dcell.z & m /* true if dcell.z==-1 */) 
      	prz[i] -= cLength.z;
    }
      
  }
	
	#pragma omp simd
	//#pragma vector aligned
	for (uint i=0; i<size; ++i)    
           mask[i] = prx[i]*prx[i] + pry[i]*pry[i] + prz[i]*prz[i] > maxRcut2;
}

