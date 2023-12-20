/// @file
/// @brief Implementations for info class for rectilinear grid


#include "globals.hpp"
#include "referenceMap.hpp"

#include "domain/domainInfo.hpp"

#include "grid/rectilinearGridInfo.hpp"

#include "utils/neighbor.hpp"


/// @brief Constructor
/// @param [in] index Domain index
RectilinearGridInfo::RectilinearGridInfo(uint index)
  : decomposition(nullptr),
    numberOfNeighbors(0), 
    neighbors(Neighbor::num_neighbors, Neighbor::init) {

	// Link decomposition
  decomposition = reinterpret_cast<RectilinearDecomposition*>(Global::domainInfo.getDecomposition());

  // Get origin and number of cell per dimension from decompostion
  origin              = decomposition->origins(index);
  numberOfCellsPerDim = decomposition->sizes  (index);

  // Get ghost thickness from globals
  ghostThickness = (int) Global::reference.getGhostThickness();

  // Store cell length and global min bound
  invCellLength = 1./Global::domainInfo.getCellLength();
  globalMin     = Global::domainInfo.getMinBounds();

  // Calculate minBounds, maxBounds and extension
  minBounds = globalMin + Global::domainInfo.getCellLength() *  origin;
  extension = Global::domainInfo.getCellLength() *  numberOfCellsPerDim;
  maxBounds = minBounds + extension;



  // Calculate number of cells
  numberOfCells = product(numberOfCellsPerDim);

  // Set the mapping tool with grid dimensions
  curve.set(numberOfCellsPerDim+2*ghostThickness);

  // Get neighbors from decomposition
  for (uint nbr=0; nbr<neighbors.size(); ++nbr) {
	neighbors[nbr] = decomposition->neighbors(index, nbr);
	if (neighbors[nbr]>=0) numberOfNeighbors++;
  }

  // Calculate local/global coordinates shift
  m_shift = ghostThickness - origin;

}


/// @brief Check if the grid is on a global edge of the system and there are periodic conditions
/// @return For each dimension, true if on the edge and periodic conditions in this dimension
vec3<bool> RectilinearGridInfo::getGlobalEdges() {

  vec3<bool> bcs(Global::domainInfo.getBoundaryConditions().x==Configuration<DomainInterface>::PERIODIC, 
		 Global::domainInfo.getBoundaryConditions().y==Configuration<DomainInterface>::PERIODIC, 
		 Global::domainInfo.getBoundaryConditions().z==Configuration<DomainInterface>::PERIODIC);

  vec3<bool> edgeMin = origin <= zeros();
  vec3<bool> edgeMax = origin + numberOfCellsPerDim >= Global::domainInfo.getNumberOfCellsPerDim();

  vec3<bool> test(bcs.x && (edgeMin.x || edgeMax.x), 
		  bcs.y && (edgeMin.y || edgeMax.y), 
		  bcs.z && (edgeMin.z || edgeMax.z));

  return test;

}


/// @brief Get the index of a neighbor domain after checking that the considered cell is on domain edge
/// @param [in] mov Coordinates of the neighbor domain
/// @param [in] cellCoords Coordianates of the cell
/// @return Index of the neighbor domain or null if cell is not on domain edge
int RectilinearGridInfo::getNeighborDomain(const vec3<int>& mov, const vec3<int>& cellCoords) {

  static const vec3<int> max = origin + (int) Global::reference.getGhostThickness();
  static const vec3<int> min = origin + numberOfCellsPerDim - (int) Global::reference.getGhostThickness();

  vec3<int> domainMov(0);

  for (int dim=0; dim<VEC3_NDIMS; ++dim) {
    if      (mov[dim]==-1 && cellCoords[dim] <  max[dim]) domainMov[dim]=-1;
    else if (mov[dim]==+1 && cellCoords[dim] >= min[dim]) domainMov[dim]=+1;
  }

  int index = Neighbor::getIndex(domainMov);

  if (index>=0) return neighbors[index];
  else          return index;

  return 0;

}

/// @brief Debug print in specified flux (do nothing)
/// @param [in,out] flux Print flux
void RectilinearGridInfo::print(std::ostream& flux) {
}
