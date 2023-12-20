/// @file
/// @brief Implementations for info class for any grid


#include <algorithm>
#include <vector>

#include "referenceMap.hpp"

#include "grid/anyGridInfo.hpp"
#include "grid/moduler.hpp"

#include "utils/find.hpp"
#include "mpi.h"

/// @brief Constructor
/// @param [in] index Domain index
AnyGridInfo::AnyGridInfo(uint index) 
  : m_domainIndex(index), 
    m_globalCurve(),
    m_globalCurveNoGhost(),
    m_ghostThickness((int) Global::reference.getGhostThickness()),
    m_cellInf(), 
    m_cellSup(),
    m_limInf(),  
    m_limSup(),
    m_numberOfCells(), 
    m_numberOfGhostCells(),
    m_decomposition(nullptr),
    m_neighbors(),
    m_myGlobalIndexes(),
    m_myNeighborDomains() {

  vec3<int> globalOrgn(m_ghostThickness);

  vec3<int> totalNumberOfRealCellsPerDim = Global::domainInfo.getNumberOfCellsPerDim();
  vec3<int> totalNumberOfCellsPerDim     = 2*m_ghostThickness+Global::domainInfo.getNumberOfCellsPerDim();

  // Set the mapping tool of the system + global ghost
  m_globalCurve.set(totalNumberOfCellsPerDim);

	// Link decomposition
  m_decomposition = reinterpret_cast<AnyDecomposition*>(Global::domainInfo.getDecomposition());

  // Set the mapping tool of the system
  m_globalCurveNoGhost.set(totalNumberOfRealCellsPerDim);

  // Fill global indexes
  {
    uint* ptr;
    uint size;
    m_decomposition->myGlobalIndexes(ptr, size);

    /* Raphaël */
    if(size == 0)
    { 
      std::cout<< " This Forest has 0 octrees " << std::endl; 
      MPI_Abort(MPI_COMM_WORLD,2);
    }

    for (uint i=0; i<size; ++i)
      m_myGlobalIndexes.push_back( m_globalCurve.convert(globalOrgn + m_globalCurveNoGhost.convert( ptr[i] )) );
  }


  // Set number of real cells
  // >> init when m_myGlobalIndexes has only real cells !
  m_numberOfCells = m_myGlobalIndexes.size();

  // Get inf/sup domain limits
  // >> init when m_myGlobalIndexes has only real cells !
  {
    auto cellCoords0 = m_globalCurve.convert(m_myGlobalIndexes[0]) - globalOrgn;
    m_cellInf = cellCoords0;
    m_cellSup = cellCoords0;
    for (const auto& elem : m_myGlobalIndexes) {
      auto cellCoordsElem = m_globalCurve.convert(elem) - globalOrgn;
      m_cellInf = auxMin(m_cellInf, cellCoordsElem);
      m_cellSup = auxMax(m_cellSup, cellCoordsElem+1);
    }
    m_limInf = Global::domainInfo.getMinBounds() + m_cellInf * Global::domainInfo.getCellLength();
    m_limSup = Global::domainInfo.getMinBounds() + m_cellSup * Global::domainInfo.getCellLength();
  }

  // Get ghost cells
  // >> only works with m_ghostThickness==1
  if (m_ghostThickness>1) exit(-1);

  vec3<bool> bcs(Global::domainInfo.getBoundaryConditions().x==Configuration<DomainInterface>::PERIODIC,
		 Global::domainInfo.getBoundaryConditions().y==Configuration<DomainInterface>::PERIODIC,
		 Global::domainInfo.getBoundaryConditions().z==Configuration<DomainInterface>::PERIODIC);

  vec3<int> globalCellInf = globalOrgn;
  vec3<int> globalCellSup = globalOrgn + Global::domainInfo.getNumberOfCellsPerDim();

  std::vector<int> ghostCells;
  
  for (uint i=0; i<m_numberOfCells; ++i) {
      
    auto coords = m_globalCurve.convert(m_myGlobalIndexes[i]);

    for (uint nb=0; nb<Neighbor::num_neighbors; ++nb) {

      auto nbrCoords = coords + Neighbor::getCoords(nb);

      auto test = (nbrCoords>=globalCellInf) && (nbrCoords<globalCellSup);

	
      if ( /* in global grid*/ countNumberOf(true, test)==3) {

      	auto nbrIndex = m_globalCurveNoGhost.convert(nbrCoords - globalOrgn);
      	if (/* is not mine */ !m_decomposition->isMine(nbrIndex))
      		ghostCells.push_back( m_globalCurve.convert(globalOrgn + m_globalCurveNoGhost.convert(nbrIndex)) );

      }
      else { /* out */

      	vec3<bool> test2( (!test.x && bcs.x) || test.x,
      			(!test.y && bcs.y) || test.y,
      			(!test.z && bcs.z) || test.z);

      	if (/* periodic in out direction */ countNumberOf(true, test2)==3) {
      		ghostCells.push_back( m_globalCurve.convert(nbrCoords) );
      	}
      }
    } // end for nb
  }

  std::sort(ghostCells.begin(), ghostCells.end());
  auto it = std::unique(ghostCells.begin(), ghostCells.end());
  ghostCells.resize(std::distance(ghostCells.begin(), it));
  
  m_numberOfGhostCells = ghostCells.size();  

  m_myGlobalIndexes.insert(m_myGlobalIndexes.end(), ghostCells.begin(), ghostCells.end());

  // Sort my cells
  std::sort(m_myGlobalIndexes.begin(), m_myGlobalIndexes.end());

}


/// @brief Check if the grid is on a global edge of the system and there are periodic conditions
/// @return For each dimension, true if on the edge and periodic conditions in this dimension
vec3<bool> AnyGridInfo::getGlobalEdges() {

  vec3<bool> bcs(Global::domainInfo.getBoundaryConditions().x==Configuration<DomainInterface>::PERIODIC,
		 Global::domainInfo.getBoundaryConditions().y==Configuration<DomainInterface>::PERIODIC,
		 Global::domainInfo.getBoundaryConditions().z==Configuration<DomainInterface>::PERIODIC);

  vec3<bool> edgeMin = m_cellInf <= zeros();
  vec3<bool> edgeMax = m_cellSup >= Global::domainInfo.getNumberOfCellsPerDim();

  vec3<bool> test(bcs.x && (edgeMin.x || edgeMax.x), 
		  bcs.y && (edgeMin.y || edgeMax.y), 
		  bcs.z && (edgeMin.z || edgeMax.z));

  return test;

}


/// @brief Set the neighbors of specified cells
/// @param [in] cells Cells
/// @param [in] neighbors Cells neighbors
/// @param [in] modI System moduler
void AnyGridInfo::setNeighbors(const Array<uint>& cells, const std::vector< Array<int> >& neighbors, const Moduler<int>& modI) {

  const vec3<int> globalOrgn(m_ghostThickness);

  const vec3<int> globalCellInf = globalOrgn;
  const vec3<int> globalCellSup = globalOrgn + Global::domainInfo.getNumberOfCellsPerDim();

  const vec3<bool> bcs(Global::domainInfo.getBoundaryConditions().x==Configuration<DomainInterface>::PERIODIC,
		       Global::domainInfo.getBoundaryConditions().y==Configuration<DomainInterface>::PERIODIC,
		       Global::domainInfo.getBoundaryConditions().z==Configuration<DomainInterface>::PERIODIC);

  m_myNeighborDomains.resize( m_myGlobalIndexes.size() );

  std::vector<uint> allNeighbors;

  // Loop on cells
  for (uint i=0; i<cells.size(); ++i) {

    const uint cell = cells[i];
    assert( cell>=0 && cell < uint(-3) );
    const vec3<int>  coords2 = m_globalCurve.convert(m_myGlobalIndexes[cell]);
    m_myNeighborDomains[cell].resize(Neighbor::num_neighbors, Neighbor::init);

    // Loop on neighborhoods
    for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr)
    {

    	// Get neighbor cell
   /* 	if( neighbors[cell][nbr]<0 || neighbors[cell][nbr]>=m_myGlobalIndexes.size() )
    	{
    	  std::cerr<<"Warning: neighbors["<<cell<<"]["<<nbr<<"]="<<neighbors[cell][nbr]<<std::endl;
    	}*/
    	
    	if( neighbors[cell][nbr] == -1) continue;
    	
      auto nbrCellGlobalIndex = m_myGlobalIndexes[ neighbors[cell][nbr] ]; // Index on global curve with ghost
      auto nbrCellGlobalCoords = m_globalCurve.convert(nbrCellGlobalIndex);

      // Check if neighbor cell is on global ghost
      auto test = (nbrCellGlobalCoords>=globalCellInf) && (nbrCellGlobalCoords<globalCellSup);

      // If not in global ghost
      if ( /* in global grid*/ countNumberOf(true, test)==3) {

      	// Get neighbor domain
      	auto nbrDomainIndex = m_decomposition->neighborIndex(m_globalCurveNoGhost.convert(coords2 - globalOrgn), nbr);

      	// Add neighbor domain to list if not own
      	if (nbrDomainIndex!=(int)m_domainIndex) {
      		m_myNeighborDomains[cell][nbr] = nbrDomainIndex;
      		allNeighbors.push_back(nbrDomainIndex);

      	}
      	else {

      		m_myNeighborDomains[cell][nbr] = Neighbor::null;
      	}

      }
      else { // If in global ghost

      	// Check boundaries, modulate and add add neighbor domain to list if not own
      	if (countNumberOf(true, test || bcs)==3) {
      		nbrCellGlobalCoords -= globalOrgn;
      		modI.module(nbrCellGlobalCoords);

      		m_myNeighborDomains[cell][nbr] = m_decomposition->neighborIndex(m_globalCurveNoGhost.convert(coords2 - globalOrgn), nbr);
      		allNeighbors.push_back(m_myNeighborDomains[cell][nbr]);
      	}
      	else {
      		m_myNeighborDomains[cell][nbr] = Neighbor::null;
      	}
      }
    }
  }

  // Sort and filter all neighbor domains to get grid neighbors
  std::sort(allNeighbors.begin(), allNeighbors.end());
  auto it = std::unique(allNeighbors.begin(), allNeighbors.end());
  allNeighbors.resize(std::distance(allNeighbors.begin(), it));

  m_neighbors = Array<int>(allNeighbors.size());

  for (uint i=0; i<m_neighbors.size(); ++i) {
    m_neighbors[i] = allNeighbors[i];
  }

}


/// @brief Search a cell in the grid from its coordinates
/// @param [in] coord Coordinates of the cell to find
/// @param [out] found Indicates if cell was found
/// @param [out] idx Index of the cell if found
/// @warning m_myGlobalIndexes must be sorted !
void AnyGridInfo::find(const vec3<int>& coord, bool& found, uint& idx) {
   
   size_t index = 0;

   ::find( (uint) m_globalCurve.convert(m_ghostThickness + coord), 
	   m_myGlobalIndexes.data(), 
	   m_myGlobalIndexes.size(), 
	   found, 
	   index);

   idx = index;

}


/// @brief Search a cell in the grid from its global index (not used)
/// @param [in] globalIndex Global Index of the cell to find
/// @param [out] found Indicates if cell was found
/// @param [out] idx Index of the cell if found
/// @warning m_myGlobalIndexes must be sorted !
void AnyGridInfo::find(const uint globalIndex, bool& found, uint& idx) {
   size_t index = 0;
   ::find( globalIndex, m_myGlobalIndexes.data(), m_myGlobalIndexes.size(), found, index);
   idx = index;
}


/// @brief Debug print in specified flux (do nothing)
/// @param [in,out] flux Print flux
void AnyGridInfo::print(std::ostream& flux) {
}
