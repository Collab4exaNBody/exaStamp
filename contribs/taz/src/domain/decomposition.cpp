/// @file 
/// @brief Implementation of Decomposition subclasses


#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

#include "globals.hpp"
#include "referenceMap.hpp"

#include "domain/decomposition.hpp"
#include "domain/domainInfo.hpp"

#include "parallel/thread/thread.hpp"

#include "utils/curve.hpp"
#include "utils/find.hpp"
#include "utils/neighbor.hpp"
#include "utils/split.hpp"


/// @brief Constructor
/// @param [in] configuration Domain configuration
RectilinearDecomposition::RectilinearDecomposition(Configuration<DomainInterface>& configuration)
  : Decomposition(product(configuration.numberOfDomainsPerDim)), 
    m_numberOfDomainsPerDim(), 
    m_origins(numberOfDomains), m_sizes(numberOfDomains), 
    m_neighbors(numberOfDomains, Array<int>(Neighbor::num_neighbors, Neighbor::init)) {

  const auto& globalNumberOfCellsPerDim =Global::domainInfo.getNumberOfCellsPerDim();

  m_numberOfDomainsPerDim.x = (int) configuration.numberOfDomainsPerDim.x;
  m_numberOfDomainsPerDim.y = (int) configuration.numberOfDomainsPerDim.y;
  m_numberOfDomainsPerDim.z = (int) configuration.numberOfDomainsPerDim.z;

  const int min = 1+2*Global::reference.getGhostThickness();

  const vec3<bool> criticalTest = min * m_numberOfDomainsPerDim > globalNumberOfCellsPerDim;
  if (countNumberOf(true, criticalTest)>0) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'RectilinearDecomposition::RectilinearDecomposition(Configuration<DomainInterface>&)' : Too many procs for domain size or potential. STOP." 
	     << std::endl;
    exit(1);
  }

  // Get a function to convert 3D index into 1D index
  LCurve<int> crv(m_numberOfDomainsPerDim);

  // Split in a regular way
  for (uint domain=0; domain<numberOfDomains; ++domain) {

    vec3<int> domainCoords = crv.convert(domain);
    
    // split in each dimension
    split<int>(globalNumberOfCellsPerDim.x, m_numberOfDomainsPerDim.x, domainCoords.x, m_origins[domain].x, m_sizes[domain].x);
    split<int>(globalNumberOfCellsPerDim.y, m_numberOfDomainsPerDim.y, domainCoords.y, m_origins[domain].y, m_sizes[domain].y);
    split<int>(globalNumberOfCellsPerDim.z, m_numberOfDomainsPerDim.z, domainCoords.z, m_origins[domain].z, m_sizes[domain].z);

    // compute neighboring domains
    for (uint nb=0; nb<Neighbor::num_neighbors; ++nb) {

      // coords (3D) of neighboring domain
      vec3<int> nbCoords = domainCoords + Neighbor::getCoords(nb);

      // case of an edge domain, regarding boundary conditions
      int flag = 0;
      for (int dim=0; dim<VEC3_NDIMS; ++dim) {

	if (nbCoords[dim]<0 || nbCoords[dim]>=m_numberOfDomainsPerDim[dim]) {
	  switch (configuration.boundaryConditions[dim]) {

	  case Configuration<DomainInterface>::FREE:
	  case Configuration<DomainInterface>::WALL:
	  case Configuration<DomainInterface>::FREE_WALL:
	  case Configuration<DomainInterface>::WALL_FREE:
	    flag++;
	    break;
	    
	  case Configuration<DomainInterface>::PERIODIC:
	    nbCoords[dim] = auxMod(nbCoords[dim], m_numberOfDomainsPerDim[dim]);
	    break;
	    
	  }
	}
      }

      // if not periodic and flag is set, no neighbor
      m_neighbors[domain][nb] = (flag==0) ? crv.convert(nbCoords) : Neighbor::null;
	
    }
            
  }

}


/// @brief Destructor (nothing to do)
RectilinearDecomposition::~RectilinearDecomposition() {
}


/// @brief Print the type of decomposition and number of domains (per dimension and total) in specified flux
/// @param [in,out] flux Print flux
void RectilinearDecomposition::print(std::ostream& flux) {

  flux<< "  " << std::setw(22) << std::left << "Decomposition" << " : " << "RECTILINEAR" << std::endl
      << "  " << std::setw(22) << std::left << "Domains"       << " : " << "[" << m_numberOfDomainsPerDim << "] " << product(m_numberOfDomainsPerDim)  << std::endl;

}


/// @brief Constructor
/// @param [in] configuration Domain configuration
/// @param [in] rank Domain rank
AnyDecomposition::AnyDecomposition(Configuration<DomainInterface>& configuration, int rank)
  : Decomposition(product(configuration.numberOfDomainsPerDim)) {

  m_rank = rank;

  // Split all cells between domains

  const vec3<int>& globalNumCells = Global::domainInfo.getNumberOfCellsPerDim();

  std::vector<int> start(numberOfDomains);
  std::vector<int> size (numberOfDomains);

  for (uint i=0; i<numberOfDomains; ++i)
    split<int>(product(globalNumCells), numberOfDomains, i, start[i], size[i]);

  // Get my indexes

  m_myGlobalIndexes.resize(size[m_rank]);
  m_neighborOwners.resize(size[m_rank]);

  for (int i=0; i<size[m_rank]; ++i)
    m_myGlobalIndexes[i] = start[m_rank]+i;

  // Get my neighbors

  LCurve<int> crv(globalNumCells);
  
  std::vector<uint> neighborOwners;

  for (uint i=0; i<m_myGlobalIndexes.size(); ++i) {

    neighborOwners.assign(Neighbor::num_neighbors, Neighbor::init);

    const auto myIndex  = m_myGlobalIndexes[i];
    const auto myCoords = crv.convert(myIndex);

    for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr) {

      auto nbrCoords = myCoords + Neighbor::getCoords(nbr);

      // if nbr coords are outside global domain, we must check
      // boundary conditions
      int flag = 0;
      for (int dim=0; dim<VEC3_NDIMS; ++dim) {
	if (nbrCoords[dim]<0 || nbrCoords[dim]>=globalNumCells[dim]) {
	  switch (configuration.boundaryConditions[dim]) {
	  case Configuration<DomainInterface>::FREE:
	  case Configuration<DomainInterface>::WALL:
	  case Configuration<DomainInterface>::FREE_WALL:
	  case Configuration<DomainInterface>::WALL_FREE:
	    flag++;
	    break;
	  case Configuration<DomainInterface>::PERIODIC:
	    nbrCoords[dim] = auxMod(nbrCoords[dim], globalNumCells[dim]);
	    break;
	  }
	}
      }

      // if not periodic and flag is set, no neighbor
      int nbrIndex = (flag==0) ? crv.convert(nbrCoords) : Neighbor::null;

      // now find owner of nbrIndex
      if (nbrIndex == Neighbor::null) {
	neighborOwners[nbr] = Neighbor::null;
      }
      else {
	for (uint j=0; j<numberOfDomains; ++j) {
	  if (nbrIndex<start[j]+size[j]) {
	    neighborOwners[nbr] = j;
	    break;
	  }
	}
      }

    } // END FOR nbr

    m_neighborOwners[i] = neighborOwners;

  } // END FOR i

}


/// @brief Search a cell in the domain
/// @param [in] index Cell index
/// @return True if cell is there
bool AnyDecomposition::isMine(const uint index) { 
  bool test; 
  size_t pos; 
  find(index, m_myGlobalIndexes.data(), m_myGlobalIndexes.size(), test, pos);
  return test;
}


/// @brief Get one neighbor owner of specified cell
/// @param [in] index Index of the cell
/// @param [in] nbrDir Number of the neighbor to get
/// @return Neighbor index
int AnyDecomposition::neighborIndex(const uint index, const uint nbrDir) { 

  bool found; 
  size_t pos; 
  find(index, m_myGlobalIndexes.data(), m_myGlobalIndexes.size(), found, pos);
  
  if (found) {
    return m_neighborOwners[pos][nbrDir];
  }
  else {
    return Neighbor::null;
  }

}


/// @brief Update decomposition part 1 : export and import cells and import neighbors of imported cells
/// @param [in] importIndexes Indexes of imported cells
/// @param [in] exportIndexes Indexes of exported cells
/// @param [in] importEnv Neighbors for imported cells
void AnyDecomposition::update_part_1(const std::vector<uint>& importIndexes, const std::vector<uint>& exportIndexes, const std::vector<uint32_t>& importEnv) {

  bool found; 
  size_t min=0, pos=0; 

  // mark exportations
  for (uint i=0; i<exportIndexes.size(); ++i) {
    find(exportIndexes[i], m_myGlobalIndexes.data()+min, m_myGlobalIndexes.size()-min, found, pos);
    pos+= min;
    min = pos+1;
    m_myGlobalIndexes[pos] = -1;
  }

  // add importations
  m_myGlobalIndexes.insert(m_myGlobalIndexes.end(), importIndexes.begin(), importIndexes.end());

  // sort
  std::sort(m_myGlobalIndexes.begin(), m_myGlobalIndexes.end());

  // get position of first element to remove, then remove
  auto it = find(m_myGlobalIndexes.begin(), m_myGlobalIndexes.end(), -1);
  m_myGlobalIndexes.erase(it, m_myGlobalIndexes.end());

  // Reset and resize neighbors
  m_neighborOwners.clear();
  m_neighborOwners.resize(m_myGlobalIndexes.size());
  
  // Set env for import cells

  parallel_region(0, importIndexes.size(), [&](const uint begin, const uint end) {

      size_t idx=0;
      
      for(uint i=begin; i<end; ++i) {
	
	    find(importIndexes[i], m_myGlobalIndexes.data(), m_myGlobalIndexes.size(), found, idx);
	
	    m_neighborOwners[idx].resize(Neighbor::num_neighbors, Neighbor::init);
	
	    const auto env_start = &importEnv[ i*Neighbor::num_neighbors ];
	    for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr) {
	      m_neighborOwners[idx][nbr] = env_start[nbr];
	    }
	
      }
      
    });

}


/// @brief Update decomposition part 2 : update neighbors for specified cell
/// @param [in] index Cell index
/// @param [in] neighbors Neighboring cells
/// @param [in] allCellOwners Owners for all cells
void AnyDecomposition::update_part_2(uint index, const Array<int>& neighbors, const Array<int>& allCellOwners) {

  bool found; 
  size_t pos; 
  find(index, m_myGlobalIndexes.data(), m_myGlobalIndexes.size(), found, pos);

  if (!found) 
    std::cerr<< "[fatal error] " << __FILE__ << ":" << __LINE__ << ": in AnyDecomposition::update_part_2" << std::endl;

  m_neighborOwners[pos].resize(Neighbor::num_neighbors, Neighbor::init);
  for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr) {
    if (neighbors[nbr]!=Neighbor::null)
    {
      m_neighborOwners[pos][nbr] = allCellOwners[ neighbors[nbr] ]; /// RP

    }
    else
      m_neighborOwners[pos][nbr] = Neighbor::null;
  }

}


/// @brief Print the type of decomposition and number of domains in specified flux
/// @param [in,out] flux Print flux
void AnyDecomposition::print(std::ostream& flux) {

  flux<< "  " << std::setw(22) << std::left << "Decomposition" << " : " << "ANY" << std::endl
      << "  " << std::setw(22) << std::left << "Domains"       << " : " << numberOfDomains << std::endl;

}


/// @brief Create a decomposition of the system in domains from a configuration
/// @param [in] configuration Configuration of the domains
/// @param [in] numberOfNodes Number of nodes (not used)
/// @param [in] rank Node rank
/// @return Created decomposition
Decomposition* Decompose(Configuration<DomainInterface>& configuration, int numberOfNodes, int rank) {

  // Get number of domains in function of 'decoupage', and construct the
  // decomposition

  Decomposition* decomposition = nullptr;

  switch (configuration.decoupage) {

  case (Configuration<DomainInterface>::RECTILINEAR) :
    decomposition = new RectilinearDecomposition(configuration);
    break;

  case (Configuration<DomainInterface>::ANY) :
    decomposition = new AnyDecomposition(configuration, rank);
    break;

  }

  return decomposition;

}
