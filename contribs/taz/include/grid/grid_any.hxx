/// @file
/// @brief Implementation of the methods used to build the grid in the case of an anygrid

/// @brief Build the grid (case rectilinear grid)
/// @tparam Grid_impl The proper Grid subclass
/// @param [in,out] domain_ Domain
/// @param [in] info Grid info
TMPLG void TMPL_Grid::build(Domain<Grid_impl>* domain_, AnyGridInfo* info) {

  const auto edgeTest = info->getGlobalEdges();
  const auto modTest  = info->getCellSup()-info->getCellInf() >= Global::domainInfo.getNumberOfCellsPerDim()/2;

  // Set the correcter and modulers
  correcter.set(edgeTest || modTest, Global::domainInfo.getExtension());
  modD.set     (edgeTest           , Global::domainInfo.getMinBounds(), Global::domainInfo.getMaxBounds());
  modI.set     (edgeTest           , zeros()                          , Global::domainInfo.getNumberOfCellsPerDim());

  // Link to the domain
  domain = domain_;

  const int  ghostThickness     = info->getGhostThickness();
  const uint totalNumberOfCells = info->getTotalNumberOfCells();

  // Reset mutexes
  mtx.check(totalNumberOfCells);

  // Reset size of the vectors
  coords     .resize(totalNumberOfCells);
  neighbors  .resize(totalNumberOfCells, Array<int>(Neighbor::num_neighbors, Neighbor::init));
  ghostLayer .resize(totalNumberOfCells, -1);
  destDomains.resize(totalNumberOfCells);
  destCells  .resize(totalNumberOfCells);
  raphael_test.resize(totalNumberOfCells); // warning

  // Get cells
  for (uint cell=0; cell<totalNumberOfCells; ++cell) {
    coords[cell] = info->coords(cell);
  }

  // Create a traversal Database
  makeTraversals(info, coords);

  // Set ghost layers
  const Array<uint>& ghostCells = getTraversal(TraversalManager::GHOST);
  for (uint i=0; i<ghostCells.size(); ++i) { ghostLayer[ghostCells[i]] = (uint8_t) ghostThickness;
  }

  const Array<uint>& realCells= getTraversal(TraversalManager::REAL);
  for (uint i=0; i<realCells.size(); ++i)  {ghostLayer[realCells[i]] = (uint8_t) 0;}

  // grep rp j'ai rajouté la double boucle suivante car j'ai besoin de calculer dans les ghosts

  // Set neighbor data
    if(Global::reference.isMEAM())//MEAM version -> we compute on reals and ghosts cells
    {
        const Array<uint>& all = getTraversal(TraversalManager::ALL);
        for (uint i=0; i<all.size(); ++ i) 
        {
            const auto& cellG = all[i];
            for (uint nb=0; nb<Neighbor::num_neighbors; ++nb) 
            {
                // compute coordinates of neighboring cell
                vec3<int> neighborCoords1 = coords[cellG] + Neighbor::getCoords(nb);
                neighbors[cellG][nb] = info->index(neighborCoords1);
            }
        }
    }
    else
    {
      // Set neighbor data
      const Array<uint>& all_but_one = getTraversal(TraversalManager::ALL_BUT_ONE);

      for (uint i=0; i<all_but_one.size(); ++ i) {
    
        const auto& cell = all_but_one[i];
        
        for (uint nb=0; nb<Neighbor::num_neighbors; ++nb) {
          // compute coordinates of neighboring cell
          vec3<int> neighborCoords = coords[cell] + Neighbor::getCoords(nb);
          neighbors[cell][nb] = info->index(neighborCoords);
        }   
      }
    }

  coords.shrink_to_fit();
  neighbors.shrink_to_fit();

  // Set neighbors of the edge cells
  info->setNeighbors(getTraversal(TraversalManager::EDGE_ONE), neighbors, modI);

}


#include <tuple>
#include <assert.h>

/// @brief Fill the list of recipient domains and cells for ghost exchange
/// @param [in] info Grid info
TMPLG void TMPL_Grid::setGhostData(AnyGridInfo* info) {

  // WARNING : will work only with one ghost layer

  const auto& globalNumberOfCellsPerDim = Global::domainInfo.getNumberOfCellsPerDim();

  destCells.clear();

  const auto& cells = getTraversal(TraversalManager::EDGE_ONE);

  // For each cell
  for (uint i=0; i<cells.size(); ++i) {

    const auto& cellIndex = cells[i];
    const auto& cellCoord = this->coords[cellIndex];



    // test if cell is on a global edge
    vec3<int> globalEdge(0);
    for (int dim=0; dim<VEC3_NDIMS; ++dim) {
      if (cellCoord[dim]==0) 
	globalEdge[dim] -= 1;
      else if (cellCoord[dim]==globalNumberOfCellsPerDim[dim]-1)
	globalEdge[dim] += 1;
    }
    // END 

    // raphael
    std::vector< std::tuple<int, vec3<int>, vec3<int>> > tmp1;
    
    // For each neighbor
    for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr) {
      
      auto nbrMov = Neighbor::getCoords(nbr);

      // Get recipient domain
      auto destDomain = info->getNeighborDomain(nbrMov, cellCoord);

      if (destDomain >= 0) {



      	vec3<int> correct(0);
      	// Get correction for periodic boundary case
      	for (int dim=0; dim<VEC3_NDIMS; ++dim)
      		correct[dim] = nbrMov[dim]==globalEdge[dim] ? nbrMov[dim] : 0;

      	// Add recipient domain and cell to tmp1
      	tmp1.push_back(std::make_tuple(destDomain, cellCoord - correct * Global::domainInfo.getNumberOfCellsPerDim(), nbrMov));



      }

    }

    // Sort tmp1
    std::sort(tmp1.begin(), tmp1.end(), [] (const std::tuple<int,vec3<int>, vec3<int>>& a, const std::tuple<int, vec3<int>, vec3<int> >& b) -> bool { 
	return std::get<0>(a)< std::get<0>(b) || ( std::get<0>(a)== std::get<0>(b) && ( std::get<1>(a).x < std::get<1>(b).x || ( std::get<1>(a).x == std::get<1>(b).x && ( std::get<1>(a).y < std::get<1>(b).y || ( std::get<1>(a).y== std::get<1>(b).y && ( std::get<1>(a).z< std::get<1>(b).z ) ) ) ) ) );
      });

    // Remove duplicates
    auto it = std::unique(tmp1.begin(), tmp1.end());
    tmp1.resize(std::distance(tmp1.begin(), it));

    // Fill recipient lists
    for (auto& elem : tmp1) {

      destDomains[cellIndex].push_back(std::get<0>(elem));
      destCells  [cellIndex].push_back(std::get<1>(elem));
      raphael_test[cellIndex].push_back(std::get<2>(elem));
      

    }

    destDomains[cellIndex].shrink_to_fit();
    destCells  [cellIndex].shrink_to_fit();

  }

}
