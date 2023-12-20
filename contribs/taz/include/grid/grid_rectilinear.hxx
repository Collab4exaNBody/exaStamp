/// @file
/// @brief Implementation of the methods used to build the grid in the case of a rectilinear grid

/// @brief Build the grid (case rectilinear grid)
/// @param [in,out] domain_ Domain
/// @param [in] info Grid info
TMPLG void TMPL_Grid::build(Domain<Grid_impl>* domain_, RectilinearGridInfo* info) {

  const auto edgeTest = info->getGlobalEdges();
  const auto modTest  = info->getNumberOfCellsPerDim() >= Global::domainInfo.getNumberOfCellsPerDim()/2;

  // Set the correcter and modulers
  correcter.set(edgeTest || modTest, Global::domainInfo.getExtension());
  modD.set     (edgeTest, Global::domainInfo.getMinBounds(), Global::domainInfo.getMaxBounds());
  modI.set     (edgeTest, zeros(), Global::domainInfo.getNumberOfCellsPerDim());

  // Link to the domain
  domain = domain_;

  const vec3<int> origin              = info->getOrigin();
  const vec3<int> numberOfCellsPerDim = info->getNumberOfCellsPerDim();
  const int ghostThickness            = info->getGhostThickness();
  const int totalNumberOfCells        = info->getTotalNumberOfCells();

  // Reset mutexes
  mtx.check(totalNumberOfCells);

  // Reset size of the vectors
  coords.resize(totalNumberOfCells);
  neighbors.resize(totalNumberOfCells, Array<int>(Neighbor::num_neighbors, Neighbor::init));
  ghostLayer.resize(totalNumberOfCells, -1);
  destDomains.resize(totalNumberOfCells);
  destCells.resize(totalNumberOfCells);
  raphael_test.resize(totalNumberOfCells); // warning
  const vec3<int> offset = origin - ghostThickness;


  // Loop over all cells (ghost included) : initializing cell coords (physical)
  for (int cell=0; cell<totalNumberOfCells; ++cell) {
    coords[cell] = offset + info->convert(cell);
  }


  // Create a traversal Database
  makeTraversals(info, coords);

  // Set ghost layers
  if (ghostThickness > 2) std::cout<< "There will be a problem here (Grid<>::build)" << std::endl;

  const Array<uint>& ghostCells = getTraversal(TraversalManager::GHOST);
  for (uint i=0; i<ghostCells.size(); ++i) ghostLayer[ghostCells[i]] = (uint8_t) ghostThickness;

  if (ghostThickness>1) {
    const Array<uint>& ghostCells2 = getTraversal(TraversalManager::ALL_BUT_ONE);
    for (uint i=0; i<ghostCells2.size(); ++i) ghostLayer[ghostCells2[i]] = (uint8_t) ghostThickness-1;
  }

  const Array<uint>& realCells = getTraversal(TraversalManager::REAL);
  for (uint i=0; i<realCells.size(); ++i) ghostLayer[realCells[i]] = (uint8_t) 0;

  // Loop over all cells except one layer (of ghost) and init neighbor data
  const Array<uint>& all_but_one = getTraversal(TraversalManager::ALL_BUT_ONE);

  for (uint i=0; i<all_but_one.size(); ++ i) {

    const auto& cell = all_but_one[i];
    
    for (uint nb=0; nb<Neighbor::num_neighbors; ++nb) {
      // compute coordinates of neighboring cell
      vec3<int> neighborCoords = coords[cell] - offset + Neighbor::getCoords(nb);
      neighbors[cell][nb] = info->convert(neighborCoords);
    }
  } // end for cell

  if(Global::reference.isMEAM())//MEAM version -> we compute on real and ghost cells
  {
    const Array<uint>& ghost = getTraversal(TraversalManager::GHOST);
    for (uint i=0; i<ghost.size(); ++ i) 
    {
      const auto& cellG = ghost[i];
      for (uint nb=0; nb<Neighbor::num_neighbors; ++nb) 
      {
        // compute coordinates of neighboring cell for ghost cells
        vec3<int> neighborCoords = coords[cellG] - offset + Neighbor::getCoords(nb);
        //if cell are in the domaine (id between 0 and the number of cells
        if(info->convert(neighborCoords)>=0 && info->convert(neighborCoords)< (int)(totalNumberOfCells))  neighbors[cellG][nb] = info->convert(neighborCoords);
      }
    }
  }

  coords.shrink_to_fit();
  neighbors.shrink_to_fit();
}


/// @brief Fill the list of recipient domains and cells for ghost exchange
/// @param [in] info Grid info
TMPLG void TMPL_Grid::setGhostData(RectilinearGridInfo* info) {

  const vec3<int> origin = info->getOrigin();
  const vec3<int> numberOfCellsPerDim = info->getNumberOfCellsPerDim();
  const int ghostThickness = (int) Global::reference.getGhostThickness();

  int domainIndex = getDomainIndex();

  // Given a neighborhood, this function compute a size and an origin giving
  // the set of cells which should be sent to the corresponding neighbors
  auto getLocalData = [&] (vec3<int>& mov, vec3<int>& orgn, vec3<int>& size) -> void {

  	// Origin at a depth of 2 ghostThickness to account for passage through the ghost layer
  	// and a ghostThick layer of cells to send
    orgn = vec3<int>(2*ghostThickness);
    size = numberOfCellsPerDim-2*ghostThickness;

    for (int dim=0; dim<VEC3_NDIMS; ++dim) {
      if      (mov[dim]>0) orgn[dim] = numberOfCellsPerDim[dim]; // + ghostThickness for the ghost to pass through, - ghostThickness for the layer to send
      else if (mov[dim]<0) orgn[dim] = ghostThickness;

      if (mov[dim]!=0) size[dim] = ghostThickness;
    }

  };

  // Loop on neighborhoods of the domain (faces, edges and vertexes)
  for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr) {

    vec3<int> nbrMov, localOrgn, localSize;
    nbrMov = Neighbor::getCoords(nbr);
    // Get the set of cells concerned by this neighborhood (no edge cells for the faces or vertex cells for the edges)
    getLocalData(nbrMov, localOrgn, localSize);

    // For each neighborhoods get the list of moves to access neighbor domains
    // For faces neighborhoods, there is only one neighbor domain, but for edges there three and for vertexes there a seven
    Array< vec3<int> > destMov = __getMovList(nbrMov);
 
    // Loop on concerned cells
    LCurve<int> localCurve(localSize);
    for (int i=0; i<product(localSize); ++i) {

    	// Get the index of the cell
      int c = info->convert(localOrgn+localCurve.convert(i));

      // Loop on neighbor domains
      for (uint m=0; m<destMov.size(); ++m) {

      	// Get index of neighbor domain and coordinates of the cell
      	int destDomain = info->getNeighborDomain(destMov[m], coords[c]);
      	vec3<int> destCell = coords[c];

      	if (destDomain!=Neighbor::null) {

      		// Add domain as recipient during ghost exchange
      		destDomains[c].push_back(destDomain);
                raphael_test[c].push_back(destMov[m]); //warning

      		// Get the recipient cell
      		if (destDomain==domainIndex) {
      			// Case where the recipient domain is the sending domain (through periodic boundary)
      			destCell -= destMov[m]*Global::domainInfo.getNumberOfCellsPerDim();
      		}
      		else {
      			vec3<int> cross(0);
      			for (int dim=0; dim<VEC3_NDIMS; ++dim) {
      				// Cases of periodic boundary
      				if      (destCell[dim]+ghostThickness*destMov[m][dim]<0) cross[dim]=-1;
      				else if (destCell[dim]+ghostThickness*destMov[m][dim]>=(Global::domainInfo.getNumberOfCellsPerDim())[dim]) cross[dim]=1;
      			}
      			destCell-=cross*Global::domainInfo.getNumberOfCellsPerDim();
      		}
      		// Add recipient cell
      		destCells[c].push_back(destCell);

      	}

      }

      destDomains[c].shrink_to_fit();
      destCells[c].shrink_to_fit();

    }

  }

  destDomains.shrink_to_fit();
  destCells.shrink_to_fit();

}
