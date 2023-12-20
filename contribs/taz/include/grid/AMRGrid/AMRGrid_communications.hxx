/// @file
/// @brief Implementation of the communication related methods in class AMRGrid

/// @brief Get the particles that leave the domain and move those that just change cell
/// @tparam P Class of the particles
TMPLSG void TMPL_AMRGrid::fillWithLeavingParticles() {

  /// @brief  Shortcut for a tuple containing particle and its recipient cell or domain
  typedef std::tuple<int, P> CellAndP;
  // Temporary lists to avoid threads locks
  std::vector< std::list<CellAndP> > stayingDomain;
  std::vector< std::list<CellAndP> > leavingDomain;

  // Work on edge cells only (other will be dealt with later)
  const auto& cells = this->getTraversal(TraversalManager::EDGE_ONE);

  // Parallelize the work by distributing the cells between the threads
  parallel_region(0, cells.size(),
		  [&](const uint begin, const uint end) {

		    // local lists
		    std::list<CellAndP> stayingDomain_;
		    std::list<CellAndP> leavingDomain_;
		    
		    // Fill local staying and leaving particle arrays for each cell
		    for(uint i=begin; i<end; ++i) {
		      fillWithLeaving_cell(cells[i], stayingDomain_, leavingDomain_);
		    }
		    
		    // Gather staying and leaving arrays
		    this->lock(0);
		    stayingDomain.push_back(std::move(stayingDomain_));
		    this->unlock(0);
		    
		    this->lock(1);
		    leavingDomain.push_back(std::move(leavingDomain_));
		    this->unlock(1);
		    
		  }, comm_partitioner);

  // Put staying particles in right cells
  // Parallelize the work by distributing the cells between the threads (same partitioner as above)
  parallel_region(0, (uint) stayingDomain.size(),
		  [&](const uint begin, const uint end) {

		    for(uint i=begin; i!=end; ++i) {

		      auto& theList = stayingDomain[i];

		      for (auto it=theList.begin(); it!=theList.end(); ++it) {

		      	const uint& c = std::get<0>(*it);

		      	this->lock(c);
		      	particles[c].add(std::get<1>(*it));
		      	this->unlock(c);

		      }

		      theList.clear();

		    }

		  }, comm_partitioner);

  // Where all leaving particles are going to be put
  MessageSend<P>* message = particleExchange.getMessageSend();

  const uint n = leavingDomain.size();

  // Fill the particles message with leaving particles
  for (uint i=0; i<n; ++i) {

    auto& theList = leavingDomain[i];

    for (auto it=theList.begin(); it!=theList.end(); ++it) {
      message->push(std::get<0>(*it), std::get<1>(*it));
    }

    // Change number of particles
    this->numberOfParticles -= theList.size();

    //
    theList.clear();

  }

}


/// @brief Fill arrays with the particles leaving the cell
/// @tparam P Class of the particles
/// @param [in] cellIndex Index of the cell
/// @param [out] staying Particles staying in the domain
/// @param [out] leaving Particles leaving the domain
TMPLSG void TMPL_AMRGrid::fillWithLeaving_cell(uint cellIndex, std::list< std::tuple<int, P> >& staying, std::list< std::tuple<int, P> >& leaving) {

  // References the cell, its neighbors and its coordinates
  auto& cell = this->particles[cellIndex];
  auto& nbrs = this->neighbors[cellIndex];
  auto& crds = this->coords[cellIndex];

  // Boundaries of the cell
  const vec3<double> cMinBounds = Global::domainInfo.getMinBounds() + crds*Global::domainInfo.getCellLength();
  const vec3<double> cMaxBounds = cMinBounds + Global::domainInfo.getCellLength();

  P p;

  // Loop on each particle of the cell
  for (uint i=0, size=cell.getNumberOfParticles(); i<size; ++i) {

  	// Get the movement of the particle
    const int8_t m = cell.computeNeighborMov(i, cMinBounds, cMaxBounds);
	
    // If particle is moving
    if (m != Neighbor::null) {
      vec3<int> destCell = crds + Neighbor::getCoords(m);

      // Convert the particles into MPI type and mark to remove
      cell.getExchange(p, i);
      cell.mark(i);

      uint32_t destCellIndex = info.index(destCell);

      // If moving to a real cell put in staying array
      if (destCellIndex<this->ghostLayer.size() && this->ghostLayer[destCellIndex]==0)
      	staying.emplace_back(nbrs[m], p);
      else {

      	// Do a modulo on cell and coords, then test again
      	this->module(destCell);
      	this->module(p.r);
      	destCellIndex = info.index(destCell);

      	if (info.inGrid(destCell))
      		staying.emplace_back(info.index(destCell), p);
      	else // Else put in leaving array
      		leaving.emplace_back(info.getNeighborDomain(Neighbor::getCoords(m), crds), p);
      }

    }

  }

  // Remove marked particles
    // Remove each particle whose index is in the marked list
  for (uint i=cell.marks.size()-1; i!=(uint)-1; --i) {
    cell.remove(cell.marks[i]);
  }
  
  cell.clear_marked();

}


/// @brief Send the particles to the other domains
///
///
TMPLSG inline void TMPL_AMRGrid::sendLeavingParticles() {
  particleExchange.send();
}


/// @brief Collect incoming particles and add them to the grid
/// @tparam P Class of the particles
TMPLSG inline void TMPL_AMRGrid::collectParticles() {

  static const vec3<double> invCellLength = 1./Global::domainInfo.getCellLength();
  const auto& minBounds     = Global::domainInfo.getMinBounds();

  particleExchange.collect();
  
  MessageRecv<P>* message = particleExchange.getMessageRecv();

  // NOTE : do not change numberOfParticles, it is done in addParticles ...
  
  addParticles(message->getData(), message->getSize(), true, [&] (const P& p) -> vec3<int> { return auxFloor<int>(invCellLength*(p.r-minBounds) ); });

  particleExchange.clean();

  #ifndef NDEBUG
  checkParticlesPosition();
  #endif
}

/// @brief check the position of all atoms
TMPLSG void TMPL_AMRGrid::checkParticlesPosition() {

  // Work on edge cells only (other will be dealt with later)
  const auto& cells = this->getTraversal(TraversalManager::EDGE_ONE);

  // Parallelize the work by distributing the cells between the threads
  parallel_region(0, cells.size(),
		  [&](const uint begin, const uint end) {
		    
		    for(uint i=begin; i<end; ++i) {

		      int cellIndex = cells[i];
                      auto& cell = this->particles[cellIndex];
                      auto& crds = this->coords[cellIndex];

                      const vec3<double> cMinBounds = Global::domainInfo.getMinBounds() + crds*Global::domainInfo.getCellLength();
                      const vec3<double> cMaxBounds = cMinBounds + Global::domainInfo.getCellLength();

                      for (uint j=0; j<cell.size; ++j) 
                      {
                        /* Is in the good cell */
                        if(cell.computeNeighborMov(j, cMinBounds, cMaxBounds) != Neighbor::null) 
                          std::cout << int(j) << " " << cMinBounds + Global::domainInfo.getMinBounds() 
                                    << " " << cell.rx[j] <<  " " << cell.ry[j] <<  " " << cell.rz[j] 
                                    << std::endl;

                        assert ( cell.computeNeighborMov(j, cMinBounds, cMaxBounds) == Neighbor::null );

                        for (uint k=0; k<cell.size; ++k)
                          /* doublon */
                          assert(j==k || cell.id[j] != cell.id[k]); 
                      }
                     
                    }
  
		    
		  }, comm_partitioner);
}



/// @brief Reset cells in the ghost
///
///
TMPLSG inline void TMPL_AMRGrid::clearGhost() {

	// Work on ghost cells only
  const auto& cells = this->getTraversal(TraversalManager::GHOST);

  // Parallelize the work by distributing the cells between the threads
  parallel_region(0, cells.size(),
		  [&](const uint begin, const uint end) {

		    for(uint i=begin; i!=end; ++i) 
		      particles[cells[i]].clear();

		  }, comm_partitioner);

}

extern bool updateVerletLists; // moche, a changer

/// @brief Get the particles that will be send the neighbor domains as ghost
TMPLSG inline void TMPL_AMRGrid::fillWithGhosts() {

   // Simplification of writing
  typedef std::vector<std::pair<size_t,infoOctreeGhost>> tmpTypeVectOctreeInfo;

  // Work on edge cells only
  const auto& cells = this->getTraversal(TraversalManager::EDGE_ONE);

  // If the verlet lists are not updated, shift and size of leafcell are not changed.
  if(updateVerletLists)  // extern value
  {
    m_mapForGhostOctree.clear();
    m_idxLeafsInOctree.clear();
    ghostComm.clear();
    size_t nbThread;

    # pragma omp parallel shared(nbThread)
    {
      nbThread=omp_get_num_threads();
    }

    // Temporary lists to avoid threads locks
    tmpTypeVectOctreeInfo leavingGhost[nbThread];

    
    parallel_region(0, cells.size(),
		  [&](const uint begin, const uint end) {

                    for(int i = begin ; i<end ; i++)
                    {
                      tmpTypeVectOctreeInfo & leavingGhostCell_ = leavingGhost[omp_get_thread_num()];
                      fillWithGhostOtherDomain(cells[i], leavingGhostCell_);
                    }

     });
      
    // Fill the particles message with leaving particles
    for (size_t i=0; i<nbThread; ++i) 
    {
      tmpTypeVectOctreeInfo& theVect = leavingGhost[i];
      for (size_t j=0; j< theVect.size(); ++j)
      {	
	assert(std::get<0>(theVect[j])>=0);
        m_mapForGhostOctree[std::get<0>(theVect[j])].push_back ( std::get<1>(theVect[j]) );
      }
    }
  }
  else
  {

    parallel_region(0, cells.size(),
		  [&](const uint begin, const uint end) {

                     for(size_t i =begin ; i<end ; i++)
                       fillWithGhostSelfDomain(cells[i]);
    });
  }
}

// tmp <- destDomain, destCell, rootIndex,
TMPLSG inline void TMPL_AMRGrid::fillRecWithGhostOtherDomain(std::pair<size_t,infoOctreeGhost> &tmp, std::set<size_t>& tab, Cell *cell, std::vector<std::pair<size_t,infoOctreeGhost>>& leaving)
 {
    if(cell->getIsLeaf()==1)
    {
      if(tmp.second.shift+tmp.second.nbElem != cell->shift)
      {
        if(tmp.second.nbElem != 0) leaving.push_back(tmp);
        tmp.second.nbElem    = int(cell->size);
        tmp.second.shift     = int(cell->shift);
      }
      else
        tmp.second.nbElem += int(cell->size);
    }
    else if(cell->getLevel() == levelMax-1 && Global::reference.isMEAM())
    {
      if(tmp.second.shift+tmp.second.nbElem != cell->shift)
      {
        if(tmp.second.nbElem != 0)  leaving.push_back(tmp);
        tmp.second.nbElem    = int(cell->size);
        tmp.second.shift     = int(cell->shift);
      }
      else
        tmp.second.nbElem += int(cell->size);
    }
    else 
      for(size_t Child:tab)
        fillRecWithGhostOtherDomain(tmp, tab, cell->getDaughterCell(Child), leaving);
 }


/// Get the domains/cells where the cell must be send as a ghost
/// and build the corresponding MPI type for each particle of the cell
/// @param [in] cellIndex Index of the cell
/// @param [out] leaving Ghost particles array
TMPLSG void TMPL_AMRGrid::fillWithGhostSelfDomain(uint cellIndex) {


  const auto domainIndex = this->getDomainIndex();
  const auto ghostThickness = info.getGhostThickness();
  
  auto linear = [] (int x, int y, int z) -> int { return x+1 + 3 * (y+1) + 9 * (z+1); };

  // Get the recipient cells and domain for the cell
  auto& cell            = this->particles[cellIndex];

  if( cell.size == 0) return;

  auto& cellDestDomains = this->destDomains[cellIndex];
  auto& cellDestCells   = this->destCells[cellIndex];

  const uint numDest = cellDestDomains.size();

  assert(cellDestDomains.size() == cellDestCells.size());

  std::set<size_t> index;
  std::vector<int>::iterator it;

  // For each recipient cell
  for (uint k=0; k<numDest; ++k) {

    auto& destDomain = cellDestDomains[k];
    auto& destCell   = cellDestCells[k];

    // If the recipient is in the domain
    if (destDomain == (int)domainIndex) {

      auto destCellIndex = info.index(destCell);

      particles[destCellIndex].clear();

      vec3<int> dcell = this->coords[cellIndex] - destCell;
           
      for(int dim=0; dim < 3 ; dim++)
        if(dcell[dim] < 0) dcell[dim]=-1;
        else if(dcell[dim]>0) dcell[dim]=1;

      // Check if dcell_i{i=x,y,z} belongs to {-1,1,0}
      assert(dcell.x*dcell.x <= 1);
      assert(dcell.y*dcell.y <= 1);    
      assert(dcell.z*dcell.z <= 1);   
      
      int getIndexes = linear(dcell.x,dcell.y,dcell.z);
      std::set<size_t>& t= recDirectionInOctree[getIndexes];
      index.insert(t.cbegin(), t.cend());

      if(k != numDest-1)
      {
        if(cellDestCells[k] != cellDestCells[k+1] || cellDestDomains[k] != cellDestDomains[k+1])
        {
          fillRecWithGhostSelfDomain( destCellIndex, index, &cell);
          index.clear();
        }
      }
      else
      {
        fillRecWithGhostSelfDomain( destCellIndex, index, &cell);
        index.clear();
      }  
    } 
  }
}

/// Get the domains/cells where the cell must be send as a ghost
/// and build the corresponding MPI type for each particle of the cell
/// @param [in] cellIndex Index of the cell
/// @param [out] leaving Ghost particles array
TMPLSG void TMPL_AMRGrid::fillWithGhostOtherDomain(uint cellIndex, std::vector< std::pair<size_t,infoOctreeGhost>>& leaving) {


  const auto domainIndex = this->getDomainIndex();
  const auto ghostThickness = info.getGhostThickness();
  
  auto linear = [] (int x, int y, int z) -> int { return x+1 + 3 * (y+1) + 9 * (z+1); };

  // Get the recipient cells and domain for the cell
  auto& cell            = this->particles[cellIndex];

  if( cell.size == 0) return;

  auto& cellDestDomains = this->destDomains[cellIndex];
  auto& cellDestCells   = this->destCells[cellIndex];

  const uint numDest = cellDestDomains.size();

  assert(cellDestDomains.size() == cellDestCells.size());

  std::set<size_t> index;
  std::vector<int>::iterator it;

  std::pair<size_t,infoOctreeGhost> preFill;

  /* some checks */
  assert(cell.atomDoublon());
  assert(cell.keepRightNumberOfAtoms());


  // For each recipient cell
  for (uint k=0; k<numDest; ++k) {

    auto& destDomain = cellDestDomains[k];
    auto& destCell   = cellDestCells[k];

    // If the recipient is in the domain
    if (destDomain == (int)domainIndex) {

      auto destCellIndex = info.index(destCell);

      particles[destCellIndex].clear();

      assert(particles[destCellIndex].getIsGhost());
      assert(particles[destCellIndex].size == 0);

      for(int i = 0; i < cell.size ; i++)
        for(int j = 0; j < particles[destCellIndex].size ; j++)
          assert( cell.id[i]!=particles[destCellIndex].id[j]);

       // In this case, the number of cells per dimension corresponds to the total number of cells for this dimension 
       
       // For rectilinear
       // Total_number of cells per dims == numberOfCells  :if the "decoupage" is set to 1 in this dimension -> dcell=1,0,-1
       // Total_number of cells per dims  > numberOfCells  :else -> dcell=0;
      vec3<int> dcell = this->coords[cellIndex] - destCell;
           
      for(int dim=0; dim < 3 ; dim++)
        if(dcell[dim] < 0) dcell[dim]=-1;
        else if(dcell[dim]>0) dcell[dim]=1;

      // Check if dcell_i{i=x,y,z} belongs to {-1,1,0}
      assert(dcell.x*dcell.x <= 1);
      assert(dcell.y*dcell.y <= 1);    
      assert(dcell.z*dcell.z <= 1);   
      
      int getIndexes = linear(dcell.x,dcell.y,dcell.z);
      std::set<size_t>& t= recDirectionInOctree[getIndexes];
      index.insert(t.cbegin(), t.cend());

      if(k != numDest-1)
      {
        if(cellDestCells[k] != cellDestCells[k+1] || cellDestDomains[k] != cellDestDomains[k+1])
        {
          fillRecWithGhostSelfDomain( destCellIndex, index, &cell);
          index.clear();
        }
      }
      else
      {
        fillRecWithGhostSelfDomain( destCellIndex, index, &cell);
        index.clear();
      }  
    }
    // Else add the particles to the ghost array
    else if (destDomain != Neighbor::null) {

      vec3<int> dcell = this->raphael_test[cellIndex][k];  ; 

      int getIndexes = linear(dcell.x, dcell.y, dcell.z);
      std::set<size_t>& t= recDirectionInOctree[getIndexes];
      index.insert(t.cbegin(), t.cend());

      /* free free free */
      /* assert(destCell == this->coords[cellIndex]); */


      if(k != numDest-1)
      {
        if(cellDestCells[k] != cellDestCells[k+1] || cellDestDomains[k] != cellDestDomains[k+1])
        {
          preFill.first = destDomain;
          preFill.second.shift=0;
          preFill.second.nbElem=0;
          preFill.second.destCell = destCell;
          preFill.second.cellIndex =cellIndex;
 
          fillRecWithGhostOtherDomain( preFill, index, &cell, leaving);

          leaving.push_back(preFill);
          index.clear();
        }
      }
      else
      {
        preFill.first = destDomain;
        preFill.second.shift=0;
        preFill.second.nbElem=0;
        preFill.second.destCell = destCell;
        preFill.second.cellIndex =cellIndex;

        fillRecWithGhostOtherDomain( preFill, index, &cell, leaving);

        leaving.push_back(preFill);
        index.clear();
      }      
    }
    
  }
  
}

/// @brief Send the ghost particles to the other domains
///
///
TMPLSG inline void TMPL_AMRGrid::sendGhosts() {

  if(updateVerletLists)
  {
    std::set<size_t> domain; 
    Array<int> domainNeigh = info.getNeighborArray();

    for(size_t numDomain=0; numDomain < domainNeigh.size(); numDomain++)
	if(domainNeigh[numDomain] >= 0)
 		domain.insert(domainNeigh[numDomain]);

    octreeToCell(ghostComm.getMapData(), m_mapForGhostOctree, m_idxLeafsInOctree, domain);


    ghostComm.defineCommunication(domain); //V2.2
    ghostComm.adjustCommunication();

    if(Global::reference.isEAM())
    {
	embComm.clear();
        octreeToCell(embComm.getMapData(), m_mapForGhostOctree, m_idxLeafsInOctree, domain);
    	embComm.defineCommunication(domain); 
    	embComm.adjustCommunication();
        assert(embComm.assertRecvInfo(info));
    }

    assert(ghostComm.assertRecvInfo(info));
  }

  packGhost();
  ghostComm.send();
}

/// @brief Collect incoming ghost particles and add them to the grid
TMPLSG inline void TMPL_AMRGrid::collectGhost() {

  ghostComm.recv(); 
  ghostComm.WaitAll();

  unPackGhost();

  // Work on inside cells only (edge cells are already done)
  auto& cells = this->getTraversal(TraversalManager::GHOST);
  
  parallel_region(0, cells.size() , [&](const size_t begin, const size_t end) {

    for(int i = begin ; i < end  ; i++)
    {
      particles[cells[i]].ghostAdjustAMR(this->correcter);
      assert(particles[cells[i]].shift == 0);
    }

  });
}




TMPLSG void TMPL_AMRGrid::fillRecWithGhostSelfDomain(uint destCellIndex, std::set<size_t>& tab, Cell *cell)
 {

    if(cell->getIsLeaf()==1) 
    {
      // If the cell is in the last layer of ghost
      if (this->ghostLayer[destCellIndex]==info.getGhostThickness())
       cell-> template ghostCopyAMR<false, false>(this->particles[destCellIndex],  this->correcter); // Make a local copy
      else // Else, this is an EAM case
        cell-> template ghostCopyAMR<false, true>(this->particles[destCellIndex],  this->correcter); // Make a local copy with a force storage
    }
    else if(cell->getLevel() == levelMax-1 && Global::reference.isMEAM()) /* intercells do not contain the right number of atoms after the internal reorganisation */ 
    {
       // If the cell is in the last layer of ghost
       for(uint8_t Child =0; Child<8 ; ++Child) 
         cell->getDaughterCell(Child)-> template ghostCopyAMR<false, false>(this->particles[destCellIndex],  this->correcter); // Make a local copy
    }
    else 
      for(size_t Child:tab)
        fillRecWithGhostSelfDomain(destCellIndex, tab, cell->getDaughterCell(Child));
 }


/// @brief copy information in the octree storage into the buffer  
/// 
TMPLSG inline void TMPL_AMRGrid::packGhost()
{
  size_t numberOfMsgRecv = ghostComm.getNbMsgSend();
  mapInfo & sendInfo     = ghostComm.getMapData();
  

  m_nbOfParticlesSend +=ghostComm.getNumberOfParticlesSend();


  for(size_t msgId = 0;  msgId<numberOfMsgRecv ; msgId++)
  {
    size_t idDestNode = ghostComm.getDestNode(msgId); // get the node number corresponding to the message number

    parallel_region(0, m_idxLeafsInOctree[msgId].size() , 
                    [&](const size_t begin, const size_t end) {

      for(size_t j = begin; j< end ; j++)
      {
        size_t shift = sendInfo[idDestNode][j].shift;

        for(int it = m_idxLeafsInOctree[msgId][j].first ; it < m_idxLeafsInOctree[msgId][j].second ; it++)
        {
          size_t octreeIndex      = m_mapForGhostOctree[idDestNode][it]. cellIndex;
          Octree & octree        = particles[octreeIndex];
          size_t octreeShiftBegin = m_mapForGhostOctree[idDestNode][it]. shift;
          size_t octreeShiftEnd   = m_mapForGhostOctree[idDestNode][it]. nbElem
                                   +m_mapForGhostOctree[idDestNode][it]. shift;
                                   
          assert(octreeShiftEnd>=octreeShiftBegin);                         
          if(m_mapForGhostOctree[idDestNode][it]. nbElem > 0)        
            ghostComm.pack(msgId, shift, octreeShiftBegin, octreeShiftEnd, octree.rx.data(), octree.ry.data(), octree.rz.data(), octree.id.data(), octree.ti.data()); 
          shift += m_mapForGhostOctree[idDestNode][it]. nbElem;

        }
      }
    });
  }
}


TMPLSG inline uint64_t TMPL_AMRGrid::getTotalNbOfParticlesSend() {return m_nbOfParticlesSend;}
TMPLSG inline uint64_t TMPL_AMRGrid::getTotalNbOfParticlesRecv() {return m_nbOfParticlesRecv;}


/// @brief copy information in the buffer into the octree storage
/// 
TMPLSG inline void TMPL_AMRGrid::unPackGhost()
{

  size_t numberOfMsgRecv = ghostComm.getNbMsgRecv();
  recvInfo & recvInfo = ghostComm.getRecvInfo();

  m_nbOfParticlesRecv +=ghostComm.getNumberOfParticlesRecv();

  for(size_t msgId = 0; msgId < numberOfMsgRecv ; msgId++)
  {

    const size_t beginIndexBuffer = ghostComm.getBeginIndexBuffer(msgId);
    const size_t endIndexBuffer = ghostComm.getEndIndexBuffer(msgId);

    assert(endIndexBuffer >= beginIndexBuffer);

    parallel_region(beginIndexBuffer, endIndexBuffer, 
                    [&](const size_t begin, const size_t end) 
      {

        for(size_t j = begin ; j < end; j++)
        {

	        auto & octree  = particles[info. index(recvInfo[j]. destCell)];
          size_t octreeShift  = (size_t) recvInfo[j]. shift;
          size_t nbElemOctree = (size_t) recvInfo[j]. nbElem;

	        assert(info. index(recvInfo[j]. destCell) >=0);
	        assert(recvInfo[j]. nbElem >=0);
	        assert(recvInfo[j]. shift >=0);
	        assert(nbElemOctree < 10000000); // too large

          octree.resize(nbElemOctree); // Resize ghost octree here

          // buffer [xxxxyyyyzzzzidididtititi]msg0[xxxyyyzzzidididtititi]msg1 ....
          if(nbElemOctree>0)
            ghostComm.unpack( msgId 
                           , octreeShift
                           , nbElemOctree //elem per octree
                           , octree.rx.data() // arrays to fill
                           , octree.ry.data()
                           , octree.rz.data()
                           , octree.id.data()
                           , octree.ti.data());
    }

    });
  }
}  


/// @brief copy information in the octree storage into the buffer  
/// 
TMPLSG inline void TMPL_AMRGrid::packEmb()
{
  size_t numberOfMsgRecv = embComm.getNbMsgSend();
  mapInfo & sendInfo     = embComm.getMapData();
  

  for(size_t msgId = 0;  msgId<numberOfMsgRecv ; msgId++)
  {
    size_t idDestNode = embComm.getDestNode(msgId); // get the node number corresponding to the message number

    parallel_region(0, m_idxLeafsInOctree[msgId].size() , 
                    [&](const size_t begin, const size_t end) {

      for(size_t j = begin; j< end ; j++)
      {
        size_t shift = sendInfo[idDestNode][j].shift;

        for(int it = m_idxLeafsInOctree[msgId][j].first ; it < m_idxLeafsInOctree[msgId][j].second ; it++)
        {
          size_t octreeIndex      = m_mapForGhostOctree[idDestNode][it]. cellIndex;
          Octree & octree        = particles[octreeIndex];
          size_t octreeShiftBegin = m_mapForGhostOctree[idDestNode][it]. shift;
          size_t octreeShiftEnd   = m_mapForGhostOctree[idDestNode][it]. nbElem
                                   +m_mapForGhostOctree[idDestNode][it]. shift;
                                   
          assert(octreeShiftEnd>=octreeShiftBegin);                         
          if(m_mapForGhostOctree[idDestNode][it]. nbElem > 0)        
            embComm.pack(msgId, shift, octreeShiftBegin, octreeShiftEnd, octree.m_eamStorage.emb()); 
          shift += m_mapForGhostOctree[idDestNode][it]. nbElem;

        }
      }
    });
  }
}


/// @brief copy information in the buffer into the octree storage
/// 
TMPLSG inline void TMPL_AMRGrid::unPackEmb()
{

  size_t numberOfMsgRecv = embComm.getNbMsgRecv();
  recvInfo & recvInfo = embComm.getRecvInfo();

  for(size_t msgId = 0; msgId < numberOfMsgRecv ; msgId++)
  {

    const size_t beginIndexBuffer = embComm.getBeginIndexBuffer(msgId);
    const size_t endIndexBuffer = embComm.getEndIndexBuffer(msgId);

    assert(endIndexBuffer >= beginIndexBuffer);

    parallel_region(beginIndexBuffer, endIndexBuffer, 
                    [&](const size_t begin, const size_t end) 
      {

        for(size_t j = begin ; j < end; j++)
        {

	        auto & octree  = particles[info. index(recvInfo[j]. destCell)];
          	size_t octreeShift  = (size_t) recvInfo[j]. shift;
          	size_t nbElemOctree = (size_t) recvInfo[j]. nbElem;

	        assert(info. index(recvInfo[j]. destCell) >=0);
	        assert(recvInfo[j]. nbElem >=0);
	        assert(recvInfo[j]. shift >=0);
	        assert(nbElemOctree < 10000000); // too large

          	octree.m_eamStorage.check(nbElemOctree); // Resize emb octree storage

          	// buffer [xxxxyyyyzzzzidididtititi]msg0[xxxyyyzzzidididtititi]msg1 ....
          	if(nbElemOctree>0)
            	embComm.unpack( msgId 
                           , octreeShift
                           , nbElemOctree //elem per octree
                           , octree.m_eamStorage.emb());
    }

    });
  }
}  



TMPLSG void TMPL_AMRGrid::fillRecWithEmbSelfDomain(size_t destCellIndex, std::set<size_t>& tab, Cell *cell)
 {
    if(cell->getIsLeaf()==1) 
    {
        cell->embCopyAMR(this->particles[destCellIndex]); // Make a local copy with a force storage
    }
    else 
      for(size_t Child:tab)
        fillRecWithEmbSelfDomain(destCellIndex, tab, cell->getDaughterCell(Child));
 }

/// Get the domains/cells where the cell must be send as a ghost
/// and build the corresponding MPI type for each particle of the cell
/// @param [in] cellIndex Index of the cell
/// @param [out] leaving Ghost particles array
TMPLSG void TMPL_AMRGrid::fillWithEmbSelfDomain(size_t cellIndex) {


  const auto domainIndex = this->getDomainIndex();
  const auto ghostThickness = info.getGhostThickness();
  
  auto linear = [] (int x, int y, int z) -> int { return x+1 + 3 * (y+1) + 9 * (z+1); };

  // Get the recipient cells and domain for the cell
  auto& cell            = this->particles[cellIndex];

  if( cell.size == 0) return;

  auto& cellDestDomains = this->destDomains[cellIndex];
  auto& cellDestCells   = this->destCells[cellIndex];

  const uint numDest = cellDestDomains.size();

  assert(cellDestDomains.size() == cellDestCells.size());

  std::set<size_t> index;
  std::vector<int>::iterator it;

  // For each recipient cell
  for (uint k=0; k<numDest; ++k) {

    auto& destDomain = cellDestDomains[k];
    auto& destCell   = cellDestCells[k];

    // If the recipient is in the domain
    if (destDomain == (int)domainIndex) {

      auto destCellIndex = info.index(destCell);

      particles[destCellIndex].resetEAMData();
      //particles[destCellIndex].checkEAMData();

      vec3<int> dcell = this->coords[cellIndex] - destCell;
           
      for(int dim=0; dim < 3 ; dim++)
        if(dcell[dim] < 0) dcell[dim]=-1;
        else if(dcell[dim]>0) dcell[dim]=1;

      // Check if dcell_i{i=x,y,z} belongs to {-1,1,0}
      assert(dcell.x*dcell.x <= 1);
      assert(dcell.y*dcell.y <= 1);    
      assert(dcell.z*dcell.z <= 1);   
      
      int getIndexes = linear(dcell.x,dcell.y,dcell.z);
      std::set<size_t>& t= recDirectionInOctree[getIndexes];
      index.insert(t.cbegin(), t.cend());

      if(k != numDest-1)
      {
        if(cellDestCells[k] != cellDestCells[k+1] || cellDestDomains[k] != cellDestDomains[k+1])
        {
          fillRecWithEmbSelfDomain( destCellIndex, index, &cell);
          index.clear();
        }
      }
      else
      {
        fillRecWithEmbSelfDomain( destCellIndex, index, &cell);
        index.clear();
      }  
    } 
  }
}


TMPLSG inline void TMPL_AMRGrid::fillEmbSelf()
{
	// Work on edge cells only
	const auto& cells = this->getTraversal(TraversalManager::EDGE_ONE);
	#pragma omp parallel for schedule(runtime)
	for(ssize_t i=0; i<cells.size() ; i++)
	{
		fillWithEmbSelfDomain(cells[i]);
	}
}

/// @brief Send embedding terms of particles to other domains
///
///
TMPLSG inline void TMPL_AMRGrid::sendEmb() {
  packEmb();
  embComm.send();
}


/// @brief Collect incoming embedding terms
TMPLSG inline void TMPL_AMRGrid::collectEmb() {

  embComm.recv();
  embComm.WaitAll();
  unPackEmb();
}

/// @brief Reorganize moving particles between the cells of the domain
/// @tparam P Class of the particles
TMPLSG void TMPL_AMRGrid::internalReorganization() {

  const auto& gMinBounds = Global::domainInfo.getMinBounds();
  const auto& cellLength = Global::domainInfo.getCellLength();

  typedef std::tuple<int, P> CellAndP;

  std::vector< std::list<CellAndP> > leavingP;

  // Work on inside cells only (edge cells are already done)
  const auto& cells = this->getTraversal(TraversalManager::INSIDE_ONE);
  
  // Parallelize the work by distributing the cells between the threads
  parallel_region(0, cells.size(),
		  [&](const uint begin, const uint end) {

  			// Fill a local array with leaving particles
		    std::list<CellAndP> li;
		    for(uint i=begin; i!=end; ++i) {
		      const vec3<double> cMinBounds = gMinBounds + this->coords[cells[i]]*cellLength;
		      const vec3<double> cMaxBounds = cMinBounds + cellLength;
                      if(particles[cells[i]].getNumberOfParticles() != 0)
		        particles[cells[i]].internalReorganization(li, this->neighbors[cells[i]], cMinBounds, cMaxBounds);
		    }

		    // Gather local arrays
		    this->lock(0);
		    leavingP.push_back(std::move(li));
		    this->unlock(0);

		  }, comm_partitioner);


  // Put moving particles in right cells
  // Parallelize the work by distributing the cells between the threads (same partitioner as above)
  parallel_region(0, (uint) leavingP.size(),
		  [&](const uint begin, const uint end) {

		    for(uint i=begin; i!=end; ++i) {

		      auto& theList = leavingP[i];

		      for (auto it=theList.begin(); it!=theList.end(); ++it) {

		      	const uint& c = std::get<0>(*it);

		      	this->lock(c);
		      	particles[c].add(std::get<1>(*it));
		      	this->unlock(c);

		      }

		      theList.clear();

		    }

		  }, comm_partitioner);

}





/// @brief After a load balancing, clear the grid by storing or sending all the particles
/// @tparam P Class of the particles
/// @param [out] toKeep Storage for the particles to keep
/// @param [in,out] toSend Message with the particles to send
/// @param [in] sendIndexes Indexes of the particles to send
/// @param [in] sendProcs Where to send the particles
/// @param [in] sendSize Number of particles to sendgetExchange
TMPLSG void TMPL_AMRGrid::dumpParticles(std::vector<P>& toKeep, MessageSend<P>& toSend, const uint* sendIndexes, const int* sendProcs, const uint sendSize) {

	// Evaluate number of particles to send and to keep
  uint numberOfParticlesToSend = ( 1+this->numberOfParticles/info.getNumberOfCells() ) * sendSize;
  toKeep.reserve(std::max<int>(1, (int) (this->numberOfParticles-numberOfParticlesToSend)));

  P p;

  // For each cell to send
  for (uint i=0; i<sendSize; ++i) {

    auto& cell = particles[sendIndexes[i]];

    // For each particle of the cell
    for (uint j=0; j<cell.getNumberOfParticles(); ++j) {
    	// Get the exchange data for this particle
      cell.getExchange(p, j);
      // Add it to the message
      toSend.push(sendProcs[i], p);
    }

    // Clear the cell
    cell.clear();

  }

  // Work on real particles only
  const auto& cells = this->getTraversal(TraversalManager::REAL);

  // For each cell
  for (uint i=0; i<cells.size(); ++i) {

  	// Store all the particles
    particles[ cells[i] ].dump (toKeep);
    // Clear the cell
    particles[ cells[i] ].clear();
  }

  assert(toKeep.size() + toSend.getTotalSize() == this->numberOfParticles );
  this->numberOfParticles = 0;

}


/// @brief copy information in the octree storage into the buffer  
/// 
TMPLSG inline void TMPL_AMRGrid::packBalance(exchangeGhost &balanceComm)
{
  size_t numberOfMsgRecv = balanceComm.getNbMsgSend();
  mapInfo & sendInfo     = balanceComm.getMapData();
  

  for(size_t msgId = 0;  msgId<numberOfMsgRecv ; msgId++)
  {
    size_t idDestNode = balanceComm.getDestNode(msgId); // get the node number corresponding to the message number
    
    #pragma omp parallel for schedule(runtime)
    for(size_t j = 0; j< sendInfo[idDestNode].size() ; j++)
    {
      size_t shift = sendInfo[idDestNode][j].shift;

      size_t octreeIndex      = info.index(sendInfo[idDestNode][j].destCell); 
      Octree & octree        = particles[octreeIndex];
      size_t octreeShiftBegin = 0;
      assert(sendInfo[idDestNode][j].nbElem == particles[octreeIndex].size);
      size_t octreeShiftEnd   = sendInfo[idDestNode][j].nbElem;
                                 
      assert(octreeShiftEnd>=octreeShiftBegin);                         
      if(sendInfo[idDestNode][j]. nbElem > 0)        
        balanceComm.pack(msgId, shift, octreeShiftBegin, octreeShiftEnd, 
          octree.rx.data(), 
          octree.ry.data(), 
          octree.rz.data(),
          octree.vx.data(), // arrays to fill
          octree.vy.data(),
          octree.vz.data(),
 /*         octree.fx.data(), // arrays to fill
          octree.fy.data(),
          octree.fz.data(),
          octree.ep.data(),  */        
          octree.id.data(), 
          octree.ti.data()
        ); 
          
      octree.clear();

    }

  }
}

TMPLSG inline void TMPL_AMRGrid::packBalanceKeep(keepAtom &toKeep)
{
  std::vector<infoCellGhost>& infoKeep = toKeep.getInfo();
  
  toKeep.allocateBuffer();
  
  #pragma omp parallel for schedule(runtime)
  for(int i = 0 ; i < infoKeep.size() ; i++) 
  {
    size_t shift = infoKeep[i].shift;
    size_t nbElem = infoKeep[i].nbElem;
    int octreeIndex = info.index(infoKeep[i].destCell);
    Octree & octree = particles[octreeIndex];
    
    toKeep.fillBuffer(shift,
      nbElem,
      octree.rx.data(), // arrays to fill
      octree.ry.data(),
      octree.rz.data(),
      octree.vx.data(), // arrays to fill
      octree.vy.data(),
      octree.vz.data(),
/*      octree.fx.data(), // arrays to fill
      octree.fy.data(),
      octree.fz.data(),
      octree.ep.data(),      */
      octree.id.data(),
      octree.ti.data());  
      
    octree.clear();      
    
  }
}


TMPLSG inline void TMPL_AMRGrid::unPackBalanceKeep(keepAtom &toKeep)
{
  std::vector<infoCellGhost>& infoKeep = toKeep.getInfo();
  
  #pragma omp parallel for schedule(runtime)
  for(int i = 0 ; i < infoKeep.size() ; i++) 
  {
    size_t shift = infoKeep[i].shift;
    size_t nbElem = infoKeep[i].nbElem;
    int octreeIndex = info.index(infoKeep[i].destCell);
    Octree & octree = particles[octreeIndex];
    
    octree.resize(infoKeep[i].nbElem);
    
    toKeep.decodeBuffer(shift,
      nbElem,
      octree.rx.data(), // arrays to fill
      octree.ry.data(),
      octree.rz.data(),
      octree.vx.data(), // arrays to fill
      octree.vy.data(),
      octree.vz.data(),
/*      octree.fx.data(), // arrays to fill
      octree.fy.data(),
      octree.fz.data(),
      octree.ep.data(),   */         
      octree.id.data(),
      octree.ti.data());      
    
  }

  this->numberOfParticles += toKeep.m_nbElem;
}

TMPLSG inline void TMPL_AMRGrid::unPackBalance(exchangeGhost &balanceComm)
{

  size_t numberOfMsgRecv = balanceComm.getNbMsgRecv();
  recvInfo & recvInfo = balanceComm.getRecvInfo();

  assert(balanceComm.assertRecvInfo(info));

  for(size_t msgId = 0; msgId < numberOfMsgRecv ; msgId++)
  {

    const size_t beginIndexBuffer = balanceComm.getBeginIndexBuffer(msgId);
    const size_t endIndexBuffer = balanceComm.getEndIndexBuffer(msgId);

    assert(endIndexBuffer >= beginIndexBuffer);

    #pragma omp parallel for schedule(runtime)
    for(size_t j = beginIndexBuffer; j < endIndexBuffer; j++)
    {
      auto & octree  = particles[info. index(recvInfo[j]. destCell)];
      size_t octreeShift  = (size_t) recvInfo[j]. shift;
      size_t nbElemOctree = (size_t) recvInfo[j]. nbElem;
      
      assert(info. index(recvInfo[j]. destCell) >=0);
      assert(recvInfo[j]. nbElem >=0);
      assert(recvInfo[j]. shift >=0);
      assert(nbElemOctree < 10000000); // too large

      octree.resize(nbElemOctree); // Resize ghost octree here

      // buffer [xxxxyyyyzzzzidididtititi]msg0[xxxyyyzzzidididtititi]msg1 ....
      if(nbElemOctree>0)
        balanceComm.unpack( msgId 
                       , octreeShift
                       , nbElemOctree //elem per octree
                       , octree.rx.data() // arrays to fill
                       , octree.ry.data()
                       , octree.rz.data()
                       , octree.vx.data() // arrays to fill
                       , octree.vy.data()
                       , octree.vz.data()
          /*             , octree.fx.data() // arrays to fill
                       , octree.fy.data()
                       , octree.fz.data() 
                       , octree.ep.data()  */                      
                       , octree.id.data()
                       , octree.ti.data());
    }

  }
   
  this->numberOfParticles += balanceComm.numberOfElementsRecv();
}

/// @brief After a load balancing, clear the grid by storing or sending all the particles
/// @tparam P Class of the particles
/// @param [out] toKeep Storage for the particles to keep
/// @param [in,out] toSend Message with the particles to send
/// @param [in] sendIndexes Indexes of the particles to send
/// @param [in] sendProcs Where to send the particles
/// @param [in] sendSize Number of particles to send
TMPLSG void TMPL_AMRGrid::dumpParticles(keepAtom& toKeep, exchangeGhost &balanceComm, const uint* sendIndexes, const int* sendProcs, const uint sendSize) {

  auto& mapToSend = balanceComm.getMapData();
  infoCellGhost tmp;
  
  tmp.shift = 0;   
  
  // For each cell to send
  for (uint i=0; i<sendSize; ++i) {
    tmp.nbElem = particles[sendIndexes[i]].size;

    tmp.destCell = this->coords[sendIndexes[i]];
    
    if(tmp.nbElem > 0)
      mapToSend[sendProcs[i]].push_back(tmp);
  } 


  std::set<size_t> domain;
  
  // Communication all to all
  for(size_t i=0 ; i< balanceComm.numberOfNodes ; ++i) 
  {
    if(i != balanceComm.rank)
      domain.insert(i);
  }

  // On adapte la valeur shift correspondant à l'endroit ou la cellule / octree va être copier dans le buffer d'envois / réception
  for(auto& it : mapToSend)
  {
    std::vector<infoCellGhost>& vecElem = it.second;
    for(int i = 1; i < vecElem.size() ; i ++)
    {
      vecElem[i].shift = vecElem[i-1].shift + vecElem[i-1].nbElem;
    }
  }

  balanceComm.defineCommunication(domain); //V2.2
  balanceComm.adjustCommunication();
  
  assert(balanceComm.assertSendInfo()); // Vérifie si le nombre à envoyer et adapté avec les informations transmise par les infCellGhosts

  packBalance(balanceComm); // Pack les atomes dans un buffer d'envois
  balanceComm.send(); // Mpi Send
  balanceComm.recv(); // Mpi Recieve
  balanceComm.WaitAll(); // Attente de la fin des communications
  

  // Work on real particles only
  const auto& cells = this->getTraversal(TraversalManager::REAL);

  auto& infoCells = toKeep.getInfo(); 

  tmp.shift = 0;

  // For each cell
  for (size_t i=0; i<cells.size(); ++i) {

    if(particles[cells[i]].size>0)
    {
      tmp.nbElem = particles[cells[i]].size;
      tmp.destCell = this->coords[cells[i]];
      infoCells.push_back(tmp);
    }
  }
  
  toKeep.computeShift();
  packBalanceKeep(toKeep); 
    
  this->numberOfParticles = 0;

}




  

