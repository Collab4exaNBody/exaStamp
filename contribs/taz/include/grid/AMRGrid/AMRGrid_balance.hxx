/// @file
/// @brief Implementation of the balance related methods in class AMRGrid

#if __use_lib_zoltan


/// @brief Get the number of cells with a function adapted to Zoltan load balancing
/// @param [in] data Pointer to the grid
/// @param [out] ierr Error signal
/// @return Number of cells
TMPLSG inline int TMPL_AMRGrid::numberOfCells(void* data, int* ierr) {

	*ierr = ZOLTAN_OK;

	auto grd = static_cast<self_type*>(data);
	auto inf = dynamic_cast<AnyGridInfo*>(&(grd->info));

	return inf->getNumberOfCells();

}


/// @brief Fill global IDs, local IDs and weights of the cells with a function adapted to Zoltan load balancing
/// @param [in] data Pointer to the grid
/// @param [in] num_gid_entries Number of integers used to index a cell globally (not used)
/// @param [in] num_lid_entries Number of integers used to index a cell locally (not used)
/// @param [out] global_ids Global IDs of the cells
/// @param [out] local_ids Local IDs of the cells
/// @param [in] wgt_dim Number of weights associated to each cell
/// @param [out] obj_wgts The to weights of each cell (workload and memory)
/// @param [out] ierr Error signal
TMPLSG void TMPL_AMRGrid::listObj(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr) {

	*ierr = ZOLTAN_OK;

	auto grd = static_cast<self_type*>(data);
	auto inf = dynamic_cast<AnyGridInfo*>(&(grd->info));

	const bool fill_weight = wgt_dim>0;

	const auto& cells = grd->getTraversal(TraversalManager::REAL);
	
  //#pragma omp parallel for // Pas franchement utile
	for (size_t i=0; i<cells.size(); ++i) {

		const auto& cellIndex = cells[i];

		global_ids[i] = inf->getGlobalId(grd->coords[cellIndex]);
		local_ids [i] = cellIndex;

		if (fill_weight) {

			obj_wgts[2*i  ] = (float) grd->particles[cellIndex].computeWorkloadAMR();
			obj_wgts[2*i+1] = (float) grd->particles[cellIndex].computeMemory();

		
		}

	}

}


/// @brief Get the number of dimensions with a function adapted to Zoltan load balancing
/// @param [in] data Pointer to the grid (not used)
/// @param [out] ierr Error signal
/// @return Number of dimensions
TMPLSG inline int TMPL_AMRGrid::numberOfDim(void* data, int* ierr) {

	*ierr = ZOLTAN_OK;

	return 3;

}


/// @brief Set the coordinates specified cell with a function adapted to Zoltan load balancing
/// @param [in] data Pointer to the grid
/// @param [in] num_gid_entries Number of integers used to index a cell globally (not used)
/// @param [in] num_lid_entries Number of integers used to index a cell locally (not used)
/// @param [in] global_id Global ID of the cell (not used)
/// @param [in] local_id Local ID of the cell
/// @param [out] geom_vec Coordinates of the cell
/// @param [out] ierr Error signal
TMPLSG inline void TMPL_AMRGrid::fillCoords(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr) {

	*ierr = ZOLTAN_OK;
	auto grd = static_cast<self_type*>(data);

	geom_vec[0] = grd->coords[*local_id].x;
	geom_vec[1] = grd->coords[*local_id].y;
	geom_vec[2] = grd->coords[*local_id].z;

}


/// @brief Get the number of neighbors of a specified cell with a function adapted to Zoltan load balancing
/// @param [in] data Pointer to the grid
/// @param [in] num_gid_entries Number of integers used to index a cell globally (not used)
/// @param [in] num_lid_entries Number of integers used to index a cell locally (not used)
/// @param [in] global_id Global ID of the cell (not used)
/// @param [in] local_id Local ID of the cell
/// @param [out] ierr Error signal
/// @return Number of edges
TMPLSG inline int TMPL_AMRGrid::numberOfEdges(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr) {

	int numberOfEdges = 0;

	*ierr = ZOLTAN_OK;
	auto grd = static_cast<self_type*>(data);

	const auto& neighborArray = grd->neighbors[*local_id];

	for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr)
	if (neighborArray[nbr]>-1) ++numberOfEdges;

	return numberOfEdges;

}


/// @brief Set the edges/neighbors of a specified cell with a function adapted to Zoltan load balancing
/// @param [in] data Pointer to the grid
/// @param [in] num_gid_entries Number of integers used to index a cell globally (not used)
/// @param [in] num_lid_entries Number of integers used to index a cell locally (not used)
/// @param [in] global_id Global ID of the cell (not used)
/// @param [in] local_id Local ID of the cell
/// @param [out] nbor_global_id Global ID of the cell neighbor
/// @param [out] nbor_procs Processor of the cell neighbor
/// @param [in] wgt_dim Number of weights for each edge (not used)
/// @param [out] ewgts Weight for each edge
/// @param [out] ierr Error signal
/// @return Number of edges
TMPLSG void TMPL_AMRGrid::fillEdges(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs, int wgt_dim, float *ewgts, int *ierr) {

	*ierr = ZOLTAN_OK;
	auto grd = static_cast<self_type*>(data);
	auto inf = dynamic_cast<AnyGridInfo*>(&(grd->info));

	const auto& neighborArray = grd->neighbors[*local_id];
	uint nxt = 0;

	for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr) {

		if (neighborArray[nbr]>-1) {

			auto nbrCoords = grd->coords[*local_id] + Neighbor::getCoords(nbr);
			grd->module(nbrCoords);

			nbor_global_id[nxt] = inf->getGlobalId(nbrCoords);
			nbor_procs    [nxt] = inf->inGrid(nbrCoords) ? grd->getDomainIndex() : inf->getNeighborDomain(Neighbor::getCoords(nbr), grd->coords[*local_id]);

			ewgts         [nxt] = 1 + (grd->particles[*local_id].getNumberOfParticles() + grd->particles[ neighborArray[nbr] ].getNumberOfParticles());
    //  ewgts         [nxt] = 0.5*(grd->particles[*local_id].getNumberOfParticles() * grd->particles[ neighborArray[nbr] ].getNumberOfParticles());
			++nxt;

		}

	}

}


#endif // __use_lib_zoltan


/// @brief Update the owners of real and ghost cells after load balancing
/// @param [out] allCellsOwners Cells owners
/// @param [in] numExport Number of cells exported to other domains
/// @param [in] exportLocalGids Local ID of the cells exported to other domains
/// @param [in] exportProcs Domains where the cells are exported
/// @param [in] numImport Number of cells to import from other domains (not used)
/// @param [in] importGlobalGids Global ID of the cells imported from other domains (not used)
TMPLSG void TMPL_AMRGrid::updateOwners(Array<int>& allCellsOwners, int numExport, uint* exportLocalGids, int* exportProcs, int numImport, uint* importGlobalGids) {

	// const data used later
	const auto domainIndex = this->getDomainIndex();
	const auto& edgeCells  = this->getTraversal(TraversalManager::EDGE_GTHICK);
	const auto& realCells  = this->getTraversal(TraversalManager::REAL);

	// init array, fill with my index
	allCellsOwners = Array<int>(info.getTotalNumberOfCells(), Neighbor::init);

	parallel_region(0, realCells.size(), [&](const uint begin, const uint end) {
		for(size_t i=begin; i<end; ++i)
		{
			assert(domainIndex>=0);
			allCellsOwners[ realCells[i] ] = domainIndex;
		}
	});

	// update array with exports 
	parallel_region(0, numExport, [&](const uint begin, const uint end) {
		for(size_t i=begin; i<end; ++i)
		{
			assert(exportProcs[i]>=0);
			allCellsOwners[ exportLocalGids[i] ] = exportProcs[i];
		}
	});

	// struct to communicate edges ownership
	MessageCenter<MPI__GhostCellOwner> ghostCellsOwners(this->domain->getCommManager(), ISession::GHOST_OWNER); 
	ghostCellsOwners.init(ghostUpdate);
	auto messageSend = ghostCellsOwners.getMessageSend();

	// fill struct with edge cell owners
	for (size_t i=0; i<edgeCells.size(); ++i) {

		auto  cellIndex       = edgeCells[i];
		auto& cellDestDomains = this->destDomains[cellIndex];
		auto& cellDestCells   = this->destCells[cellIndex];

		const uint numDest = cellDestDomains.size();

		for (size_t k=0; k<numDest; ++k) {

			auto& destDomain = cellDestDomains[k];
			auto& destCell   = cellDestCells[k];

		      if (destDomain == (int)domainIndex)
		      {
		      	allCellsOwners[ info.index(destCell) ] = allCellsOwners[cellIndex];
      			assert( info.index(destCell) >= 0);
		      }
			else if (destDomain != Neighbor::null)
				messageSend->push(destDomain, std::move(MPI__GhostCellOwner(destCell, allCellsOwners[cellIndex])));

		}

	}

	// send & collect
	ghostCellsOwners.send();
	ghostCellsOwners.collect();

	// process collected data
	auto messageRecv = ghostCellsOwners.getMessageRecv();
	auto messageSize = messageRecv->getSize();
	auto messageData = messageRecv->getData();

	Array<uint> indexes(messageSize, -1);

	parallel_region(0, messageSize, 
		[&](const uint begin, const uint end) 
		{
			info.findCellIndex(messageData + begin, indexes.data() + begin, end-begin, 
			[&] (const MPI__GhostCellOwner& gco) -> vec3<int> { return gco.cell; });
			for(uint i=begin; i<end; ++i)
			{
				allCellsOwners[ indexes[i] ] = messageData[i].base;
				assert(allCellsOwners[ indexes[i] ]>=0);
			}
		}
	);

	// clean struct
	ghostCellsOwners.clean();

}


struct MyPairItem
{
  ssize_t index;
  ssize_t value;
  inline bool operator < (const MyPairItem& rhs) const
  {
    return value < rhs.value;
  }
};

/// @brief Exchange neighbors of displaced cells
/// @param [in] procsToComm Ranks of the domains with which current domain exchanges particles
/// @param [in] allCellsOwners Cells owners
/// @param [in] numImport Number of cells to import from other domains
/// @param [out] importEnv Neighbors of imported cells
/// @param [in] importGlobalGids Global ID of the cells imported from other domains (not used)
/// @param [in] numExport Number of cells exported to other domains
/// @param [in] exportGlobalGids Global ID of the cells exported to other domains
/// @param [in] exportLocalGids Local ID of the cells exported to other domains
/// @param [in] exportProcs Domains where the cells are exported
TMPLSG void TMPL_AMRGrid::exchangeMovingCellEnv(const Array<uint>& procsToComm, const Array<int>& allCellsOwners, int numImport, std::vector<uint32_t>& importEnv, uint* importGlobalGids, int numExport, uint* exportGlobalGids, uint* exportLocalGids, int* exportProcs)
{

	std::vector<uint32_t> exportEnv;
	std::vector<uint32_t> importMap;

	// init array, fill with my index
	importMap.resize(numImport, -2);
	importEnv.resize(numImport*Neighbor::num_neighbors, -2);
	exportEnv.resize(numExport*Neighbor::num_neighbors, -2);

	// parallel fill
	parallel_region(0, numExport, 
		[&](const uint begin, const uint end) 
		{
			for(uint i=begin; i<end; ++i) {

				auto exportEnv_start = &exportEnv[i*Neighbor::num_neighbors];

				for (uint nbr=0; nbr<Neighbor::num_neighbors; ++nbr) 
				{
					int nbrIdx = this->neighbors[ exportLocalGids[i] ][nbr];
					if (nbrIdx>=0)
					{
						assert(allCellsOwners[nbrIdx]>=0);
						exportEnv_start[nbr] = allCellsOwners[nbrIdx];
					}
					else
					{
						exportEnv_start[nbr] = Neighbor::null;
					}
				}

			}

		}
	);

	//
	MessageCenter<MPI__CellEnv> movingCellEnv(this->domain->getCommManager(), ISession::CELL_SYNC); 
	movingCellEnv.init(this->getDomainIndex(), procsToComm);
	MessageSend<MPI__CellEnv>* messageSend = movingCellEnv.getMessageSend();

	for (int i=0; i<numExport; ++i) 
	{
		MPI__CellEnv tmp;
		tmp.gid = exportGlobalGids[i];
		auxMemCpy(&tmp.neighborOwner[0], &exportEnv[i*Neighbor::num_neighbors], Neighbor::num_neighbors);
		messageSend->push(exportProcs[i], std::move(tmp));
	}

	// send 
	movingCellEnv.send();
	
		// collect
	movingCellEnv.collect();

	// init map (start)
	parallel_region(0, numImport, 
		[&](const uint32_t begin, const uint32_t end) 
		{
			for(uint32_t i=begin; i<end; ++i) importMap[i] = i;
		}
	);

	// process colected data
	auto messageRecv = movingCellEnv.getMessageRecv();
	size_t messageSize = messageRecv->getSize(); // std::vector->size()
	MPI__CellEnv* messageData = messageRecv->getData(); // std::vector->data()

  int dom = this->getDomainIndex();

  std::vector<MyPairItem> tableau_de_tri(messageSize);
  for(size_t i=0;i<messageSize;i++)
  {
    tableau_de_tri[i].index = i;
    tableau_de_tri[i].value = messageData[i].gid;
  }  
  std::sort( tableau_de_tri.begin() , tableau_de_tri.end() );
  
  assert( std::is_sorted(tableau_de_tri.begin() , tableau_de_tri.end()) );
  assert( importMap.size() == messageSize );
  
  for(size_t i=0;i<messageSize;i++) { importMap[i] = tableau_de_tri[i].index; }
  

	assert(std::is_sorted(importMap.begin(), importMap.end(), [messageData,messageSize](const uint32_t& a, const uint32_t& b) -> bool
	{
	  assert( a>=0 && a<messageSize );
		assert( b>=0 && b<messageSize );
		return messageData[a].gid<messageData[b].gid;
	}));

	//
	parallel_region(0, messageSize, 
		[&](const uint begin, const uint end) 
		{
			for(uint i=begin; i<end; ++i)
			{
				auxMemCpy(&(importEnv[i*Neighbor::num_neighbors]), &(messageData[importMap[i]].neighborOwner[0]), Neighbor::num_neighbors);
		  }
		}
	);

	//
	movingCellEnv.clean();
}


/// @brief Updates the neighbors of the cells in the decomposition
/// @param [in] allCellsOwners Cells owners
TMPLSG inline void TMPL_AMRGrid::updateDecomposition(const Array<int>& allCellsOwners) {

	auto inf    = dynamic_cast<AnyGridInfo*>(&info);
	auto decomp = inf->getDecomposition();

	const auto domainIndex = this->getDomainIndex();
	// Work on real cells only
	const auto& cells      = this->getTraversal(TraversalManager::REAL);

	// Parallelize the work by distributing the cells between the threads
	parallel_region(0, cells.size(), 
		[&](const uint begin, const uint end) 
		{

			for(uint i=begin; i<end; ++i) 
			{
				const auto cell = cells[i];
				if (allCellsOwners[cell]!=(int)domainIndex) continue;
				decomp->update_part_2(inf->getGlobalId(this->coords[cell]), this->neighbors[cell], allCellsOwners);
			}
		}
	);




}
