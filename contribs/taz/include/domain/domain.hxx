/// @file
/// @brief Functions in class Domain

/// @brief Accessor to index
inline uint LocalInfo::getIndex() { 
  return index;
}


/// @brief Get symmetrization and convert to bool
/// @return True if symmetrization value is 1
inline bool LocalInfo::symmetrize() { 
  return symmetrization==1;
}


/// @brief Print symmetrization parameter of the LocalInfo class in specified flux
/// @param [in,out] flux Print flux
inline void LocalInfo::print(std::ostream& flux) {

  std::string str, bol;

  if      (symmetrization==1) bol="yes";
  else if (symmetrization==0) bol="no";
  else bol="__";

  flux << "  " << std::setw(22) << std::left << "Symetrize (def)" << " : " << bol << std::endl;

}




/// @brief Constructor with arguments
/// @tparam Grid_impl The proper Grid subclass
/// @param [in,out] node Node
/// @param [in] index Domain index
/// @param [in] config Configuration of the domain
/// @param [in] numberOfParticles Number of particles on the domain
TMPLG TMPL_Domain::Domain(Node* node, uint index, const Configuration<DomainInterface>& config, uint64_t numberOfParticles)
  : LocalInfo(index, config.symmetrization),
    m_inputOutputManager(node->inputOutput()),
    ref(&node->getManager()),
    m_balanceComm(this->getCommManager()->getCommunicator()), m_toKeep(),
    grid(nullptr) {

  grid = new Grid_impl(this, index, numberOfParticles);
  
  void * tmpPtr = nullptr;
  
  // storage of atoms that remains on this domain
  m_toKeep.defineSizeOfData(
    (double*)   (tmpPtr),	//rx
    (double*)   (tmpPtr),	//ry
    (double*)   (tmpPtr),	//rz
    (double*)   (tmpPtr),	//vx
    (double*)   (tmpPtr),	//vy
    (double*)   (tmpPtr),	//vz
    (uint64_t*) (tmpPtr),	//id
    (uint8_t*)  (tmpPtr)	//ti
  );

  // storage of atoms that leaves this domain
  m_balanceComm.defineSizeOfData(
    (double*)   (tmpPtr),	//rx
    (double*)   (tmpPtr),	//ry
    (double*)   (tmpPtr),	//rz
    (double*)   (tmpPtr),	//vx
    (double*)   (tmpPtr),	//vy
    (double*)   (tmpPtr),	//vz
    (uint64_t*) (tmpPtr),	//id
    (uint8_t*)  (tmpPtr)	//ti
  );
  
  // Describe the maximal number of atom send per mpi message
   m_balanceComm.defineSizeOfDataPerMsg(1000000);
}

/// @brief Destructor
///
///
TMPLG TMPL_Domain::~Domain() {
  if (grid!=nullptr) delete grid;
}


/// @brief Initialize device if running on GPU or something else
///
/// Nothing to do yet
TMPLG inline void TMPL_Domain::initDevice() {
}


/// @brief Initialize a chunk of the particles from a configuration
/// @param [in] particles Configuration of all the particles
/// @param [in] cBegin First particle to initialize
/// @param [in] cEnd End of the particle to initialize
TMPLG inline void TMPL_Domain::initParticles(Configuration<MPI__Particle>& particles, uint64_t cBegin, uint64_t cEnd) {
  grid->initParticles(particles, cBegin, cEnd);
}



TMPLG void TMPL_Domain::aferrgreg(uint64_t& send, uint64_t & recv) {

  send = grid->getTotalNbOfParticlesRecv();
  recv = grid->getTotalNbOfParticlesSend();

}



/// @brief Initialize a chunk of the particles from an array of particles
/// @param [in] particles Particles
/// @param [in] initialT Initial temperature
/// @param [in] numP Number of particles
/// @param [in] nCells Number of duplicates of the base cell in each dimension, default 1
/// @param [in] cellSize Size of the base cell, default 0 (if not used)
/// @param [in] numMol Number of molecules in the base cell, default 0
TMPLG inline void TMPL_Domain::initParticles(MPI__Particle** particles, double initialT, uint64_t numP, vec3<int> nCells, vec3<double> cellSize, uint64_t numMol) {
  grid->initParticles(particles, initialT, numP, nCells, cellSize, numMol);
}


/// @brief Update particles positions using to first order approximation of the movement equations
/// @param [in] time Time step
TMPLG inline void TMPL_Domain::pushPositions1stOrder(double time) {
  grid->pushPositions1stOrder(time);
}


/// @brief Update particles positions using to second order approximation of the movement equations
/// @param [in] time Time step
TMPLG inline void TMPL_Domain::pushPositions2ndOrder(double time) {
  grid->pushPositions2ndOrder(time);
}


/// @brief Update particles velocities using to first order approximation of the movement equations
/// @param [in] time Time step
TMPLG inline void TMPL_Domain::pushVelocities1stOrder(double time) {
  grid->pushVelocities1stOrder(time);
}


/// @brief Update particles velocities by applying Langevin dissipation
/// @param [in] time Time step
/// @param [in] gamma Gamma parameter of the Langevin scheme
/// @param [in] beta Beta parameter of the Langevin scheme
TMPLG inline void TMPL_Domain::pushDissipationLangevin(double time, double gamma, double beta) {
  grid->pushDissipationLangevin(time, gamma, beta);
}


/// @brief Make neighbor lists
TMPLG void TMPL_Domain::makeNeighborLists() {
  grid->checkLocalBuffers();
  grid->makeNeighborLists();
  grid->sortNeighbors(Traversal::ALL);
}

/// @brief Make neighbor lists inside the domain
TMPLG void TMPL_Domain::makeNeighborListsInside() {
  grid->checkLocalBuffers();
  grid->makeNeighborListsInside();
  grid->sortNeighbors(Traversal::INSIDE);

}

/// @brief Make neighbor lists on the domain edges
TMPLG void TMPL_Domain::makeNeighborListsOnEdges() {
  grid->checkLocalBuffers();
  grid->makeNeighborListsOnEdges();
  grid->sortNeighbors(Traversal::EDGE);

}

/// @brief Clear neighbor lists
TMPLG void TMPL_Domain::refineCells() {
  grid->refineCells();
}

/// @brief Clear neighbor lists
TMPLG void TMPL_Domain::clearNeighborLists() {
  grid->clearNeighborList();
}


extern bool updateVerletLists; ///< boolean value to need if the Verlet list will be updated ot not

/// @brief Initialize forces, energies and neighbors
///
/// @param doInitEint If true, initialize the internal energies
/// @param [in] initialTint Initial internal temperature
TMPLG void TMPL_Domain::initForces(bool doInitEint, double initialTint) {

	clearGhost();
	grid->clearNeighborList();
	grid->refineCells();

	updateVerletLists = true;
	updateGhost();
	collectGhost();
	grid->makeNeighborLists();
	doComputeForcesVerlet();

	updateVerletLists = false;
}

/// @brief Compute forces
///
///
TMPLG void TMPL_Domain::doComputeForcesVerlet() {

  grid->resetForceAndEnergy();
  
  computeForcesVerlet(Traversal::ALL);

}


/// @brief execute checks
///
///
TMPLG void TMPL_Domain::ctest(std::set<ParticleReferenceValue> &reference_values_set, double &ae, double &re) {

  grid->ctest(reference_values_set, ae, re);
}



/// @brief Compute force on inside cells with Verlet Lists
TMPLG inline void TMPL_Domain::doComputeForcesInsideVerlet() {

  grid->resetForceAndEnergy();
  computeForcesVerlet(Traversal::INSIDE);

}


/// @brief Compute force on edge cells
TMPLG inline void TMPL_Domain::doComputeForcesOnEdgesVerlet() {

  computeForcesVerlet(Traversal::EDGE);

}



/// @brief Exchange particles between the domains and update cells
TMPLG void TMPL_Domain::updateCells() {

  exchangeParticles();
  internalReorganization();
  collectParticles();

  grid->updateGhost();
  grid->collectGhost();
  
}



/// @brief Update cells
TMPLG void TMPL_Domain::updateCellsInside() {

  exchangeParticles();
  internalReorganization();
  collectParticles();
  updateVerletLists = true;
  
}


/// @brief Get the number of real particles on the domain
/// @return Number of particles
TMPLG inline uint64_t TMPL_Domain::getNumberOfParticles() {
  return grid->getNumberOfParticles();
}


/// @brief Get the number of real cells on the domain
/// @return Number of particles
TMPLG inline uint64_t TMPL_Domain::getNumberOfRealCells() {
  return grid->getNumberOfRealCells();
}


/// @brief Get the total energy (kinetic plus potential) of the particles on the domain
/// @return Total energy
TMPLG inline double TMPL_Domain::getTotalEnergy() {
  return grid->getTotalEnergy(); 
}


/// @brief Get the kinetic energy of the particles on the domain
/// @return Kinetic energy
TMPLG inline double TMPL_Domain::getKineticEnergy() {
  return grid->getKineticEnergy(); 
}


/// @brief Get the potential energy of the particles on the domain
/// @return Potential energy
TMPLG inline double TMPL_Domain::getPotentialEnergy() {
  return grid->getPotentialEnergy(); 
}


/// @brief Get internal energy
TMPLG inline double TMPL_Domain::getInternalEnergy() {
  return grid->getInternalEnergy(); 
}


/// @brief Get the chemical energy of the particles on the domain
/// @return chemical energy
TMPLG inline double TMPL_Domain::getChemicalEnergy() {
  return grid->getChemicalEnergy();
}


/// @brief Get the total momentum of the particles on the domain
/// @return Total momentum
TMPLG inline vec3<double> TMPL_Domain::getTotalMomentum() {
  return grid->getTotalMomentum();
}


/// @brief Get the total mass of the particles on the domain
/// @return Total mass
TMPLG inline double TMPL_Domain::getTotalMass() {
  return grid->getTotalMass();
}


/// @brief Get the total kinetic energy on the domain in the center of momentum frame
///
/// @param [in] vshift Global velocity of the system
/// @return Kinetic energy in the center of momentum frame
TMPLG inline double TMPL_Domain::getShiftedKineticEnergy(const vec3<double>& vshift) {
  return grid->getShiftedKineticEnergy(vshift);
}


/// @brief Get the total pressure tensor on the domain
///
/// @param [in] vshift Global velocity of the system
/// @return Pressure tensor
TMPLG inline mat3<double> TMPL_Domain::getPressure(const vec3<double>& vshift) {
  return grid->getPressure(vshift);
}


/// @brief Accessor to ref (communication Manager)
TMPLG inline CommManager* TMPL_Domain::getCommManager() {
  return ref;
}


/// @brief Accessor to I/O manager
TMPLG inline InputOutputManager* TMPL_Domain::getInputOutputManager(){
  return m_inputOutputManager;
}


/// @brief Write domain data at current step in specified flux (not used)
///
/// Write domain index, number of particles and energy
/// @param [in,out] flux Print flux
TMPLG void TMPL_Domain::writeStep(std::ostream& flux) {

  // Basic : index -- number of particles -- energy
  flux<< "  " << "Domain " 
      << std::setw(4) << getIndex() << " | " 
      << std::setw(6) << getNumberOfParticles() << " particle(s) | "
      << "energy " << std::setprecision(16) << getTotalEnergy()
      << std::endl;

}


/// @brief Check if there is escaped particles at free boundaries
/// @param [out] outOfFreeBounds Presence of escaped particles for each free boundary
TMPLG inline void TMPL_Domain::checkFreeBoundaries(Array<int>& outOfFreeBounds) {
  grid->checkFreeBoundaries(outOfFreeBounds);
}


/// @brief Stop the wall particles
TMPLG inline void TMPL_Domain::stopWalls() {
  grid->stopWalls();
}


/// @brief Fill a buffer with all the particles of the node, case of a ParticleOutput buffer
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Domain::fillBuffer(ParticleOutput* buffer) {
  grid->fillBuffer(buffer);
}

/// @brief Fill a buffer with all the particles of the node, case of a ParticleInSitu buffer
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Domain::fillBuffer(ParticleInSitu* buffer) {
  grid->fillBuffer(buffer);
}


/// @brief Fill a buffer with the ghost information of the nodes, case of a ParticleInSitu buffer
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Domain::fillGhostBuffer(ParticleInSitu* buffer) {
  grid->fillGhostBuffer(buffer);
}


/// @brief Fill a buffer with all the particles of the node, case of a legacy Stamp dump buffer
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Domain::fillBuffer(Array<LegacyParticleIOStruct>& buffer) {
  grid->fillBuffer(buffer);
}


/// @brief Fill a buffer with all the particles of the node, case of a legacy Stamp dump buffer for DPDE
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Domain::fillBuffer(Array<LegacyDPDEParticleIOStruct>& buffer) {
  grid->fillBuffer(buffer);
}


/// @brief Fill a buffer with all the particles of the node, case of a Hercule dump buffer
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Domain::fillBuffer(HerculeParticleIODumpStruct& buffer) {
  grid->fillBuffer(buffer);
}


/// @brief Fill a buffer with all the particles of the node, case of a Hercule dump buffer for DPDE
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Domain::fillBuffer(HerculeDPDEParticleIODumpStruct& buffer) {
  grid->fillBuffer(buffer);
}


/// @brief Fill a buffer with all the cells of the node
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Domain::fillBuffer(CellOutput* buffer) {
  grid->fillBuffer(buffer);
}


/// @brief Print domains local info (like symmetrization) in specified flux
/// @param [in,out] flux Print flux
TMPLG inline void TMPL_Domain::printInfoBase(std::ostream& flux) {
  LocalInfo::print(flux);
}


/// @brief Print grid info in specified flux
/// Call non implemented or useless function (not used)
/// @param [in,out] flux Print flux
TMPLG inline void TMPL_Domain::printInfoSpec(std::ostream& flux) {
  grid->printInfo(flux);
}


/// @brief Debug print for all cells in specified flux (not used)
/// @param [in,out] flux Print flux
TMPLG inline void TMPL_Domain::print(std::ostream& flux) {
  grid->print(flux);
  ref->barrier();
}


/// @brief Move particles inside domain
///
///
TMPLG inline void TMPL_Domain::internalReorganization() {
  grid->internalReorganization();
}


/// @brief Move particles between domains
///
///
TMPLG inline void TMPL_Domain::exchangeParticles() {
  grid->exchangeParticles();
}


/// @brief Collect particles after a communication
///
///
TMPLG inline void TMPL_Domain::collectParticles() {
  grid->collectParticles();
}


/// @brief Update ghost cells
///
///
TMPLG inline void TMPL_Domain::updateGhost() {
  grid->updateGhost();
}


/// @brief Clear the ghost cells
///
///
TMPLG inline void TMPL_Domain::clearGhost() {
	grid->clearGhost();
}


/// @brief Collect ghost after a communication
///
///
TMPLG inline void TMPL_Domain::collectGhost() {
  grid->collectGhost();
}


/// @brief Compute forces on a set of cells given by traversal
/// @param [in] traversal Indicates on which cells forces must be computed
TMPLG void TMPL_Domain::computeForces(Traversal traversal) {
	grid->computeForces(traversal);
}

/// @brief Compute forces on a set of cells given by traversal with the Verlet lists
/// @param [in] traversal Indicates on which cells forces must be computed
TMPLG void TMPL_Domain::computeForcesVerlet(Traversal traversal) {
	grid->computeForcesVerlet(traversal);
}


/// @brief Compute workload on each particle
/// @return Total workload on the domain
TMPLG inline double TMPL_Domain::workload() {
  return grid->workload();
}


/// @brief Set Zoltan load balancing parameters
/// @tparam Grid_impl The proper Grid subclass
/// @param [in,out] loadBalancer Load balancer (pointer)
TMPLG inline void TMPL_Domain::setCallbackQueryFunctions(LBS::LoadBalancer* loadBalancer) {

#if __use_lib_zoltan

  auto ptr = loadBalancer->zoltanStruct();

  Zoltan_Set_Num_Obj_Fn    (ptr, Grid_impl::numberOfCells, grid);
  Zoltan_Set_Obj_List_Fn   (ptr, Grid_impl::listObj      , grid);
  Zoltan_Set_Num_Geom_Fn   (ptr, Grid_impl::numberOfDim  , grid);
  Zoltan_Set_Geom_Fn       (ptr, Grid_impl::fillCoords   , grid);
  Zoltan_Set_Num_Edges_Fn  (ptr, Grid_impl::numberOfEdges, grid);
  Zoltan_Set_Edge_List_Fn  (ptr, Grid_impl::fillEdges    , grid);

#endif

}


/// @brief Balance workload on domains
/// @tparam Grid_impl The proper Grid subclass
/// @param [in] loadBalancer Load balancer (pointer)
/// @return State of the load balancer at the end
TMPLG LBS::State TMPL_Domain::balance(LBS::LoadBalancer* loadBalancer) {

#if __use_lib_zoltan

  /// @brief Shortcut for exchange_t type (from Grid_impl)
  typedef typename Grid_impl::exchange_t exchange_t;

  Zoltan_Struct* ptr = loadBalancer->zoltanStruct();

  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
  int *importProcs, *importToPart, *exportProcs, *exportToPart;

  Zoltan_LB_Partition(ptr,        // input (all remaining fields are output)
		      &changes,           // 1 if partitioning was changed, 0 otherwise 
		      &numGidEntries,     // Number of integers used for a global ID
		      &numLidEntries,     // Number of integers used for a local ID
		      &numImport,         // Number of vertices to be sent to me
		      &importGlobalGids,  // Global IDs of vertices to be sent to me
		      &importLocalGids,   // Local IDs of vertices to be sent to me
		      &importProcs,       // Process rank for source of each incoming vertex
		      &importToPart,      // New partition for each incoming vertex
		      &numExport,         // Number of vertices I must send to other processes
		      &exportGlobalGids,  // Global IDs of the vertices I must send
		      &exportLocalGids,   // Local IDs of the vertices I must send
		      &exportProcs,       // Process to which I send each of the vertices
		      &exportToPart);     // Partition to which each vertex will belong 
  
  if (changes) {

    //updateCells();
 
    // amr
    m_balanceComm.clear();
    m_toKeep.clear();
    
    // (0) Get domains I'll communicate with
    Array<uint> procsToComm;
    getProcToComm(procsToComm, importProcs, numImport, exportProcs, numExport);
    
    // (*) Dump particles to keep and particles to send, send leaving particles
    
    grid->dumpParticles(m_toKeep, m_balanceComm, exportLocalGids, exportProcs, numExport);

    // While sending particles, get other stuff done


    // (*)
    Array<int> allCellsOwners;
    grid->updateOwners(allCellsOwners, numExport, exportLocalGids, exportProcs, numImport, importGlobalGids);
    
    for(size_t i = 0 ; i < allCellsOwners.size() ; i++)
	    assert(allCellsOwners[i] >= 0 );

    // (*)
    std::vector<uint32_t> importEnv;
    grid->exchangeMovingCellEnv(procsToComm, allCellsOwners, numImport, importEnv, importGlobalGids, numExport, exportGlobalGids, exportLocalGids, exportProcs);
    
    // (*) Get zoltan content in vector

    std::vector<uint> importIndexes(numImport);
    parallel_region(0, numImport, [&](const uint begin, const uint end) {
	    for(uint i=begin; i<end; ++i) 
	      importIndexes[i] = importGlobalGids[i];
      });
    std::stable_sort(importIndexes.begin(), importIndexes.end());

    std::vector<uint> exportIndexes(numExport);
    parallel_region(0, numExport, [&](const uint begin, const uint end) {
	    for(uint i=begin; i<end; ++i) 
	      exportIndexes[i] = exportGlobalGids[i];
    });
    std::stable_sort(exportIndexes.begin(), exportIndexes.end());

    // (*) Free Zoltan Stuff 
    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, &exportProcs, &exportToPart);

    // (*) Update decomposition
    static_cast<AnyDecomposition*>(Global::domainInfo.getDecomposition())->update_part_1(importIndexes, exportIndexes, importEnv);
    importIndexes.clear();
    importEnv.clear();
    exportIndexes.clear();
    static_cast<Grid_impl*>(grid)->updateDecomposition(allCellsOwners);

    // (*) Build new empty grid, put particles back
    delete grid;
    grid = new Grid_impl(this, this->index, 0);
    myreintegrateParticles();
    
    // clear
    m_balanceComm.clear();
    m_toKeep.clear();
    
        
    //build of amr and neighbors lists --> needed for imbalance computation
  /*  updateVerletLists = true;
    grid->refineCells();
    
    updateGhost();
	  collectGhost();
    
    grid->checkLocalBuffers();
  	grid->clearNeighborList();
  	grid->makeNeighborLists();

  	
  	grid->sortNeighbors(Traversal::ALL);
  	updateVerletLists = false;*/

    
    return LBS::BALANCE;

  }
  else {
  
    return LBS::NOTHING_TO_DO;

  }

#else

  return LBS::VOID;

#endif

}


/// @brief Gather the ranks of the domains that exchange (import or export) particles with the current domain during load balancing
/// @param [out] procsToComm Ranks of the domains with which current domain exchanges particles
/// @param [in] importProcs Ranks of the domains from which particles must be imported
/// @param [in] numImport Number of domains for import
/// @param [in] exportProcs Ranks of the domains to which particles must be exported
/// @param [in] numExport Number of domains for export
TMPLG void TMPL_Domain::getProcToComm(Array<uint>& procsToComm, int* importProcs, const uint numImport, int* exportProcs, const uint numExport) {

  std::set<uint> tmp;
      
  for (uint i=0; i<numImport; ++i) tmp.insert(importProcs[i]);
  for (uint i=0; i<numExport; ++i) tmp.insert(exportProcs[i]);
      
  procsToComm = Array<uint>(tmp.size());

  uint nxt = 0;

  for (auto& elem : tmp) 
    procsToComm[nxt++] = elem;

}


/// @brief Reintegrate received and not moved particles to the grid
/// @tparam Grid_impl The proper Grid subclass
TMPLG void TMPL_Domain::myreintegrateParticles() {

  static_cast<Grid_impl*>(grid)->unPackBalance(m_balanceComm);
  static_cast<Grid_impl*>(grid)->unPackBalanceKeep(m_toKeep);
  
  this->initForces(false,0.);
}

/// @brief Determines whether an atom has moved more than 1/2 of the verlet radius.
/// @return Indicates if the neighbor lists must be updated
TMPLG inline bool TMPL_Domain::checkVerlet() {

  return grid->checkVerlet();
  
}

/// @brief Determines whether an atom has moved more than 1/2 of the verlet radius.
/// @return Indicates if the neighbor lists must be updated
TMPLG inline void TMPL_Domain::refineBalance() {

}
  
