/// @file
/// @brief Implementation of class Grid templated methods

/// @brief Default constructor
///
///
TMPLG TMPL_Grid::Grid()
  : numberOfParticles(0), totalNumberOfParticles(0),
    energyFlag(false), energyKin(0.), energyPot(0.),
    energyInt(0.), energyChm(0.), 
    m_workloadFlag(false), m_workload(0.),
    correcter(), modD(), modI(), 
    tManager(),
    domain(nullptr),
    mtx(), coords(), neighbors(), ghostLayer(),
    raphael_test(),
    destDomains(), destCells()
    {}


/// @brief Destructor
///
///
TMPLG TMPL_Grid::~Grid() {
}


Array< vec3<int> > __getMovList(const vec3<int>& mov);


/// @brief Update particles positions using to first order approximation of the movement equations
///
/// Call specialized function from Grid_impl
/// @param [in] time Time step
TMPLG inline uint64_t TMPL_Grid::getTotalNbOfParticlesSend() {
  return impl().getTotalNbOfParticlesSend();
}

/// @brief Update particles positions using to first order approximation of the movement equations
///
/// Call specialized function from Grid_impl
/// @param [in] time Time step
TMPLG inline uint64_t TMPL_Grid::getTotalNbOfParticlesRecv() {
  return impl().getTotalNbOfParticlesRecv();
}

/// @brief Update particles positions using to first order approximation of the movement equations
///
/// Call specialized function from Grid_impl
/// @param [in] time Time step
TMPLG inline void TMPL_Grid::pushPositions1stOrder(double time) {
  impl().pushPositions1stOrder(time);
}



/// @brief Update particles positions using to second order approximation of the movement equations
///
/// Call specialized function from Grid_impl
/// @param [in] time Time step
TMPLG inline void TMPL_Grid::pushPositions2ndOrder(double time) {
  impl().pushPositions2ndOrder(time);
}


/// @brief Update particles velocities using to first order approximation of the movement equations
///
/// Call specialized function from Grid_impl
/// @param [in] time Time step
TMPLG inline void TMPL_Grid::pushVelocities1stOrder(double time) {
  impl().pushVelocities1stOrder(time);
  energyFlag = false;
}

/// @brief Update particles velocities by applying Langevin dissipation
///
/// Call specialized function from Grid_impl
/// @param [in] time Time step
/// @param [in] gamma Gamma parameter of the Langevin scheme
/// @param [in] beta Beta parameter of the Langevin scheme
TMPLG inline void TMPL_Grid::pushDissipationLangevin(double time, double gamma, double beta) {
  impl().pushDissipationLangevin(time, gamma, beta);
  energyFlag = false;
}

/// @brief Reset forces of all the particles of the grid
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::resetForce() {
  impl().resetForce();
}


/// @brief Reset fluctuation forces of all the particles of the grid
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::resetForceFluct() {
  impl().resetForceFluct();
}


/// @brief Reset energies of all the particles of the grid
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::resetEnergy() {
  impl().resetEnergy();
  energyFlag = false;
}


/// @brief Reset forces and energies of all the particles of the grid
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::resetForceAndEnergy() {
  impl().resetForceAndEnergy();
  energyFlag = false;
}


/// @brief Reset densities of all the particles of the grid
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::resetDensity() {
  impl().resetDensity();
}


/// @brief Reset reaction evolution of all the particles of the grid
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::resetChemistry() {
  impl().resetChemistry();
}


/// @brief Make the neighbors lists for all the particles
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::makeNeighborLists() {
  impl().makeNeighborLists();
}


/// @brief Make the neighbors lists for the inside cells
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::makeNeighborListsInside() {
  impl().makeNeighborListsInside();
}


/// @brief Make the neighbors lists for the edge cells
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::makeNeighborListsOnEdges() {
  impl().makeNeighborListsOnEdges();
}


/// @brief Reset neighbors lists
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::clearNeighborList() {
  impl().clearNeighborList();
}


/// @brief Sort neighbor list for each particle of a traversal
///
/// Call specialized function from Grid_impl
/// @param [in] T Traversal considered
TMPLG inline void TMPL_Grid::sortNeighbors(::Traversal T) {
  impl().sortNeighbors(T);
}



/// @brief Compute all the forces on the specified traversal (potentials and boundary conditions)
/// @param [in] T Base traversal specifying where the forces are computed
TMPLG inline void TMPL_Grid::computeForcesVerlet(::Traversal T) {

  // Loop on all potentials : for each type of potential, cast it and compute
  // interaction between particles of type 'first' and type 'second'
  Global::reference.for_all_potentials([&](TypedPotential pot){

      switch(pot.traversal) {

      case IPotential::IDEAL_GAS : // ideal gas = no interactions
	break;

      case IPotential::PAIR :
      	impl().computeForceVerlet(static_cast<PairPotential*>(pot.potential),
			  pot.first, pot.second, T);
	break;

      case IPotential::EAM :
      	impl().computeForceVerlet(static_cast<EAMPotential*>(pot.potential),
			  pot.first, pot.second, T);
	break;

      }

    });
  

  // Add boundary conditions
  const auto& bcs = Global::domainInfo.getBoundaryConditions();
  for (uint dim=0; dim<VEC3_NDIMS; ++dim) {

    switch (bcs[dim]) {

    case Configuration<DomainInterface>::FREE:
      // some day ...
      break;

    case Configuration<DomainInterface>::PERIODIC:
      // nothing to do here
      break;

    case Configuration<DomainInterface>::WALL:
      this->applyWallCondition(dim, -1);
      this->applyWallCondition(dim, +1);
      break;

    case Configuration<DomainInterface>::FREE_WALL:
      this->applyWallCondition(dim, +1);
      break;

    case Configuration<DomainInterface>::WALL_FREE:
      this->applyWallCondition(dim, -1);
      break;

    }
  }

}


/// @brief Get moving particles and send them to the other domains
///
/// Call specialized functions from Grid_impl
TMPLG inline void TMPL_Grid::exchangeParticles() {

  energyFlag = false;

  impl().fillWithLeavingParticles();
  impl().sendLeavingParticles();

}


/// @brief Collect and add received particles
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::collectParticles() {

  energyFlag = false;

  impl().collectParticles();

}

/// @brief Update ghost to send to other domains
///
/// Call specialized functions from Grid_impl
TMPLG inline void TMPL_Grid::refineCells() {
  impl().refineCells();
}

/// @brief Update ghost to send to other domains
///
/// Call specialized functions from Grid_impl
TMPLG inline void TMPL_Grid::updateGhost() {
  impl().fillWithGhosts();
  impl().sendGhosts();
}

/// @brief Clear the ghost
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::clearGhost() {
  impl().clearGhost();
}


/// @brief Collect received ghost cells
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::collectGhost() {
  impl().collectGhost();
}


/// @brief Reorganize particles between the cells of the domain
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::internalReorganization() {
  impl().internalReorganization();
  energyFlag = false;
}


/// @brief Resize local buffers
///
/// Call specialized function from Grid_impl
TMPLG inline void TMPL_Grid::checkLocalBuffers() {
  impl().checkLocalBuffers();
}


/// @brief Accessor to number of real point particles
TMPLG inline uint64_t TMPL_Grid::getNumberOfParticles() {
  return numberOfParticles;
}


/// @brief Return number of real cells
TMPLG inline uint64_t TMPL_Grid::getNumberOfRealCells() {
  return this->getTraversal(TraversalManager::REAL).size();
}


/// @brief Get total energy (kinetic, potential, internal and chemical) of the particles (compute if necessary)
/// @return Total energy
TMPLG inline double TMPL_Grid::getTotalEnergy() {
  if (!energyFlag) {
    impl().computeEnergy();
    energyFlag = true;
  }
  return energyKin + energyPot + energyInt + energyChm;
}


/// @brief Get total kinetic energy of the particles (compute if necessary)
/// @return Kinetic energy
TMPLG inline double TMPL_Grid::getKineticEnergy() {
  if (!energyFlag) {
    impl().computeEnergy();
    energyFlag = true;
  }
  return energyKin;
}


/// @brief Get total potential energy of the particles (compute if necessary)
/// @return Potential energy
TMPLG inline double TMPL_Grid::getPotentialEnergy() {
  if (!energyFlag) {
    impl().computeEnergy();
    energyFlag = true;
  }
  return energyPot;
}


/// @brief Get internal energy of the particles (compute if necessary)
/// @return Internal energy
TMPLG inline double TMPL_Grid::getInternalEnergy() {
  if (!energyFlag) {
    impl().computeEnergy();
    energyFlag = true;
  }
  return energyInt;
}


/// @brief Get chemical energy of the particles (compute if necessary)
/// @return Chemical energy
TMPLG inline double TMPL_Grid::getChemicalEnergy() {
  if (!energyFlag) {
    impl().computeEnergy();
    energyFlag = true;
  }
  return energyChm;
}


/// @brief Get workload on the domain (compute if necessary)
/// @return Workload
TMPLG inline double TMPL_Grid::workload() {

//  if (!m_workloadFlag) {
    impl().computeWorkload();
//    m_workloadFlag = true;
//  }
    return m_workload;
}


/// @brief Get total momentum of the particles
///
/// Call specialized function from Grid_impl
/// @return Momentum
TMPLG inline vec3<double> TMPL_Grid::getTotalMomentum() {
  return impl().computeTotalMomentum();
}


/// @brief @brief Get total mass of the particles
///
/// Call specialized function from Grid_impl
/// @return Mass
TMPLG inline double TMPL_Grid::getTotalMass() {
  return impl().computeTotalMass();
}


/// @brief Get the total kinetic energy in the center of momentum frame
///
/// Call specialized function from Grid_impl
/// @param [in] vshift Global velocity of the system
/// @return Kinetic energy in the center of momentum frame
TMPLG inline double TMPL_Grid::getShiftedKineticEnergy(const vec3<double>& vshift) {
  return impl().computeShiftedKineticEnergy(vshift);
}

/// @brief Get the pressure tensor
///
/// Call specialized function from Grid_impl
/// @param [in] vshift Global velocity of the system
/// @return Pressure tensor
TMPLG inline mat3<double> TMPL_Grid::getPressure(const vec3<double>& vshift) {
  return impl().computePressure(vshift);
}


/// @brief Get a chunk of particles from the configuration and initialize those that are in the grid
///
/// Call specialized function from Grid_impl
/// @param [in] particles Configuration of the particles
/// @param [in] cBegin Start of the set
/// @param [in] cEnd End of the set
TMPLG inline void TMPL_Grid::initParticles(Configuration<Particle>& particles, uint64_t cBegin, uint64_t cEnd) {
  impl().initParticles(particles, cBegin, cEnd);
  energyFlag = false;
}


/// @brief Get the particles from an array of particle pointers and initialize those in the grid
///
/// Call specialized function from Grid_impl
/// @param [in] particles Pointer to the system particles
/// @param [in] initialT Initial temperature
/// @param [in] numP Number of particles
/// @param [in] nCells Number of duplicates of the base cell in each dimension
/// @param [in] cellSize Size of the base cell
/// @param [in] numMol USELESS FOR AMR
TMPLG inline void TMPL_Grid::initParticles(MPI__Particle** particles, double initialT, uint64_t numP, vec3<int> nCells, vec3<double> cellSize, uint64_t numMol) {
  energyFlag = false;
}


/// @brief Apply wall conditions to the particles on specified boundary
///
/// Call specialized function from Grid_impl
/// @param [in] dim Dimension of the wall : 0 for x, 1 for y, 2 for z
/// @param [in] direction Direction of the wall (-1 for inf bound, 1 for sup bound)
TMPLG inline void TMPL_Grid::applyWallCondition(const int dim, const int direction) {
  impl().applyWallCondition(dim, direction);
}


/// @brief Check if there is escaped particles at free boundaries
///
/// Call specialized function from Grid_impl
/// @param [out] outOfFreeBounds Presence of escaped particles for each free boundary
TMPLG inline void TMPL_Grid::checkFreeBoundaries(Array<int>& outOfFreeBounds) {
  impl().checkFreeBoundaries(outOfFreeBounds);
}


/// @brief Stop the wall particles
TMPLG inline void TMPL_Grid::stopWalls() {

  impl().stopWalls();
  
}

/// @brief Fill a buffer with all the particles of the node, case of a ParticleOutput buffer
///
/// Call specialized function from Grid_impl
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Grid::fillBuffer(ParticleOutput* buffer) {
  impl().fillBuffer(buffer);
}

/// @brief Fill a buffer with all the particles of the node, case of a ParticleInSitu buffer
///
/// Call specialized function from Grid_impl
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Grid::fillBuffer(ParticleInSitu* buffer) {
  impl().fillBuffer(buffer);
}


/// @brief Fill a buffer with the ghost informations of the node, case of a ParticleInSitu buffer
///
/// Call specialized function from Grid_impl
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Grid::fillGhostBuffer(ParticleInSitu* buffer) {
  impl().fillGhostBuffer(buffer);
}


/// @brief Fill a buffer with all the particles of the node, case of a legacy Stamp dump buffer
///
/// Call specialized function from Grid_impl
/// @tparam DumpStruct Dump structure (legacy or legacy for DDPE)
/// @brief [out] Buffer to fill
TMPLG template <class DumpStruct>
inline void TMPL_Grid::fillBuffer(Array<DumpStruct>& buffer) {
  impl().fillBuffer(buffer);
}


/// @brief Fill a buffer with all the particles of the node, case of a Hercule dump buffer
///
/// Call specialized function from Grid_impl
/// @tparam DumpStruct Dump structure (Hercule or Hercule for DDPE)
/// @brief [out] Buffer to fill
TMPLG template <class DumpStruct>
inline void TMPL_Grid::fillBuffer(DumpStruct& buffer) {
  impl().fillBuffer(buffer);
}


/// @brief Fill a buffer with all the cells of the node
///
/// Call specialized function from Grid_impl
/// @brief [out] Buffer to fill
TMPLG inline void TMPL_Grid::fillBuffer(CellOutput* buffer) {
  impl().fillBuffer(buffer);
}


/// @brief Update the owners of real and ghost cells after load balancing
///
/// Call specialized function from Grid_impl
/// @param [out] allCellsOwners Cells owners
/// @param [in] numExport Number of cells exported to other domains
/// @param [in] exportLocalGids Local ID of the cells exported to other domains
/// @param [in] exportProcs Domains where the cells are exported
/// @param [in] numImport Number of cells to import from other domains
/// @param [in] importGlobalGids Global ID of the cells imported from other domains
TMPLG inline void TMPL_Grid::updateOwners(Array<int>& allCellsOwners, int numExport, uint* exportLocalGids, int* exportProcs, int numImport, uint* importGlobalGids) {
  impl().updateOwners(allCellsOwners, numExport, exportLocalGids, exportProcs, numImport, importGlobalGids);
}


/// @brief Exchange neighbors of displaced cells
///
/// Call specialized function from Grid_impl
/// @param [in] procsToComm Ranks of the domains with which current domain exchanges particles
/// @param [in] allCellsOwners Cells owners
/// @param [in] numImport Number of cells to import from other domains
/// @param [out] importEnv Neighbors of imported cells
/// @param [in] importGlobalGids Global ID of the cells imported from other domains (not used)
/// @param [in] numExport Number of cells exported to other domains
/// @param [in] exportGlobalGids Global ID of the cells exported to other domains
/// @param [in] exportLocalGids Local ID of the cells exported to other domains
/// @param [in] exportProcs Domains where the cells are exported
TMPLG inline void TMPL_Grid::exchangeMovingCellEnv(const Array<uint>& procsToComm, const Array<int>& allCellsOwners, int numImport, std::vector<uint32_t>& importEnv, uint* importGlobalGids, int numExport, uint* exportGlobalGids, uint* exportLocalGids, int* exportProcs) {
  impl().exchangeMovingCellEnv(procsToComm, allCellsOwners, numImport, importEnv, importGlobalGids, numExport, exportGlobalGids, exportLocalGids, exportProcs);
}


/// @brief Print grid info in specified flux (debug)
///
/// Call specialized function from Grid_impl
/// @warning Not implemented for SingleSpecGrid
/// @param [in,out] flux Print flux
TMPLG void TMPL_Grid::printInfo(std::ostream& flux) {
  impl().printInfo(flux);
}


/// @brief Debug print for all cells in specified flux
///
/// Call specialized function from Grid_impl
/// @param [in,out] flux Print flux
TMPLG void TMPL_Grid::print(std::ostream& flux) {
  impl().print(flux);
}


/// @brief Update the velocities of the ghost particles
TMPLG inline void TMPL_Grid::updateGhostVel() {
  impl().updateGhostVel();
}


/// @brief Update the fluctuation forces of the ghost particles
TMPLG inline void TMPL_Grid::updateGhostFluct() {
  impl().updateGhostFluct();
}


/// @brief Update the densities of the ghost particles
TMPLG inline void TMPL_Grid::updateDensities() {
  impl().updateDensities();
}


/// @brief Get index from the domain
/// @return Domain index
TMPLG inline uint TMPL_Grid::getDomainIndex() {
  return domain->getIndex();
}


/// @brief Initialize all the traversals from the grid info (case AnyGrid)
/// @param [in] info Grid info
/// @param [in] coords Coordinates of the cells in the grid
TMPLG inline void TMPL_Grid::makeTraversals(AnyGridInfo* info, const std::vector< vec3<int> >& coords) {
  tManager.makeTraversals(info, coords);
}


/// @brief Initialize all the traversals from the grid info (case RectilinearGrid)
/// @param [in] info Grid info
/// @param [in] coords Coordinates of the cells in the grid
TMPLG inline void TMPL_Grid::makeTraversals(RectilinearGridInfo* info, const std::vector< vec3<int> >& coords) {
  tManager.makeTraversals(info, coords);
}


/// @brief Get the cells in specified traversal
/// @param [in] traversal Traversal
/// @return Cells of the traversal
TMPLG inline const Array<uint>& TMPL_Grid::getTraversal(TraversalManager::CellTraversal traversal) {
  return tManager.getTraversal(traversal);
}


/// @brief Identify traversal and get the cells from a base traversal and a level
/// @param [in] T Base traversal
/// @param [in] level Level depending on ghost thickness used for the identification, dafault=0
/// @return Cells of the traversal
TMPLG inline const Array<uint>& TMPL_Grid::getTraversal(::Traversal T, uint level) {
  return tManager.getTraversal(tManager.getCellTraversal(T, level));
}


/// @brief Correct a distance to account for boundary conditions
/// @param [in,out] distance Distance to correct
TMPLG inline void TMPL_Grid::correct(vec3<double>& distance) {
  return correcter.correct(distance);
}


/// @brief Correct multiple distances to account for boundary conditions
/// @param [in] size Number of distances to correct
/// @param [in,out] rx Pointer to x components
/// @param [in,out] ry Pointer to y components
/// @param [in,out] rz Pointer to z components
TMPLG inline void TMPL_Grid::correct(uint size, double* rx, double* ry, double* rz) {
  correcter.correct(size, rx, ry, rz);
}


/// @brief Modulate 3D double coordinates to set them back in the system bounds
/// @param [in,out] distance Coordinates to modulate
TMPLG inline void TMPL_Grid::module(vec3<double>& distance) {
  modD.module(distance);
}


/// @brief Modulate 3D integer coordinates to set them  back in the system bounds
/// @param [in,out] distance Coordinates to modulate
TMPLG inline void TMPL_Grid::module(vec3<int>& distance) {
  modI.module(distance);
}


/// @brief Modulate 3 one dimensional double coordinates to set them  back in the system bounds (not used)
/// @param [in,out] x X component to modulate
/// @param [in,out] y Y component to modulate
/// @param [in,out] z Z component to modulate
TMPLG inline void TMPL_Grid::module(double& x, double& y, double& z) {
  modD.module(x, y, z);
}


/// @brief Lock specified mutex
/// @param [in] i Index of the mutex to lock
TMPLG inline void TMPL_Grid::lock(uint i) {
  mtx.lock(i);
}


/// @brief Unlock specified mutex
/// @param [in] i Index of the mutex to unlock
TMPLG inline void TMPL_Grid::unlock(uint i) {
  mtx.unlock(i);
}


/// @brief Check if specified cell is in ghost
/// @param [in] i Cell index
/// @return True if in ghost
TMPLG inline bool TMPL_Grid::isGhost(uint i) {
  return ghostLayer[i] > 0;
}


/// @brief Check if specified cell is real
/// @param [in] i Cell index
/// @return True if real
TMPLG inline bool TMPL_Grid::isReal(uint i) {
  return ghostLayer[i] == 0;
}


/// @brief Call the proper Grid subclass
/// @tparam Grid_impl The proper Grid subclass
TMPLG inline Grid_impl& TMPL_Grid::impl() {
  return reinterpret_cast<Grid_impl&>(*this);
}


/// @brief Determines whether an atom has moved more than 1/2 of the verlet radius.
///
/// Call specialized function from Grid_impl
/// @return True if verlet list must be updated
TMPLG inline bool TMPL_Grid::checkVerlet() {
  return impl().checkVerlet();
}

TMPLG inline void TMPL_Grid::refineBalance() {
  return impl().Refine();
}

