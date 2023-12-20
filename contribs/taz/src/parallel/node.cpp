/// @file 
/// @brief Implementation of class Node


#include <fstream>
#include <iomanip>
#include <vector>

#include "parallel/node.hpp"

#include "particle/molecule/configInMol.hpp"

#include "time/time.hpp"
#include "time/timeIntegration.hpp"

#ifndef __xsp_init_slice
/// @brief Defines default slice value
#define __xsp_init_slice 1048576
#endif

class Metrics;
Metrics * ptM;

/// @brief A constructor from the program input data
/// @param [in] input Input to convert
Configuration<Node>::Configuration(const Input& input) {

  switch (input.nodeType) {

  case Input::FLAT :
    type = FLAT;
    break;

  }

  maxNumberOfThreads = input.maxNumberOfThreadsPerNode;

  switch (input.balancingMethod) {

  case Input::NONE :
    balancingMethod = LBS::NONE;
    break;

  case Input::BLOCK :
    balancingMethod = LBS::BLOCK;
    break;

  case Input::RANDOM :
    balancingMethod = LBS::RANDOM;
    break;

  case Input::COORD_BSC :
    balancingMethod = LBS::REC_COORD_BSC;
    break;

  case Input::INERT_BSC :
    balancingMethod = LBS::REC_INERT_BSC;
    break;

  case Input::SPC_FILL_CRV :
    balancingMethod = LBS::SPACE_FILL_CURVE;
    break;

  case Input::PHG :
    balancingMethod = LBS::PHG;
    break;

  case Input::METIS :
    balancingMethod = LBS::PARMETIS;
    break;

  case Input::SCOTCH :
    balancingMethod = LBS::SCOTCH;
    break;

  default :
    balancingMethod = LBS::NONE;
    break;

  }

  logWorkload      = input.logWorkload;
  expandFreeLimits = input.expandFreeLimits;
  expandFactor     = input.expandFactor;

  Global::seed += input.seed;
  
}


/// @brief Constructor from a MPI communicator
/// @param [in,out] comm MPI communicator
Node::Node(MPI_Comm comm) 
  : setup(),metrics(),
    commManager(comm),
    isMaster(rank()==Global::masterNode),
    scheme(nullptr),
    loadBalancer(comm),
    m_outOfFreeBounds(6,0)
{

  metrics.tic(Metrics::INIT);

  // Header
  if (isMaster) {

  std::cout<<"\e[0;"<<34<<"m"<<"           "<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                                  $$     o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                                 $$$   o$   ooooooo"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                               $ $$$  o$ o$$$$$$$ooo"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                               $$$$$$o$$$$$$$$$$$$$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                          oooo$$$$$$$$$$'$$o$$o$$$'$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                         o$$$$$$$$o$$$o$$$$$$$$$$$$$$o      $"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                         '$$$$$$$'$$$$$$$$$$$$'$$$$$$$    o$'"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                     oooo $$$$$$$$$$$$''$$$$$$$$$$o$$$  o$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                  $$$$$$$$$$$$$$'$$$o$$   '$$$$$'$$$$$$$'$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                o$$$$$'$$$'$oo   ''$''      $$$$$$$$$$$$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                              o$$$$$'$$$''           'o oo   $$'$''''$$$'$$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                           o$$$$$'$$$$'       ''$oo      oo$$$$$o    $$$$$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                        oo$$$$$'$$$$$$$'      o$$'      o$$''  $$$o$$'''$$$$$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                     oo$$$$$'$$$$$$'$$'       $$$$      '$$$$$$$$$$       ''$$$'o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                 o$$$$$$$$'$$$$$$'$$$'       $$$'$       '''$$$$$'         $o$$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                $$$$$$'$o$$$'$$$$$$$$       $$$$o '$                        $$$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                 $$$o$$$$o   $$$o$$$        $$$$$o$$  o$   o  ' $o o$        $$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                 $ '$$$$$$o  $$$$$$$        $$$$$$$' o$$$ o$  o$      o' '  o$$$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                    '$$'$$$  $$$$o$$        $$$$$ $o$$$$$$$$$$$$o   '        $$$$'$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"              $$$$$$  $'$$$$$    'o$$$$$  'oo'''$$$$$$$$$$$o o   '  $$$$$$$$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"               '$$o$o $$$$$$$oo$$$$$$$$$$o '  '$' '$'''$$$$$$ooo$   $$$$$$$$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                $$$$$o$$$$o$$$$$$$$$$$$$$$ ''''' oo$   $$$$$$$$$$o $  ''$$$$o$$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                 '$$'$$$o$$$$$$$$o$$$  ''            ''$$$$$$$$$$$o$     '$$$$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                   '$$$$$$$$$$$$$$$$$$$o     o $ooooooo   '$$$$'$$$$$   o$$$'$$$'"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                     '$$$$$$$$''''$$$$$$oo$''''  $$$$$$$$o$$$$$$$$$$  o$$$$$$$$'"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                        $$$'$$o    ''$$$'        '$$$$$$$$$$$$$$$$$ o$$$$$$$''"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                 $$$$$$                     '   $$$$$$$$$$$$$$$''"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                    oo           '$$'$$                          '$$$$o$$$$''"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                     '$$o         $$$$$                           $ $$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                       $$$        $$$$$                         $$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                        $$$        $o$$                        $$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                         $$$o      $$$$                        $$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                        oo$$$o     $$$$                      $$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                        $$$$o$$$$$$$$$$                    $$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                        '$$$$$'$$$$$o$$                   $$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                          '$o$$$$$'$$$$                 o$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                            '$$$$$$$$$$o            ooo$$$$$$oo"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                            o$$$'$$$o$$$$oo    ooo$$$$$$$$$$$$$$oo"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                          o$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                         $$$$'$$$'$o$$$$'''           ''$$$$$o$$$$$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                        o$$$$$$$$$$$''                 o$$$'$$$$$$$$'"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                        ' $'$$'$$$$oo                o$$$$$$'''  '"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                              ''$$$$$$oo            $$$$$$oo"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"               oooo$$$$$oooo$$$$$$$o$$$$$           '''$$$$$$$$$$$$$$$$$$$$$$$oo"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"             o$$$$$$$$$$$$$$$$''''                               '''''$$$$$$$$$$$"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"            $$$$$$$$o$$$'''                                             ''o$$$$$$o"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"           o$$$$$$o$''                                                      '''''"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"             '''"<<"\e[0;0m"<<std::endl;              
  std::cout<<"\e[0;"<<34<<"m"<<"                                 ___________   _____   __________ "<<"\e[0;0m"<<std::endl; 
  std::cout<<"\e[0;"<<34<<"m"<<"                                 \\__    ___/  /  _  \\  \\____    / "<<"\e[0;0m"<<std::endl; 
  std::cout<<"\e[0;"<<34<<"m"<<"                                   |    |    /  /_\\  \\   /     /  "<<"\e[0;0m"<<std::endl; 
  std::cout<<"\e[0;"<<34<<"m"<<"                                   |    |   /    |    \\ /     /_  "<<"\e[0;0m"<<std::endl; 
  std::cout<<"\e[0;"<<34<<"m"<<"                                   |____|   \\____|__  //_______ \\ "<<"\e[0;0m"<<std::endl; 
  std::cout<<"\e[0;"<<34<<"m"<<"                                                    \\/         \\/ "<<"\e[0;0m"<<std::endl;
  std::cerr<<"                                                                  "<<std::endl;
  std::cerr<<"                                     Tools for Adaptive mesh refinement and Zcurves          "<<std::endl;
  std::cerr<<"                                         by me                       "<<std::endl
	     << std::endl << "INIT " << std::endl;
  }
}


/// @brief Destructor
///
/// Finalize MPI process
Node::~Node() {

  if (time  !=nullptr) delete time;
  if (scheme!=nullptr) delete scheme;
  if (m_inputOutput!=nullptr) delete m_inputOutput;

  loadBalancer.destroy();

  MPI__Free_types();
  MPI_Finalize();

}


/// @brief Initialize MPI process and get the nodes
/// @param [in] argc_ Program's arguments number.
/// @param [in] argv_ Program's arguments.
/// @return The node of each process.
Node* createNode(int* argc_, char*** argv_) {

#if __use_orchestrator
  int provided;
  MPI_Init_thread(argc_, argv_, MPI_THREAD_MULTIPLE, &provided);
  if ( provided != MPI_THREAD_MULTIPLE ) // we need more
    {
      std::cerr << "ERROR: MPI_THREAD_MULTIPLE unsupported!" << std::endl;
      exit(-1);
    }
#else 
  MPI_Init(argc_, argv_);
#endif /* __use_orchestrator */

  MPI__Init_types();

  LBS::init(*argc_, *argv_);

  return new NodeSingleDomain(MPI_COMM_WORLD);

}


/// @brief Get rank from comm manager
///
///
int Node::rank() {
  return commManager.getRank();
}


/// @brief Get number of nodes from comm manager
///
///
int Node::numberOfNodes() {
  return commManager.getNumberOfNodes();
}


/// @brief Check if there are accelerators
/// @return True if there is at least one GPU or Mic
bool Node::hasAccelerator() {
  return setup.numberOfGPUs + setup.numberOfMics > 0;
}


/// @brief Accessor comm manager
CommManager& Node::getManager() {
  return commManager;
}


/// @brief Configure parallelization and load balancing on the node
/// @param [in] configuration The node configuration
void Node::configure(Configuration<Node>& configuration) {
  
  switch (configuration.type) {

  case Configuration<Node>::FLAT :
    setup.numberOfThreads = configuration.maxNumberOfThreads;
    setup.numberOfGPUs = 0;
    setup.numberOfMics = 0;
    break;

  }

  loadBalancer.setParameters();
  loadBalancer.setMethod(configuration.balancingMethod);

  setup.logWorkload      = configuration.logWorkload;
  setup.expandFreeLimits = configuration.expandFreeLimits;
  setup.expandFactor     = configuration.expandFactor;

}


/// @brief Create integration scheme for the node from a configuration
/// @param [in] configuration The scheme configuration
void Node::createScheme(Configuration<NumericalScheme>& configuration) {
  scheme = ::createScheme(configuration);
}


/// @brief Create time manager for the node from a configuration
/// @param [in] configuration The time manager configuration
void Node::createTime(Configuration<TimeManager>& configuration) {
  time = new TimeManager(configuration);
}

/// @brief Create I/O manager for the node from a configuration (not used)
/// @param [in] configuration The I/O manager configuration
void Node::createInputOutput(Configuration<InputOutputManager>& configuration) {
  m_inputOutput = new InputOutputManager(&commManager,configuration);
}

/// @brief Link to an existing I/O manager
/// @param [in,out] _ioMgr Existing I/O manager to use
void Node::setInputOutput(InputOutputManager* _ioMgr) {
  m_inputOutput = _ioMgr;
}

#if __use_orchestrator
/// @brief Set shared variables
/// @param [in] _shared Initialized shared variables
void Node::setSharedVariables(SharedVariables& _shared)
{
  m_shared = &_shared;
}
#endif /* __use_orchestrator */

/// @brief Destructor
///
///
NodeSingleDomain::~NodeSingleDomain() {
  if (domain!=nullptr) delete domain;
}

/// @brief Create domain from domain configuration and particles configuration
/// @param [in] configuration Configuration of domain
/// @param [in] particles Configuration of particles (if no molecules)
/// @param [in] molecules Configuration of the molecules (if there is some)
void NodeSingleDomain::createDomains(Configuration<DomainInterface>& configuration, Configuration<Particle>& particles, Configuration<MPI__InMol>& molecules) {

  setup.expandFreeLimits = setup.expandFreeLimits && (!Global::domainInfo.noFreeBoundaryConditions());

  // ===========================================================================
  if (isMaster) std::cout<< "  Building domains" << std::flush ;

#ifdef __xsp_mem_info
  if (isMaster) std::cout<< std::endl;
#endif

  // ===========================================================================

  std::vector< std::vector<int> > domainsPerNode(numberOfNodes());

  for (uint i=0; i<domainsPerNode.size(); ++i) 
    domainsPerNode[i] = std::vector<int>(1, i);

  commManager.setDomainToNode(domainsPerNode);

  domain = buildDomain    (this, this->rank(), configuration, -1);

  // If ANY is not used, do not Balance ! Else, set callback functions
  if (configuration.decoupage!=Configuration<DomainInterface>::ANY) {
    loadBalancer.setMethod(LBS::NONE);
  }
  else {
    loadBalancer.setCallbackQueryFunctions(this->getNumberOfDomains(), this->getDomains());
  }

  // ===========================================================================
  commManager.barrier();
  if (isMaster) std::cout<< "\r  Building domains ................. ok" << std::endl;
  // ===========================================================================

  bool dumpStart=false;
  uint64_t slice=0, size=0;

  switch (particles.getInitType()) {
  case InitType::INIT_DEFAULT :
    dumpStart = false;
    size      = particles.totalNumberOfLatticeCells();
    slice     = __xsp_init_slice/particles.numberOfAtomPerCell();
    break;
  case InitType::INIT_STAMPV4 :
  	dumpStart = false;
  	size			= molecules.m_particles.size();
  	break;
  case InitType::INIT_STAMP_LEGACY_DUMP :
	  //std::cout << "PreInit createDomains legacy ... " <<  std::endl;
    dumpStart = true;
    size      = particles.totalNumberOfParticles();
    slice     = __xsp_init_slice;
    m_inputOutput->openIOFile(particles.getFilename(), "r", particles.getFileId());
    break;
  case InitType::INIT_HERCULE_DUMP :
	  //std::cout << "PreInit createDomains hercule ... " <<  std::endl;
    dumpStart = true;
    size      = particles.totalNumberOfParticles();
    slice     = size;
    //if (isMaster)  std::cout << __FILE__ << " INIT_HERCULE_DUMP " << std::endl;
   // commManager.openIOFile(particles.getFilename(), "r", particles.getFileId());
    break;
  }

  // ===========================================================================
  if (isMaster) std::cout<< "  Initializing particles " << std::flush;

#ifdef __xsp_mem_info
  if (isMaster && dumpStart) std::cout<< std::endl;
#endif

  // ===========================================================================

  switch (particles.getInitType()) {
  case InitType::INIT_DEFAULT :
  case InitType::INIT_STAMP_LEGACY_DUMP :
  {

	  uint64_t begin=0, end=0;
	  while (end<size) {

		  begin = end;
		  end   = auxMin(end+slice, size);
		  domain->initParticles(particles, begin, end);
#ifdef __xsp_init_barrier
		  commManager.barrier();
#endif

		  // ===========================================================================
		  if (isMaster) {
			  std::cout<< "\r  Initializing particles ........... "
			 << std::setw(6) << std::fixed << std::setprecision(2) << 100. * (double)end / (double)size << " % "
			 << std::flush;
		  }
		  // ===========================================================================

	  }
	  break;
  }
  case InitType::INIT_STAMPV4 :
  {
		domain->initParticles(molecules.m_particles.data(), molecules.m_initialTemperature, size, molecules.m_nCells, molecules.m_cellSize, molecules.m_numMol);
	  break;
  }
	  case InitType::INIT_HERCULE_DUMP :
	  {
#ifdef __use_lib_hercule
		  int initStep=particles.getInitStep();
		  if (initStep == -2){
			  vector<int> steps = this->m_inputOutput->herculeDumpR_steps();
			  if (steps.size()== 0){
				  std::cout << "Error no hercule dump base" << std::endl;
				  abort();
			  }
			  initStep=steps[steps.size()-1];
		  }
		  int nbSSDom=m_inputOutput->herculeDumpR_nbSSDom(initStep);

		  for (int iSSDom=0;iSSDom<nbSSDom;iSSDom++) {
			  if (isMaster){
				  //std::cout << "" << commManager.getRank() << ":herculeDumpR_openStep " << iSSDom << " / " << nbSSDom << "\n";
			  }
				this->m_inputOutput->herculeDumpR_openStep(initStep,iSSDom);
		  }
		  int modeInitHercule=this->m_inputOutput->herculeDumpR_modeInitHercule();
		  if (isMaster){
			  std::cout << "mode init hercule : " << modeInitHercule << std::endl;
		  }
		  int nbSSDomPerPhase=1;
		  switch(modeInitHercule){
		  case 0:
		  case 2:
			  // prot 0 lit et diffuse par groupe de 16
			  nbSSDomPerPhase=16;
			  break;
		  default:
		  case 1:
			  // mode proc 0 lit
			  nbSSDomPerPhase=1;
			  break;
		  }
		  for (int iSSDom=0;iSSDom<nbSSDom;iSSDom+=nbSSDomPerPhase) {
				  if (isMaster){
					  //std::cout << "" << commManager.getRank() << ":initParticles " << iSSDom << " / " << nbSSDom << "\n";
				  }
				  domain->initParticles(particles, iSSDom, HER_MIN(iSSDom+nbSSDomPerPhase,nbSSDom));

	//#ifdef __xsp_init_barrier
//				  commManager.barrier();
	//#endif

				  // ===========================================================================
				  if (isMaster) {
					  std::cout<< "\n  Initializing particles " << initStep << " .. " << iSSDom << " ........... "
					 << std::setw(6) << std::fixed << std::setprecision(2) << 100. * ((double)(HER_MIN(iSSDom+nbSSDomPerPhase,nbSSDom)) / (double)nbSSDom) << " % "
					 << std::flush;
				  }
				  // ===========================================================================

		  }

		  for (int iSSDom=0;iSSDom<nbSSDom;iSSDom++) {
				this->m_inputOutput->herculeDumpR_closeStep(iSSDom);
		  }
#else
#endif
		  break;
	  }
  }
  switch (particles.getInitType()) {
  case InitType::INIT_DEFAULT :
    break;
  case InitType::INIT_STAMP_LEGACY_DUMP :
    m_inputOutput->closeIOFile(particles.getFileId());
    break;
  case InitType::INIT_HERCULE_DUMP :
#ifdef __use_lib_hercule
	    //if (isMaster)  std::cout << __FILE__ << " INIT_HERCULE_DUMP " << std::endl;
    //commManager->closeIOFile(particles.getFileId());
#endif
    break;
  case InitType::INIT_STAMPV4 :
	  while(molecules.m_particles.size()!=0) {
	  	delete molecules.m_particles[molecules.m_particles.size()-1];
	  	molecules.m_particles.pop_back();
	  }
	  break;
  }
  
#ifdef __xsp_init_barrier
  commManager.barrier();
#endif



  // ===========================================================================
  if (isMaster) std::cout<< "\r  Initializing particles ........... ok      " << std::endl;
  // ===========================================================================

  domain->initDevice();

  // ===========================================================================
  if (isMaster) std::cout<< "  Initializing forces" << std::flush;
  // ===========================================================================

  // Compute forces to have everything initialized
  domain->initForces(particles.getInitType()==InitType::INIT_DEFAULT, particles.getInitTint());

  if (isMaster) {
    std::cout<< "\r  Initializing forces .............. ok" << std::endl;
    std::cout << std::endl;
  }

  if (isMaster && dumpStart) {
    std::cout<< "  Init from dump file " << particles.getFilename() << std::endl;
    std::cout << std::endl;
  }

}


void NodeSingleDomain::aferrgreg(uint64_t& send, uint64_t & recv) {

  domain->aferrgreg(send,recv);
}



/// @brief Reduce energies on all nodes
/// @return Vector with kinetic energy, potential energy, internal energy, chemical energy and kinetic energy in the center of momentum frame
Array<double> NodeSingleDomain::reduceEnergy() {

  Array<double> totalValues(8, 0.);
  Array<double> localValues(8, 0.);

  auto totalMomentum = domain->getTotalMomentum();

  localValues[0] = domain->getKineticEnergy();
  localValues[1] = domain->getPotentialEnergy();
  localValues[2] = domain->getInternalEnergy();
  localValues[3] = domain->getChemicalEnergy();
  localValues[4] = totalMomentum.x;
  localValues[5] = totalMomentum.y;
  localValues[6] = totalMomentum.z;
  localValues[7] = domain->getTotalMass();

  commManager.allReduce(localValues, totalValues, MPI_SUM);

  Array<double> localKinShift(1, 0.);
  Array<double> totalKinShift(1, 0.);

  localKinShift[0] = domain->getShiftedKineticEnergy( vec3<double>(totalValues[4], totalValues[5], totalValues[6]) / totalValues[7] );

  commManager.reduce(localKinShift, totalKinShift, MPI_SUM);

  // [0] is eKin
  // [1] is ePot
  // [2] is eInt
  // [3] is eChm
  // [4] is eKin in the center of momentum frame

  Array<double> tmp(5, 0.);
  tmp[0] = totalValues[0];
  tmp[1] = totalValues[1];
  tmp[2] = totalValues[2];
  tmp[3] = totalValues[3];
  tmp[4] = totalKinShift[0];

  return tmp;
}


/// @brief Reduce number of particles on all nodes
/// @return Number of particles
uint64_t NodeSingleDomain::reduceNumberOfParticles() {

  uint64_t n = domain->getNumberOfParticles();

  Array<uint64_t> totalNumberOfParticles(1, 0);
  Array<uint64_t> localNumberOfParticles(1, n);
  commManager.reduce(localNumberOfParticles, totalNumberOfParticles, MPI_SUM);

  return totalNumberOfParticles[0];

}


/// @brief Reduce number of cells on all nodes
/// @return Number of particles
uint64_t NodeSingleDomain::reduceNumberOfCells() {

  uint64_t n = product(Global::domainInfo.getNumberOfCellsPerDim());

  return n;

}


/// @brief Reduce pressure on all nodes
/// @return Pressure
double NodeSingleDomain::reducePressure() {

  Array<double> totalValues(4, 0.);
  Array<double> localValues(4, 0.);

  auto totalMomentum = domain->getTotalMomentum();

  localValues[0] = totalMomentum.x;
  localValues[1] = totalMomentum.y;
  localValues[2] = totalMomentum.z;
  localValues[3] = domain->getTotalMass();
  
  commManager.allReduce(localValues, totalValues, MPI_SUM);

  double localPressure;
  double totalPressure;

  localPressure = dot(domain->getPressure( vec3<double>(totalValues[0], totalValues[1], totalValues[2]) / totalValues[3] ),diag(1.));

  commManager.reduce(localPressure, totalPressure, MPI_SUM);

  return totalPressure;

}

/// @brief Calculate imbalance rate
/// @param [in] write Print imbalance rate if true, default=false
/// @param [in] step Step where imbalance rate is printed, default=0
/// @return Imbalance rate
double NodeSingleDomain::imbalanceRate(bool write, int step) {

  static bool firstCall = true;
  static std::ofstream out;

  if (firstCall) { 
    firstCall = false;
    if (isMaster && setup.logWorkload) {
      std::string filename = "xstamp-workload.log";
      out.open(filename.c_str(), std::fstream::out | std::fstream::trunc);
      out << "# " << std::endl
  	  << "# Workload at each step : step, avg, std, min, max, imb, rate " << std::endl 
  	  << "# " << std::endl;
    }
  }

  double w = domain->workload();

  Array<double> totalWorkload(commManager.getNumberOfNodes(), 0);
  Array<double> localWorkload(1, w);

  commManager.allGather(localWorkload, totalWorkload);

  double mean, std, min, max;
  quickAnalysis(totalWorkload.data(), totalWorkload.size(), mean, std, min, max);

  if (write && isMaster && setup.logWorkload) {
    out << std::fixed << std::setw(8) << step << "   " 
	<< std::scientific << std::setprecision(3) << mean << " " << std << " " << min << " " << max << "   " << max/mean-1.0 
	<< std::endl;
  }

  return max/mean-1.0;

}
