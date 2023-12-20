/// @file 
/// @brief Definition of class AMRGrid

#ifndef __AMR_GRID_HPP_INCLUDED
#define __AMR_GRID_HPP_INCLUDED



#include "grid/correcter.hpp"
#include "cellList/AMR/cellListParticleAMR.hpp"
#include "grid/grid.hpp"
#include "grid/AMRGrid/AMRGrid_coloring.hxx"
#include "io/particleInput.hpp"
#include <utility>
#include "parallel/communications/messageCenter.hpp"
#include "parallel/thread/thread.hpp"
#include "parallel/types/MPI_cell.hpp"
#include "particle/lattice.hpp"
#include "utils/kernelFunction/allKernelFunctions.hpp"
#include "parallel/commAMR/mycommAMR.hpp"

#include "grid/AMRGrid/externAMRGrid.hpp"


/// @brief Shortcut for a template that depend on an info class, a type of particles with its ghost class, some unused integer and an alignment size
#define TMPLSG template <class Info, class P>
/// @brief Shortcut for a grid with one class of particles
#define TMPL_AMRGrid AMRGrid<Info, P>

// It should be in the SingleSpecGrid class but it does not compile ...
/// @brief An affinity partitioner for the communications
static affinity_hint comm_partitioner;
/// @brief An affinity partitioner for the force computations
static affinity_hint force_partitioner;
/// @brief An affinity partitioner for the neighbors search
static affinity_hint nbr_partitioner;

/// @brief A Grid subclass for the case when there is only one class of particles
/// @tparam Info Grid info type
/// @tparam P Type of particles
/// @tparam elem_chunk Base chunk for cell arrays
/// @tparam align Alignment size for vectorization
TMPLSG class AMRGrid : public Grid<TMPL_AMRGrid> {

  /// @brief Shortcut for the grid type
  typedef TMPL_AMRGrid self_type;
  /// @brief Shortcut for the type of a pointer to a grid type
  typedef self_type* self_type_ptr;
  /// @brief Shortcut for the type of a reference to a grid type
  typedef self_type& self_type_ref;

  static uint thread_grain; ///< Control the number of chunk in the parallel regions that run on cells

public:

  /// @brief Shortcut for the exchange type (to allow use of that type by the domain)
  typedef P exchange_t;
  /// @brief Shortcut for the cellList type (to allow use of that type by the primary Grid)
  typedef  Octree  cellList_t;
  
  typedef leafCell Cell;

  AMRGrid(Domain<self_type>* domain_, uint index, uint64_t numberOfParticles);

  ~AMRGrid();
  
  void pushPositions1stOrder (double time);
  void pushPositions2ndOrder (double time);
  void pushVelocities1stOrder(double time);

  void pushDissipationLangevin(double time, double gamma, double beta);
  void pushFluctuationDPD(double time);

  void resetForce();
  void resetEnergy();
  void resetForceAndEnergy();

  void makeNeighborLists();
  void makeNeighborListsInside();
  void makeNeighborListsOnEdges();
  void clearNeighborList();
  void sortNeighbors(::Traversal traversal);

  void computeForceVerlet(PairPotential* potential, uint8_t typeIndexA, uint8_t typeIndexB, ::Traversal T);
  void computeForceVerlet(EAMPotential* potential, uint8_t typeIndexA, uint8_t typeIndexB, ::Traversal T);

  void fillWithLeavingParticles();
  void sendLeavingParticles();
  void collectParticles();

  void clearGhost();
  void fillWithGhosts();
  void sendGhosts();
  void collectGhost();

  void internalReorganization();
  void checkLocalBuffers();
  void checkLocalBuffersGHOST();

  void initParticles(Configuration<Particle>& particleInit, uint64_t cell_begin, uint64_t cell_end);

  template <class Q, typename F>
  void addParticles(const Q* p, const uint size, const bool updateCount, F f);
  
  void unPackBalanceKeep(keepAtom &balanceComm);
  void unPackBalance(exchangeGhost &balanceComm);

  void computeEnergy();

  void computeWorkload();

  vec3<double> computeTotalMomentum();
  double computeTotalMass();
  double computeShiftedKineticEnergy(const vec3<double>& vshift);

  mat3<double> computePressure(const vec3<double>& vshift);

  cellList_t* getCells();

  void applyWallCondition(const int dim, const int direction);

  void checkFreeBoundaries(Array<int>& outOfFreeBounds);

  void stopWalls();
  
  void fillBuffer(ParticleOutput* buffer);
  template <class DumpStruct>
  void fillBuffer(Array<DumpStruct>& buffer);
  template <class DumpStruct>
  void fillBuffer(DumpStruct& buffer);

  void fillBuffer(ParticleInSitu* buffer);
  void fillGhostBuffer(ParticleInSitu* buffer);

  void fillBuffer(CellOutput* buffer);
  
  void print(std::ostream& flux);

  bool checkVerlet();


#if __use_lib_zoltan

  static int numberOfCells(void *data, int *ierr);
  static void listObj (void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr);
  static int numberOfDim (void *data, int *ierr);
  static void fillCoords (void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr);
  static int numberOfEdges(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr);
  static void fillEdges (void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs, int wgt_dim, float *ewgts, int *ierr);

#endif

  void updateOwners(Array<int>& allCellsOwners, int numExport, uint* exportLocalGids, int* exportProcs, int numImport, uint* importGlobalGids);
  void exchangeMovingCellEnv(const Array<uint>& procsToComm, const Array<int>& allCellsOwners, int numImport, std::vector<uint32_t>& importEnv, uint* importGlobalGids, int numExport, uint* exportGlobalGids, uint* exportLocalGids, int* exportProcs);
  void dumpParticles(std::vector<P>& toKeep, MessageSend<P>& toSend, const uint* sendIndexes, const int* sendProcs, const uint sendSize);
  void dumpParticles(keepAtom& toKeep, exchangeGhost &balanceComm, const uint* sendIndexes, const int* sendProcs, const uint sendSize);
  void updateDecomposition(const Array<int>& allCellsOwners);


protected:

  Info info; ///< Grid info : defines grid type and build parameters
  InputOutputManager* m_inputOutputManager; ///< Input/Output manager
  MessageCenter<P> particleExchange; ///< Message center to exchange particles
  MessageCenter<MPI__Ghost<typename P::Base> > ghostUpdate; ///< Message center to exchange ghost
  MessageCenter<ExchangeEAM> exchangeEAM; ///< Message center to exchange EAM density
  MessageCenter<ExchangeEAM> exchangeDst; ///< Message center to exchange density
  MessageCenter<ExchangeGhostV3> exchangeGhostVel; ///< Message center to exchange ghost velocities
  MessageCenter<ExchangeGhostV3> exchangeGhostFluct; ///< Message center to exchange ghost fluctuations
  
  
// TRY
  mapInfoOctree  m_mapForGhostOctree;
  vector2DofPair m_idxLeafsInOctree;
  exchangeGhost ghostComm;
  exchangeGhost embComm;
  uint64_t m_nbOfParticlesSend;
  uint64_t m_nbOfParticlesRecv;

  std::vector< Octree > particles; ///< Cells containing the particles
  const uint levelMax;                                         ///< refinement maximum
  int doRefine;
  int refineStepMax;  
  int iteration;                                        ///< number of steps simulation
  vec3<int> numberOfCellsPerDim;
  uint g_amrCriterionValue;
  int *depGraph;
  coloringOctree cWave;
  coloringOctree cWaveEdge;

  
  // TraversalAMR
  // These arrays should be in the traversal (include in gridUtils)
  std::vector< Octree* >  octrees; // root cell real
  std::vector< size_t >  octreesMEAM; // root cell real + ghost
  std::vector< leafCell* >  leafcells;    //leaf cell real
  std::vector< leafCell* >  leafcellsMEAM; //leaf cell real + ghost
  std::vector < std::set < size_t> > recDirectionInOctree;

public : 


  /* AMRGrid_graph.hxx */
  template<typename F>void launchGraph(F potential);
  template<typename F>void launchGraphEdge(F potential);
  template<typename F>void OpenMPTask(std::size_t nbDep, std::size_t uTask, int idOctree, std::size_t* pTask, F potential);
  bool isTask(vec3<int> p);


  template<class Pot_t> void AllocationNode_computeForcePair(Pot_t * pot, uint8_t typeIndexA, uint8_t typeIndexB );
  void BuildWaveMethod();

  void refineCells();
  void Refine();
  void destructorLeafCells();

  uint64_t getTotalNbOfParticlesSend();
  uint64_t getTotalNbOfParticlesRecv();

  /* Ctest */ 
  void ctest(std::set<ParticleReferenceValue> &reference_values_set, double &ae, double &re);

private:

  void initVelocities(P& p, Array<double>& sigma);
  void initParticlesHard(Configuration<Particle>& particleInit, uint64_t cell_begin, uint64_t cell_end);
  void initParticlesFile(Configuration<Particle>& particleInit, uint64_t part_begin, uint64_t part_end);
#ifdef __use_lib_hercule
  void initHerculeParticlesFile(Configuration<Particle>& particleInit, uint64_t part_begin, uint64_t part_end);
#endif
  void __makeNeighborLists_gt1(TraversalManager::CellTraversal t);

  template<typename F> void __RecApply(F function, Cell *C);
  void __makeNeighborLists_Traversal_AMR(TraversalManager::CellTraversal t);
  void __makeNeighborLists_Traversal_MEAM_AMR(TraversalManager::CellTraversal t);
  
  template <class Pot_t>
  void __computeForcePairVerlet(Pot_t* pot, ::Traversal T, const uint8_t typeIndexA, const uint8_t typeIndexB);
  
    template <class Pot_t>
  void __computeForceEAMVerlet (Pot_t* pot, ::Traversal T, const uint8_t typeIndexA, const uint8_t typeIndexB);
  
    template <class Pot_t>
  void __computeForceMEAMVerlet(Pot_t* pot, ::Traversal T, const uint8_t typeIndexA, const uint8_t typeIndexB);

  void __resetEAMData(const Array<uint>& cells);

  template <class Pot_t>
  void __computeForceEmb(Pot_t* pot, const Array<uint>& cells, const uint8_t typeIndexA, const uint8_t typeIndexB);

  template <class Pot_t>
  void __computeForcesAndPotentialMEAMVerlet(Pot_t* pot, const Array<uint>& cells, const uint8_t typeIndexA, const uint8_t typeIndexB);
  
  
  // AMR functions //
  vec3<int> get_Vec3_Val(Cell *C);
  bool testInf(vec3<int> a, vec3<int> b);
  bool testSup1(vec3<int> a);
  void fillneighborCells_List();
  void fillNeighborsCells_list(Cell *C);
  void fillNeighborsCells_list_GHOST(Cell *C);
  void CellSearch(uint32_t neighborRoot, Octree *root, Cell *C);
  void recursiveCellSearch_level_0(Cell *root, Cell *C);
  void recursiveCellSearch(Cell *root, Cell *C);
  void getDifferentLevelNeighborsCellInOtherNode(Cell *root, Cell *C);
  void cellSearch_1(int noNeighbor, Cell *C);
  void recursiveCellSeach_1(Cell *root, Cell *C);
  void buildInfoAMR();

  void __makeNeighborLists_gt1_AMR();

  template <bool perBlock> void __makeNeighborLists_gt1_AMR(Cell* C);
  template <bool perBlock> void __makeNeighborLists_gt1_AMR(size_t indexOctree, size_t indexLeafCell);

  void __makeNeighborLists_MEAM_AMR(Cell* C);
  void UpdateTraversalAMR();
  void UpdateTraversalAMRoctrees();
  void UpdateTraversalAMR_MEAM();
  void Adjust();
  void fillRecWithGhost_cell(uint destCellIndex, std::set<size_t>& tab, Cell *cell);
  void fillRecWithGhost_cell(uint destDomain, vec3<int>& desCell, std::set<size_t>& tab, Cell *cell, std::vector< std::tuple<int, MPI__Ghost<typename P::Base> > >& leaving);

  void checkParticlesPosition();
  // end //

  void fillEmbSelf();
  void sendEmb();
  void collectEmb();


  void fillWithLeaving_cell(uint cellIndex, std::list< std::tuple<int, P> >& staying, std::list< std::tuple<int, P> >& leaving);
  void fillWithGhost_cell (uint cellIndex, std::vector< std::tuple<int, MPI__Ghost<typename P::Base> > >& leaving);

  void packGhost();
  void unPackGhost();
  void packEmb();
  void unPackEmb();
  void packBalance(exchangeGhost &balanceComm);
  void packBalanceKeep(keepAtom &balanceComm);  
  void fillWithGhostSelfDomain (uint cellIndex);
  void fillWithEmbSelfDomain(size_t cellIndex);
  void fillRecWithGhostSelfDomain(uint destCellIndex, std::set<size_t>& tab, Cell *cell);
  void fillRecWithEmbSelfDomain(size_t destCellIndex, std::set<size_t>& tab, Cell *cell);
  void fillWithGhostOtherDomain (uint cellIndex, std::vector< std::pair<size_t,infoOctreeGhost > >& leaving);
  void fillRecWithGhostOtherDomain(std::pair<size_t,infoOctreeGhost> &preFill, std::set<size_t>& tab, Cell *cell, std::vector< std::pair<size_t,infoOctreeGhost> >& leaving);

};


TMPLSG uint TMPL_AMRGrid::thread_grain = 0;


/// @brief Constructor from a domain, its index and a grid info type
/// @param [in,out] domain_ Domain
/// @param [in] index Index of the domain
/// @param [in] not_used Why is that there ?
TMPLSG TMPL_AMRGrid::AMRGrid(Domain<self_type>* domain_, uint index, uint64_t not_used)
	: Grid<self_type>(),
	info(index),
	m_inputOutputManager(domain_->m_inputOutputManager),
	particleExchange(domain_->getCommManager(), ISession::PARTICLE_EXCHANGE),
	ghostUpdate (domain_->getCommManager(), ISession::GHOST_UPDATE),
	exchangeEAM (domain_->getCommManager(), ISession::EAM_EMB),
	exchangeDst (domain_->getCommManager(), ISession::UPDATE_DENSITIES),
	exchangeGhostVel (domain_->getCommManager(), ISession::GHOST_UPDATE_VELOCITIES),
	exchangeGhostFluct (domain_->getCommManager(), ISession::GHOST_UPDATE_FLUCTUATION),
	m_mapForGhostOctree() , m_idxLeafsInOctree(), ghostComm(domain_->getCommManager()->getCommunicator()), embComm(domain_->getCommManager()->getCommunicator()),
	m_nbOfParticlesSend(0), m_nbOfParticlesRecv(0), 
	particles(), levelMax(Global::reference.getDmax()), doRefine(0), iteration(0), g_amrCriterionValue(Global::reference.r_amrCriterion) 
{

	// Set size of the cell array
	particles.resize(info.getTotalNumberOfCells()); 

	numberOfCellsPerDim = info.getCellSup() - info.getCellInf();
	refineStepMax = 0;

	// Describe data in the buffer
	void* tmpPtr=nullptr;

	ghostComm.defineSizeOfData(
		(double*)   (tmpPtr),	//rx
		(double*)   (tmpPtr),	//ry
		(double*)   (tmpPtr),	//rz
		(uint64_t*) (tmpPtr),	//id
		(uint8_t*)  (tmpPtr)	//ty
	);

	embComm.defineSizeOfData(
		(double*)   (tmpPtr)	//emb
	);

	// Les messages MPI ne peuvent pas contenir plus d'1 million d'atomes
	ghostComm.defineSizeOfDataPerMsg(1000000);
	embComm.defineSizeOfDataPerMsg(1000000);
 
	defineRecDirectionInOctree(recDirectionInOctree);
      
	thread_grain = 2*Global::maxNumberOfThreads;

	// Build the grid
	this->build(domain_, &info);

	buildInfoAMR();

	// Set data for ghost exchange
	this->setGhostData(&info);

	auto tmpArray = info.getNeighborArray();
	Array<uint> tmp_neighbors (tmpArray.size()); 

	for(size_t i=0 ; i<tmpArray.size();i++)
		tmp_neighbors[i] = tmpArray[i];

	// Initialize the message centers
	particleExchange.init(this->getDomainIndex(), tmp_neighbors);
	ghostUpdate.init (this->getDomainIndex(), tmp_neighbors);

	if (Global::reference.isEAM())
		exchangeEAM.init(ghostUpdate);
}

TMPLSG inline void TMPL_AMRGrid::buildInfoAMR()
{
	const auto& cells = this->getTraversal(TraversalManager::ALL);

  #pragma omp parallel for schedule(static)
  for(int i = 0 ; i < cells.size() ; i++)
  {
		vec3<int> position = info.coords(cells[i])+1;//tmp_freeOrNot; 
		particles[cells[i]].setInfo(position,0);
		particles[cells[i]].setOffset(position, this->ghostLayer[cells[i]]);
  }
}

/// @brief Destructor
///
/// Clear the cell array
TMPLSG TMPL_AMRGrid::~AMRGrid() {
	particles.clear();
}

/// @brief Add a set of particles to the grid
/// @tparam Q The type of particle added
/// @tparam F A function type
/// @param [in] q Set of particles to add
/// @param [in] size Size of the set
/// @param [in] updateCount Indicates if the number of particles must be changed according to the addition
/// @param [in] f Function to identify the cell of each particle
TMPLSG template <class Q, typename F>
inline void TMPL_AMRGrid::addParticles(const Q* q, const uint size, const bool updateCount, F f) {

	std::vector<uint> indexes(size, -1);

	// Parallelize the work by distributing the particles between the threads
	parallel_region(0, size, [&](const uint begin, const uint end) {

	info.findCellIndex(q + begin, indexes.data() + begin, end-begin, f);

	for(uint i=begin; i<end; ++i) {


		const uint& cell = indexes[i];

		this->lock(cell);
		particles[cell].add(q[i]);
		this->unlock(cell);

	}
	});

	// increase number of particles
	if (updateCount) {
		this->lock(0);
		this->numberOfParticles += size;
		this->unlock(0);
	}

}


/// @brief Debug print for all cells in specified flux
/// @param [in,out] flux Print flux
TMPLSG void TMPL_AMRGrid::print(std::ostream& flux) {

  const auto& cells = this->getTraversal(TraversalManager::REAL);

  // basic print of real octrees
  for (uint i=0; i<cells.size(); ++i) {
    uint cell = cells[i];
    particles[cell].__debug_print(flux);
  }

}


/// @brief Accessor to the cells
/// @return Pointer to the cellLists
TMPLSG  Octree * TMPL_AMRGrid::getCells() {
	return particles.data();
}

TMPLSG void TMPL_AMRGrid::ctest(std::set<ParticleReferenceValue>  &reference_values_set, double &ae, double &re)
{
  const auto& cells = this->getTraversal(TraversalManager::REAL);

  auto fun =  [&reference_values_set, &ae, &re](
	uint64_t m_numberOfParticles, uint64_t* m_id, 
	double * m_rx, double * m_ry, double * m_rz,
	double * m_vx, double * m_vy, double * m_vz,
	double * m_fx, double * m_fy, double * m_fz)
    {
      for(uint64_t i=0; i<m_numberOfParticles; i++)
      {
        auto it = reference_values_set.find( ParticleReferenceValue{ m_id[i] } );
        //particleForce[ pa.m_start + i ].x = pa.m_fx[i];
        //particleForce[ pa.m_start + i ].y = pa.m_fy[i];
        //particleForce[ pa.m_start + i ].z = pa.m_fz[i];
        if( it != reference_values_set.end() )
        {
          double drx = m_rx[i] - it->m_r[0];
          double dry = m_ry[i] - it->m_r[1];
          double drz = m_rz[i] - it->m_r[2];
          re += drx*drx + dry*dry + drz*drz;
          
          double dax = m_fx[i] - it->m_a[0];
          double day = m_fy[i] - it->m_a[1];
          double daz = m_fz[i] - it->m_a[2];
          ae += dax*dax + day*day + daz*daz;          
        }
      }
    };


  // basic print of real octrees
  for (uint i=0; i<cells.size(); ++i) {
    uint cell = cells[i];
    auto& oct = particles[cell];
    fun(oct.size,oct.getId(),
	oct.getPositionX(),oct.getPositionY(),oct.getPositionZ(),
	oct.getVelocityX(),oct.getVelocityY(),oct.getVelocityZ(),
	oct.getForceX(),oct.getForceY(),oct.getForceZ());
  }	
}


#include "grid/AMRGrid/AMRGrid_resetParticlesFields.hxx"
#include "grid/AMRGrid/AMRGrid_wall.hxx"
#include "grid/AMRGrid/AMRGrid_buffer.hxx"
#include "grid/AMRGrid/AMRGrid_push.hxx"
#include "grid/AMRGrid/AMRGrid_reduction.hxx"
#include "grid/AMRGrid/AMRGrid_initParticles.hxx"
#include "grid/AMRGrid/AMRGrid_balance.hxx"
#include "grid/AMRGrid/AMRGrid_communications.hxx"
#include "grid/AMRGrid/AMRGrid_force_Verlet.hxx"
#include "grid/AMRGrid/AMRGrid_neighbors.hxx"
#include "grid/AMRGrid/AMRGrid_graph.hxx"
#include "grid/AMRGrid/AMRGrid_refinement.hxx"

#undef TMPLSG
#undef TMPL_AMRGrid

#endif // __AMR_GRID_HPP_INCLUDED
