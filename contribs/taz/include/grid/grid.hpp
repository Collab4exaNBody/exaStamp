/// @file 
/// @brief Definition of class Grid

#ifndef __GRID_HPP_INCLUDED
#define __GRID_HPP_INCLUDED


#include <algorithm>
#include <ostream>
#include <utility>
#include <vector>

#include "globals.hpp"

#include "grid/moduler.hpp"
#include "grid/anyGridInfo.hpp"
#include "grid/rectilinearGridInfo.hpp"


template <class Grid_impl> class Domain;
template <class P> class MessageSend;


/// @brief Shortcut for a template that depend on the grid implementation
#define TMPLG template <class Grid_impl>
/// @brief Shortcut for a grid
#define TMPL_Grid Grid<Grid_impl>


/// @brief Class Grid : handle the point particles
/// @tparam Grid_impl The proper Grid subclass
TMPLG class Grid {

public:

  Grid();

  /// @brief Destructor (nothing to do)
  virtual ~Grid();

  // Movement methods
  void pushPositions1stOrder (double time);
  void pushPositions2ndOrder (double time);
  void pushVelocities1stOrder(double time);
  void pushDissipationLangevin(double time, double gamma, double beta);
  void pushFluctuationDPD(double time);

  void resetForce();
  void resetForceFluct();
  void resetEnergy();
  void resetForceAndEnergy();
  void resetDensity();
  void resetChemistry();
  
  void makeNeighborLists();
  void makeNeighborListsInside();
  void makeNeighborListsOnEdges();

  uint64_t getTotalNbOfParticlesSend();
  uint64_t getTotalNbOfParticlesRecv();

  void clearNeighborList();

  void sortNeighbors(::Traversal T);

  void computeForcesVerlet(::Traversal T);

  void exchangeParticles();
  void collectParticles();

  void refineCells();

  void clearGhost();
  void updateGhost();
  void collectGhost();

  void internalReorganization();
  void checkLocalBuffers();

  uint64_t getNumberOfParticles();
  uint64_t getNumberOfRealCells();
  
  double getTotalEnergy();
  double getKineticEnergy();
  double getPotentialEnergy();
  double getInternalEnergy();
  double getChemicalEnergy();

  double workload();

  vec3<double> getTotalMomentum();
  double getTotalMass();
  double getShiftedKineticEnergy(const vec3<double>& vshift);

  mat3<double> getPressure(const vec3<double>& vshift);
  
  void initParticles(Configuration<Particle>& particles, uint64_t cBegin, uint64_t cEnd);
  void initParticles(MPI__Particle** particles, double initialT, uint64_t numP, vec3<int> nCells, vec3<double> cellSize, uint64_t numMol);

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
  
  void updateOwners(Array<int>& allCellsOwners, int numExport, uint* exportLocalGids, int* exportProcs, int numImport, uint* importGlobalGids);
  void exchangeMovingCellEnv(const Array<uint>& procsToComm, const Array<int>& allCellsOwners, int numImport, std::vector<uint32_t>& importEnv, uint* importGlobalGids, int numExport, uint* exportGlobalGids, uint* exportLocalGids, int* exportProcs);

  bool checkVerlet();
  void refineBalance();

  void ctest(std::set<ParticleReferenceValue> &reference_values_set, double &ae, double &re) {impl().ctest(reference_values_set, ae, re);};

  template <class Q>
  void dumpParticles(std::vector<Q>& toKeep, MessageSend<Q>& toSend, const uint* sendIndexes, const int* sendProcs, const uint sendSize) { 
    impl().dumpParticles(toKeep, toSend, sendIndexes, sendProcs, sendSize);
  }
  
  void dumpParticles(keepAtom& toKeep,exchangeGhost &balanceComm,  const uint* sendIndexes, const int* sendProcs, const uint sendSize) { 
    impl().dumpParticles(toKeep, balanceComm, sendIndexes, sendProcs, sendSize);
  }

  void printInfo(std::ostream& flux);
  void print(std::ostream& flux);

  void updateGhostVel();
  void updateGhostFluct();
  void updateDensities();
  
  inline void destructorLeafCells() {
    impl().destructorLeafCells();
  }
 
protected:


  void build(Domain<Grid_impl>* domain_, AnyGridInfo* info);
  void build(Domain<Grid_impl>* domain_, RectilinearGridInfo* info);
  // void build(Domain<Grid_impl>* domain_, VoronoiGridInfo* info) {}

  void setGhostData(AnyGridInfo* info);
  void setGhostData(RectilinearGridInfo* info);
  // void setGhostData(VoronoiGridInfo* info) {}

  void makeTraversals(AnyGridInfo* info, const std::vector< vec3<int> >& coords);
  void makeTraversals(RectilinearGridInfo* info, const std::vector< vec3<int> >& coords);
  // void makeTraversals(VoronoiGridInfo* info, const std::vector< vec3<int> >& coords) {}


  uint getDomainIndex();

  const Array<uint>& getTraversal(TraversalManager::CellTraversal traversal);
  const Array<uint>& getTraversal(::Traversal T, uint level=0);

  void correct(vec3<double>& distance);
  void correct(uint size, double* rx, double* ry, double* rz);

  void module(vec3<double>& distance);
  void module(vec3<int>& distance);

  void module(double& x, double& y, double& z);

  void lock(uint i);
  void unlock(uint i);

  bool isGhost(uint i);
  bool isReal(uint i);

  /// @brief Check if specified cell is on system min boundaries
  /// @param [in] index Index of the cell
  /// @param [in] dim Dimension to check
  /// @return True if cell is on min boundary
  bool onMinBounds(const uint index, const uint8_t dim) { return coords[index][dim]==0; }
  /// @brief Check if specified cell is on system max boundaries
  /// @param [in] index Index of the cell
  /// @param [in] dim Dimension to check
  /// @return True if cell is on max boundary
  bool onMaxBounds(const uint index, const uint8_t dim) { return (coords[index][dim]+1)==Global::domainInfo.getNumberOfCellsPerDim()[dim]; }

  uint64_t numberOfParticles;       ///< Number of particles (real only)
  uint64_t totalNumberOfParticles;  ///< Total number of particles (real and ghost)

  bool energyFlag;   ///< Flag to know if energies are up to date
  double energyKin;  ///< Kinetic energy of particles
  double energyPot;  ///< Potential energy of particles
  double energyInt;  ///< Internal energy of particles
  double energyChm;  ///< Chemical energy of particles

  bool m_workloadFlag; ///< Flag to know if workload has been computed yet
  double m_workload; ///< Workload for the grid/domain

  Correcter correcter;   ///< Tool to correct distances
  Moduler<double> modD;  ///< Tool to modulate double distances
  Moduler<int> modI;     ///< Tool to modulate integer distances

  TraversalManager tManager;  ///< Traversal manager

  Domain<Grid_impl>* domain;  ///< Pointer to the domain

  MMutex mtx;                      ///< A set of mutex to protect cells

  std::vector< vec3<int>  > coords;      ///< Cells coordinates
  std::vector< Array<int> > neighbors;   ///< Cells neighbor indexes
  std::vector< uint8_t    > ghostLayer;  ///< Ghost layer for each cell : 0 for real cells, ghost thickness if last layer of ghost and ghost thickness -1 if ghost but not last layer

  std::vector< std::vector<int> >         destDomains; ///< Recipient domains during ghost exchange for each cell
  
  std::vector< std::vector<vec3<int>> >   raphael_test; ///< Recipient domains during ghost exchange for each cell
  
  std::vector< std::vector< vec3<int> > > destCells; ///< Recipient cells during ghost exchange for each cell
  


private:

  Grid_impl& impl();

};


#include "grid/grid.hxx"
#include "grid/grid_any.hxx"
#include "grid/grid_rectilinear.hxx"


#undef TMPLG
#undef TMPL_Grid

#endif // __GRID_HPP_INCLUDED
