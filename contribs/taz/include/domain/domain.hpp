/// @file 
/// @brief Definition of domain class and associated classes

#ifndef __DOMAIN_HPP_INCLUDED
#define __DOMAIN_HPP_INCLUDED


#include <algorithm>
#include <iomanip>
#include <iostream>

#include "referenceMap.hpp"

#include "domain/domainInfo.hpp"
#include "domain/domainInterface.hpp"

#include "grid/gridUtils.hpp"

#include "parallel/communications/session.hpp"
#include "parallel/node.hpp"

#include "parallel/thread/thread.hpp"
#include "parallel/commAMR/mycommAMR.hpp"
#include "parallel/commAMR/keepAtomClass.hpp"
/// @brief Base class to handle domain local data
class LocalInfo {

public:

  /// @brief Constructor
  /// @param [in] index_ Domain index
  /// @param [in] symmetrization_ Boolean allowing symmetrization, converted to an integer
  LocalInfo(uint index_, bool symmetrization_)
    : index(index_), 
      symmetrization(symmetrization_? 1:0) {}

  /// @brief Destructor (nothing to do)
  ~LocalInfo() {}

  uint getIndex();
  bool symmetrize();

  void print(std::ostream& flux);

protected:

  uint index; ///< Domain index
  uint symmetrization; ///< Allow symmetrization of the force calculation if 1

};


/// @brief Shortcut for a template that depend on the grid class
#define TMPLG template <class Grid_impl>
/// @brief Shortcut for a domain
#define TMPL_Domain Domain<Grid_impl>


class CommManager;
class InputOutputManager;
template <class T> class MessageCenter;
template <class Grid_impl> class Grid;

/// @brief Domain class, inherit externally called methods from DomainInterface and local members from LocalInfo
/// @tparam Grid_impl The proper Grid subclass
TMPLG class Domain : public DomainInterface, public LocalInfo {

public:


  Domain(Node* node, uint index, const Configuration<DomainInterface>& config, uint64_t numberOfParticles);

  virtual ~Domain();

  virtual void initDevice();

  virtual void initParticles(Configuration<MPI__Particle>& particles, uint64_t cBegin, uint64_t cEnd);
  virtual void initParticles(MPI__Particle** particles, double initialT, uint64_t numP, vec3<int> nCells=vec3<int>(1), vec3<double> cellSize=vec3<double>(0.), uint64_t numMol=0);

  virtual void pushPositions1stOrder (double time);
  virtual void pushPositions2ndOrder (double time);
  virtual void pushVelocities1stOrder(double time);

  virtual void pushDissipationLangevin(double time, double gamma, double beta);
  virtual void aferrgreg(uint64_t& send, uint64_t & recv);

  virtual void makeNeighborLists();
  virtual void makeNeighborListsInside();
  virtual void makeNeighborListsOnEdges();
  virtual void clearNeighborLists();

  virtual void initForces(bool doInitEint, double initialTint);
  virtual void ctest(std::set<ParticleReferenceValue> &reference_values_set, double &ae, double &re);  


  virtual void doComputeForcesVerlet();
  virtual void doComputeForcesInsideVerlet();
  virtual void doComputeForcesOnEdgesVerlet();

  virtual void refineCells();
  virtual void updateCells();
  virtual void updateCellsInside();
  virtual void clearGhost();
  virtual void updateGhost() ;
  virtual void collectGhost();

  virtual uint64_t getNumberOfParticles();
  virtual uint64_t getNumberOfRealCells();
  
  virtual double getTotalEnergy();
  virtual double getKineticEnergy();
  virtual double getPotentialEnergy();
  virtual double getInternalEnergy();
  virtual double getChemicalEnergy();

  virtual vec3<double> getTotalMomentum();
  virtual double getTotalMass();
  virtual double getShiftedKineticEnergy(const vec3<double>& vshift);

  virtual mat3<double> getPressure(const vec3<double>& vshift);

  virtual CommManager* getCommManager();
  virtual InputOutputManager* getInputOutputManager();

  virtual double workload();
  virtual LBS::State balance(LBS::LoadBalancer* ptr);

  virtual void writeStep(std::ostream& flux);

  virtual void checkFreeBoundaries(Array<int>& outOfFreeBounds);

  virtual void stopWalls();
  
  virtual void fillBuffer(ParticleOutput* buffer);
  virtual void fillBuffer(ParticleInSitu* buffer);
  virtual void fillGhostBuffer(ParticleInSitu* buffer);
  virtual void fillBuffer(Array<LegacyParticleIOStruct>& buffer);
  virtual void fillBuffer(Array<LegacyDPDEParticleIOStruct>& buffer);
  virtual void fillBuffer(HerculeParticleIODumpStruct& buffer);
  virtual void fillBuffer(HerculeDPDEParticleIODumpStruct& buffer);

  virtual void fillBuffer(CellOutput* buffer);
  
  virtual void setCallbackQueryFunctions(LBS::LoadBalancer* ptr);

  virtual void printInfoBase(std::ostream& flux);
  virtual void printInfoSpec(std::ostream& flux);
  virtual void print(std::ostream& flux);

  InputOutputManager* m_inputOutputManager; ///< I/O Manager of the domain

  virtual bool checkVerlet();
  virtual void refineBalance();

private:

  void internalReorganization();

  void exchangeParticles();
  void collectParticles();

  void computeForces(Traversal traversal);
  void computeForcesVerlet(Traversal traversal);
  void initForcesSDPD();

  void getProcToComm(Array<uint>& procsToComm, int* importProcs, const uint numImport, int* exportProcs, const uint numExport);
  void myreintegrateParticles();

  static GlobalInfo& globalInfo;  ///< Global attributes about domains gathered in a class

  CommManager* ref;  ///< Communication manager

  Grid<Grid_impl>* grid;  ///< Grid gathering all particles

  keepAtom m_toKeep;
  exchangeGhost m_balanceComm;

};


TMPLG GlobalInfo& TMPL_Domain::globalInfo = Global::domainInfo;

#include "domain/domain.hxx"


#undef TMPLG
#undef TMPL_Domain

#endif // __DOMAIN_HPP_INCLUDED
