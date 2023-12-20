/// @file
/// @brief Definition of class DomainInterface

#ifndef __DOMAIN_INTERFACE_HPP_INCLUDED
#define __DOMAIN_INTERFACE_HPP_INCLUDED


#include "parallel/balance.hpp"
#include "particle/lattice.hpp"

#include "utils/mat3/mat3.hpp"

#include <yaml-cpp/yaml.h>
#include <yaml_check_particles.h>

class CellOutput;
class HerculeParticleIODumpStruct;
class LegacyParticleIOStruct;
class MPI__Particle;
class ParticleInSitu;
class ParticleOutput;


/// @brief Interface gathering all services called on a domain
class DomainInterface {

public:

  /// @brief Destructor (nothing to do)
  virtual ~DomainInterface() {}

  /// @brief Initialize device if running on GPU or something else
  virtual void initDevice() = 0;

  /// @brief Initialize a chunk of the particles from a configuration
  /// @param [in] particles Configuration of all the particles
  /// @param [in] cBegin First particle to initialize
  /// @param [in] cEnd End of the particle to initialize
  virtual void initParticles(Configuration<MPI__Particle>& particles, uint64_t cBegin, uint64_t cEnd) = 0;

  /// @brief Initialize a chunk of the particles from an array of particles
  /// @param [in] particles Particles
  /// @param [in] initialT Initial temperature
  /// @param [in] numP Number of particles
  /// @param [in] nCells Number of duplicates of the base cell in each dimension, default 1
  /// @param [in] cellSize Size of the base cell, default 0 (if not used)
  /// @param [in] numMol Number of molecules in the base cell, default 0
  virtual void initParticles(MPI__Particle** particles, double initialT, uint64_t numP, vec3<int> nCells=vec3<int>(1), vec3<double> cellSize=vec3<double>(0.), uint64_t numMol=0) = 0;

  /// @brief Update particles positions using to first order approximation of the movement equations
  /// @param [in] time Time step
  virtual void pushPositions1stOrder (double time) = 0;
  /// @brief Update particles positions using to second order approximation of the movement equations
  /// @param [in] time Time step
  virtual void pushPositions2ndOrder (double time) = 0;
  /// @brief Update particles velocities using to first order approximation of the movement equations
  /// @param [in] time Time step
  virtual void pushVelocities1stOrder(double time) = 0;

  /// @brief Update particles velocities by applying Langevin dissipation
  /// @param [in] time Time step
  /// @param [in] gamma Gamma parameter of the Langevin scheme
  /// @param [in] beta Inverse temperature for the Langevin scheme
  virtual void pushDissipationLangevin(double time, double gamma, double beta) = 0;

  /// @brief Make neighbor lists
  virtual void makeNeighborLists() = 0;
  /// @brief Make neighbor lists inside the domain
  virtual void makeNeighborListsInside() = 0;
  /// @brief Make neighbor lists on the domain edges
  virtual void makeNeighborListsOnEdges() = 0;
  /// @brief Clear neighbor lists
  virtual void clearNeighborLists() = 0;

  virtual void aferrgreg(uint64_t& send, uint64_t & recv)=0;

  /// @brief Initialize forces, energies and neighbors
  /// @param [in] doInitEint True if the internal energies need to be initialized
  /// @param [in] initialTint Initial internal temperature
  virtual void initForces(bool doInitEint, double initialTint) = 0;
  
    /// @brief Compute forces with verlet lists
  virtual void doComputeForcesVerlet() = 0;
    /// @brief Compute forces inside the domain
  virtual void doComputeForcesInsideVerlet() = 0;
  /// @brief Compute forces on the domain edges
  virtual void doComputeForcesOnEdgesVerlet() = 0;

  /// @brief Exchange particles between the domains and update cells
  virtual void refineCells() = 0;

  /// @brief Exchange particles between the domains and update cells
  virtual void updateCells() = 0;
  /// @brief Exchange particles between the domains and update cells
  virtual void updateCellsInside() = 0;
  /// @brief Clear the ghost cells
  virtual void clearGhost() = 0;
  /// @brief Update ghost cells
  virtual void updateGhost() = 0;
  /// @brief Collect ghost after a communication
  virtual void collectGhost() = 0;

  /// @brief Get number of real particles on the domain
  /// @return Number of particles
  virtual uint64_t getNumberOfParticles()      = 0;
  /// @brief Get number of real cells on the domain
  /// @return Number of particles
  virtual uint64_t getNumberOfRealCells()      = 0;

  /// @brief Get the total energy (kinetic plus potential) of the particles on the domain
  /// @return Total energy
  virtual double getTotalEnergy()     = 0;
  /// @brief Get the kinetic energy of the particles on the domain
  /// @return Kinetic energy
  virtual double getKineticEnergy()   = 0;
  /// @brief Get the potential energy of the particles on the domain
  /// @return Potential energy
  virtual double getPotentialEnergy() = 0;
  /// @brief Get the internal energy of the particles on the domain
  /// @return Internal energy
  virtual double getInternalEnergy() = 0;
  /// @brief Get the chemical energy of the particles on the domain
  /// @return chemical energy
  virtual double getChemicalEnergy() = 0;

  /// @brief Get the total momentum of the particles on the domain
  /// @return Total momentum
  virtual vec3<double> getTotalMomentum() = 0;
  /// @brief Get the total mass of the particles on the domain
  /// @return Total mass
  virtual double getTotalMass() = 0;
  /// @brief Get the total shifted kinetic energy on the domain
  /// @return Shifted kinetic energy
  virtual double getShiftedKineticEnergy(const vec3<double>& vshift) = 0;

  /// @brief Get the total pressure on the domain
  /// @return Pressure tensor
  virtual mat3<double> getPressure(const vec3<double>& vshift) = 0;

  /// @brief Compute workload on each particle
  /// @return Total workload on the domain
  virtual double workload() = 0;
  /// @brief Balance workload on domains
  /// @param [in] ptr Load balancer (pointer)
  /// @return State of the load balancer at the end
  virtual LBS::State balance(LBS::LoadBalancer* ptr) = 0;

  /// @brief Write domain data at current step in specified flux
  /// @param [in,out] flux Print flux
  virtual void writeStep(std::ostream& flux) = 0;

  /// @brief Check if there is escaped particles at free boundaries
  /// @param [out] outOfFreeBounds Number of escaped particles for each free boundary
  virtual void checkFreeBoundaries(Array<int>& outOfFreeBounds) = 0;

  /// @brief Stop the wall particles
  virtual void stopWalls() = 0;

  /// @brief Fill a buffer with all the particles of the node, case of a ParticleOutput buffer
  /// @param [out] buffer Buffer to fill
  virtual void fillBuffer(ParticleOutput* buffer) = 0;
  /// @brief Fill a buffer with all the particles of the node, case of a ParticleInSitu buffer
  /// @param [out] buffer Buffer to fill
  virtual void fillBuffer(ParticleInSitu* buffer) = 0;
  /// @brief Fill a buffer with the ghost information of the node, case of a ParticleInSitu buffer
  /// @param [out] buffer Buffer to fill
  virtual void fillGhostBuffer(ParticleInSitu* buffer) = 0;
  /// @brief Fill a buffer with all the particles of the node, case of a legacy Stamp dump buffer
  /// @param [out] buffer Buffer to fill
  virtual void fillBuffer(Array<LegacyParticleIOStruct>& buffer) = 0;
  /// @brief Fill a buffer with all the particles of the node, case of a legacy Stamp dump buffer for DPDE
  /// @param [out] buffer Buffer to fill
  virtual void fillBuffer(Array<LegacyDPDEParticleIOStruct>& buffer) = 0;
  /// @brief Fill a buffer with all the particles of the node, case of a Hercule dump buffer
  /// @param [out] buffer Buffer to fill
  virtual void fillBuffer(HerculeParticleIODumpStruct& buffer) = 0;
  /// @brief Fill a buffer with all the particles of the node, case of a Hercule dump buffer for DPDE
  /// @param [out] buffer Buffer to fill
  virtual void fillBuffer(HerculeDPDEParticleIODumpStruct& buffer) = 0;

  /// @brief Fill a buffer with all the cells of the node
  /// @param [out] buffer Buffer to fill
  virtual void fillBuffer(CellOutput* buffer) = 0;

  /// @brief Set Zoltan load balancing parameters
  /// @param [in,out] ptr Load balancer (pointer)
  virtual void setCallbackQueryFunctions(LBS::LoadBalancer* ptr) = 0;

  /// @brief Print domains local info of the domain in specified flux
  /// @param [in,out] flux Print flux
  virtual void printInfoBase(std::ostream& flux) = 0;
  /// @brief Print grid info in specified flux ?
  /// @param [in,out] flux Print flux
  virtual void printInfoSpec(std::ostream& flux) = 0;
  /// @brief Debug print for all cells in specified flux (not used)
  /// @param [in,out] flux Print flux
  virtual void print        (std::ostream& flux) = 0;

  /// @brief Check if the verlet lists have to be updated
  virtual bool checkVerlet() =0;
  virtual void refineBalance() =0;

  virtual void ctest(std::set<ParticleReferenceValue> &reference_values_set, double &ae, double &re) = 0;

};


class Node;


DomainInterface* buildDomain(Node* node, uint index, Configuration<DomainInterface>& config, uint64_t numberOfParticles);


DomainInterface* buildDomainSOTL(Node* node, uint index, Configuration<DomainInterface>& config, uint numberOfParticles);

#endif // __DOMAIN_INTERFACE_HPP_INCLUDED
