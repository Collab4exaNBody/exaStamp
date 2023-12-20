/// @file 
/// @brief Interface for Node object, and one description of it

#ifndef __NODE_HPP_INCLUDED
#define __NODE_HPP_INCLUDED


#include "parallel/balance.hpp"
#include "parallel/mympi.hpp"
#include "parallel/metrics.hpp"

#include "io/input.hpp"

#if __use_orchestrator
#include "orchestrator/shared.hpp"
#endif /* __use_orchestrator */


class CellOutput;
class DomainInterface;
class MPI__InMol;
class MPI__Particle;
class Node;
class NumericalScheme;
class ParticleOutput;
class TimeManager;

/// @brief see MPI particle
typedef MPI__Particle Particle;


/// @brief Temporary structure gathering all data to configure a node
///
/// This structure is a base to choose the domain internal structure.
template <> struct Configuration<Node> {

  /// @brief Enumeration of the types of node
  enum Type {
	  FLAT, ///< Flat node (no accelerator)
	  SOTL ///< SOTL node (with accelerators)
  };

  // Later : add other types (HYBRID, ...), maybe subtypes, number of 
  // threads, gpu, mic ... 
  // This struct is a base to choose the domain internal structure !

  /// @brief The default constructor
  ///
  ///
  Configuration() : type(FLAT) {}
  /// @brief Destructor (nothing to do)
  ~Configuration() {}

  Configuration(const Input& input);

  Type type; ///< Node type
  int maxNumberOfThreads; ///< Max number of threads

  LBS::Method balancingMethod; ///< Load balancing method
  bool logWorkload; ///< Indicates if workload is printed in a special log file

  bool expandFreeLimits; ///< Allows to expand free boundary limit
  double expandFactor; ///< Expand factor for domain in dump in case of expansion

};


/// @brief Store the parallelization setup of a node
struct NodeSetup {

  /// @brief Default constructor
  ///
  ///
  NodeSetup() : numberOfThreads(1), numberOfGPUs(0), numberOfMics(0), useSOTL(false), expandFreeLimits(false) {}
  /// @brief Destructor (nothing to do)
  ~NodeSetup() {}

  int numberOfThreads; ///< Number of threads
  int numberOfGPUs; ///< Number of GPU accelerators
  int numberOfMics; ///< Number of Mics accelerators

  bool useSOTL; ///< True if SOTL library is used
  bool logWorkload; ///< Indicates if workload is printed in a special log file
  bool expandFreeLimits; ///< Allows to expand free boundary limits
  double expandFactor; ///< Expand factor for domain in dump in case of expansion

};


/// @brief The node class.
///
/// The node is the top structure of the code. It contains the integration scheme, a list of domains and a communication manager.
class Node {

public:

  Node(MPI_Comm comm);

  virtual ~Node();

  int rank();
  int numberOfNodes();

  bool hasAccelerator();
  bool useSOTL();
  
  /// @brief Return number of domains on this node
  virtual uint getNumberOfDomains() = 0;
  /// @brief Return pointer to domains' interfaces
  virtual DomainInterface** getDomains() = 0;

  CommManager& getManager();

  // Not implemented
  //void load(const std::string& inputfile, std::string& datastream, const char comment);
  void configure(Configuration<Node>& configuration);
  void createScheme(Configuration<NumericalScheme>& configuration);
  void createTime(Configuration<TimeManager>& configuration);
  void createInputOutput(Configuration<InputOutputManager>& configuration);
  void setInputOutput(InputOutputManager* _ioMgr);

  /// @brief compare computed values with reference values stored in "check_values" file
  virtual void verifySimulationValues(uint64_t timeStep) =0;

  /// @brief Create domains from domains configuration and particles configuration
  /// @param [in] configuration Configuration of domains
  /// @param [in] particles Configuration of particles (if no molecules)
  /// @param [in] molecules Configuration of the molecules (if there is some)
  virtual void createDomains(Configuration<DomainInterface>& configuration, Configuration<Particle>& particles, Configuration<MPI__InMol>& molecules) = 0;
  /// @brief Do all the work
  virtual void doComputeWork() = 0;

  virtual void aferrgreg(uint64_t& send, uint64_t & recv)=0;

  /// @brief Reduce energies on all nodes
  /// @return Array with kinetic energy, potential energy, internal energy and shifted kinetic energy
  virtual Array<double> reduceEnergy() = 0;

  /// @brief Reduce number of particles on all nodes
  /// @return Number of particles
  virtual uint64_t reduceNumberOfParticles() = 0;

  /// @brief Reduce pressure on all nodes
  /// @return Pressure
  virtual double reducePressure() = 0;

  /// @brief Calculate imbalance rate
  /// @param [in] write Print imbalance rate if true
  /// @param [in] step Step where imbalance rate is printed
  /// @return Imbalance rate
  virtual double imbalanceRate(bool write=false, int step=0) = 0;

  ///@brief Accessor to the input/ouput manager
  InputOutputManager*   inputOutput(){return m_inputOutput;}

  ///@brief Accessor to the time manager
  TimeManager* timeManager() { return time; }

  ///@brief Accessor to the metrics
  Metrics* getMetrics() { return &metrics; }
protected:
  void print(std::ostream& flux);

  void writeDump(std::ostream& flux, int step, double dt, double physTime,const std::string& _rep, DumpType type);
  void writeEnergies(std::ostream& flux, int step, double dt, LBS::State state);
  void writeParticles(int step, double cpuTime, double physTime, const std::string& _rep,OutputType type);
  void writeCells(int step, double cpuTime, double physTime, const std::string& _rep,OutputType type);
  void dumpParticles(int step, double cpuTime, double physTime, const std::string& _rep,DumpType type);
  void writeEndLine(std::ostream& flux=std::cout, int n=1);

  /// @brief Gather particles on all nodes
  /// @param [in] step Step where the gathering is done
  /// @param [in] doTypes Indicates if the types are written
  /// @param [in] doVelocities Indicates if the velocities are written
  /// @param [in] doEints Indicates if the internal energies are written
  /// @param [in] doProg Indicate if the progress variables are written
  virtual ParticleOutput* gatherParticles(int step, bool doTypes, bool doVelocities, bool doEints, bool doProg) = 0;

  /// @brief Gather cells on all nodes
  /// @param [in] step Step where the gathering is done
  /// @return A CellOutput containing all particles
  virtual CellOutput* gatherCells(int step) = 0;
  /// @brief Create data for legacy dump writing (default version)
  /// @param [out] header Created header
  /// @param [out] particlesArray Created particles array
  /// @param [in] step Step where the dump writing is done
  /// @param [in] cpuTime CPU time
  /// @param [in] physTime Physical time
  virtual void setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,Array<LegacyParticleIOStruct>& particlesArray, int step, double cpuTime, double physTime) = 0;
  /// @brief Create data for legacy dump writing (DPDE version)
  /// @param [out] header Created header
  /// @param [out] particlesArray Created particles array
  /// @param [in] step Step where the dump writing is done
  /// @param [in] cpuTime CPU time
  /// @param [in] physTime Physical time
  virtual void setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,Array<LegacyDPDEParticleIOStruct>& particlesArray, int step, double cpuTime, double physTime) = 0;
  /// @brief Create data for legacy dump writing (Hercule version)
  /// @param [out] header Created header
  /// @param [out] particle Created particles array
  /// @param [in] step Step where the dump writing is done
  /// @param [in] cpuTime CPU time
  /// @param [in] physTime Physical time 
 virtual void setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,HerculeParticleIODumpStruct& particle, int step, double cpuTime, double physTime) = 0;
  /// @brief Create data for legacy dump writing (Hercule version for DPDE)
  /// @param [out] header Created header
  /// @param [out] particle Created particles array
  /// @param [in] step Step where the dump writing is done
  /// @param [in] cpuTime CPU time
  /// @param [in] physTime Physical time
  virtual void setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,HerculeDPDEParticleIODumpStruct& particle, int step, double cpuTime, double physTime) = 0;


  /// @brief Print domains local info (like symmetrization)
  /// @param [in,out] flux Print flux
  virtual void printLocal(std::ostream& flux) = 0;

  NodeSetup setup; ///< Parallelization setup
  Metrics metrics; ///< Time and memory usage measure tool

  CommManager commManager; ///< Communication manager
  bool isMaster; ///< True if master node

  TimeManager*     time; ///< Time manager
  NumericalScheme* scheme; ///< Integration scheme
  InputOutputManager*   m_inputOutput; ///< Input/output manager

  LBS::LoadBalancer loadBalancer; ///< Load balancer configuration
  Array<int> m_outOfFreeBounds; ///< Store the number of particles out of each free boundary

};


/// @brief The class for a node with only one domain
class NodeSingleDomain : public Node {

public:

  /// @brief Constructor from a MPI communicator
  /// @param [in,out] comm MPI communicator, default=MPI_COMM_WORLD
  NodeSingleDomain(MPI_Comm comm=MPI_COMM_WORLD) 
  : Node(comm), domain(nullptr) {}

  virtual ~NodeSingleDomain();

  /// @brief Return number of domains on this node
  virtual uint getNumberOfDomains() { return 1; };
  /// @brief Return pointer to domains' interfaces
  virtual DomainInterface** getDomains() { return &domain; }

  virtual void createDomains(Configuration<DomainInterface>& configuration, Configuration<Particle>& particles, Configuration<MPI__InMol>& molecules);
  virtual void doComputeWork();

  virtual Array<double> reduceEnergy();
  virtual void aferrgreg(uint64_t& send, uint64_t & recv);
  virtual uint64_t reduceNumberOfParticles();
  virtual uint64_t reduceNumberOfCells();
  virtual double reducePressure();
  virtual double imbalanceRate(bool write=false, int step=0);

  /// @brief compare computed values with reference values stored in "check_values" file
  void verifySimulationValues(uint64_t timeStep) override final;
  

protected:

  void updateCells(bool& upToDate);

  void oneStep();
  void balance(LBS::State& state, bool& cellsUpToDate);
  void writeIO(LBS::State& state, bool& cellsUpToDate, double& tmpTime);
  void endIteration(LBS::State& state, bool& cellsUpToDate, double& tmpTime);
  void goAheadChecking(bool& cellsUpToDate, double& tmpTime);
  void first_balance(LBS::State& state, bool& cellsUpToDate);
  void wallChecking();
  void verletChecking(bool& collect);
  void collectTimer();

  void verifySimulationYAML(uint64_t timeStep);

  virtual void printLocal(std::ostream& flux);

  virtual ParticleOutput* gatherParticles(int step, bool doTypes, bool doVelocities, bool doEints, bool doProg);
  virtual void setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,Array<LegacyParticleIOStruct>& particlesArray, int step, double cpuTime, double physTime);
  virtual void setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,Array<LegacyDPDEParticleIOStruct>& particlesArray, int step, double cpuTime, double physTime);
  virtual void setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,HerculeParticleIODumpStruct& particles, int step, double cpuTime, double physTime);
  virtual void setLegacyDumpData(Array<LegacyHeaderIOStruct>& header,HerculeDPDEParticleIODumpStruct& particles, int step, double cpuTime, double physTime);

  virtual CellOutput* gatherCells(int step);
  
  DomainInterface* domain; ///< Interface to the domain

};


/// @brief The class for a node with several domains
///
/// Not implemented yet
class NodeMultiDomain : public Node {

};


// Initialize MPI process and get the nodes
Node* createNode(int* argc_, char*** argv_);

#endif // __NODE_HPP_INCLUDED
