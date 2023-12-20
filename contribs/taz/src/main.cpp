/// @file
/// @brief Main function and associated functions
///
/// This file contain the main function, functions used to initialize the parallelization and the function that run the simulation



#include <cassert> // ASSERT AMR

#include <thread>

#include "compile.hpp"

#include "referenceMap.hpp"

#include "eos/allEOS.hpp"

#include "io/input.hpp"
#include "io/outputManager.hpp"

#include "parallel/node.hpp"
#include "parallel/thread/thread.hpp"

#include "particle/molecule/initMolecules.hpp"
#include "particle/types/allTypes.hpp"

#include "potential/allPotentials.hpp"

#include "time/time.hpp"
#include "time/timeIntegration.hpp"


bool updateVerletLists;

void run(Node* node, const std::string& inputfile);


/// @brief The main function
///
/// Check the arguments, initialize the MPI processes, run the calculation and finalize the MPI processes.
/// @param [in] argc Number of arguments.
/// @param [in] argv Arguments.
/// @return EXIT_SUCCESS
int main (int argc, char** argv) {

  // Process command-in-line arguments
  if (argc!=2) {

    std::cerr<< "Usage : mpirun -np <numprocs> " 
	     << "bin/stamp <inputfile>"
	     << std::endl;

    exit(-1);

  }

  // This variable is not used if the linked cell method is used.
  // Otherwise it must be equal to true to build neighbour lists during initialization.
  updateVerletLists=true;

  // Get input file
  std::string inputfile(argv[1]);
  
  // Initialize MPI processes and get the Node
  Node* node = createNode(&argc, &argv);

  ptM = node->getMetrics();

  // Run a simulation
  run(node, inputfile);

  // Finalize MPI processes
  delete node;

  return EXIT_SUCCESS;

}


Input& tmp_load(InputOutputManager* _ioMgr, const std::string& inputfile);
void   tmp_test(Node* node, const Input& input);
void   tmp_threads(Input& input);


class TypeParticle;
class IPotential;
class IEOS;


/// @brief Run the simulation
/// @param [in,out] node The process node
/// @param [in] inputfile Entry file
void run(Node* node, const std::string& inputfile) {

	  Configuration<InputOutputManager> outputConfig;
	  CommManager* commManager=&(node->getManager());
	  InputOutputManager *inputOutput = new InputOutputManager(commManager,outputConfig);

  // Read from file and convert to local units
  Input& input = tmp_load(inputOutput, inputfile);

  // Modification if Verlet Lists are used
  Global::reference.neighbours_method(true, input.rVerlet);
  if(input.neighMethod == Input::NeighMethod::BLOCK_VERLET) Global::reference.setIsBlockVerlet();

  // Test the node repartition
  tmp_test(node, input);

  // Initialization of OpenMP
  tmp_threads(input);

  omp_set_num_threads(Global::maxNumberOfThreads);
  omp_set_schedule(omp_sched_guided,-1); // Tests réalisée durant la thèse Raphaël PRAT

  thread_observer observer; 
  observer.init(); 
  observer.observe(input.bindThreads);
  
  // Set the IO manager
  Configuration<InputOutputManager> outputConfig2(input);
  inputOutput->changeConf(outputConfig2);
  node->setInputOutput(inputOutput);

  // Extract and configure the type of node
  Configuration<Node> nodeConfig(input);
  node->configure(nodeConfig);


  // Extract optional timers chosen
  if(input.isTimers)
  {
    (node->getMetrics())->data.usePotential(  input.isTimersPotential  );
    (node->getMetrics())->data.useNeighbours( input.isTimersNeighbours );
    (node->getMetrics())->data.useGhost(      input.isTimersGhost      );
    (node->getMetrics())->data.useRefine(     input.isTimersRefine     );
    //(node->getMetrics())->data.useRefine(     input.isTimersGeneric    );
  }

  // Extract the types, potentials, equations of state, reaction kinetics and molecules configurations from the input
  Configuration<TypeParticle> typesConfig(input);
  Configuration<IPotential> potentialConfig(input);
  Configuration<IEOS> eosConfig(input);
  Configuration<IKinetics> reactionConfig(input);
  Configuration<MPI__InMol> moleculesConfig(input);

  // If the initialization is not done from a stampv4 input file,
  // configure the reference map
  if(input.m_initType!=INIT_STAMPV4) {
    Global::reference.configure(typesConfig);
    Global::reference.configure(potentialConfig);
  }




  // Extract the domain configuration
  Configuration<DomainInterface> domainConfig(input);
  
  // Configure equation of state and other DPD stuff
  Global::reference.configure(eosConfig);
  Global::reference.configure(reactionConfig);
  Global::reference.configure(input); // Must be done after configuring the potentials

  // Configure domains global data
  Global::domainInfo.configure(domainConfig, node->numberOfNodes(), node->rank());
  // Extract the particles configuration
  Configuration<Particle> particles(input);
  // Create the domains
  node->createDomains(domainConfig, particles, moleculesConfig);

  // Configure time
  Configuration<TimeManager> timeConfig(input);
  node->createTime(timeConfig);

  // Configure the integration scheme
  Configuration<NumericalScheme> schemeConfig(input);
  node->createScheme(schemeConfig);
      
  // Free input data
  delete &input;

  // Main time loop
  node->doComputeWork();  

#if __use_orchestrator
  // Wait for orchestrator thread to complete
  orchestrator_thread.join();
#endif /* __use_orchestrator */

}


/// @brief Read the input file and create the local input data
///
/// @param _ioMgr Pointer to the Input/Output manager
/// @param inputfile Name of the input file
/// @return Pointer to the local input data
Input& tmp_load(InputOutputManager* _ioMgr, const std::string& inputfile) {

  Input* ptrInput = new Input(_ioMgr);

  std::string datastream;
  _ioMgr->load(inputfile, datastream, ptrInput->comment);
  readFromString(datastream, ptrInput);

  return *ptrInput;

}


/// @brief Test the nodes repartition
///
/// @param node The node
/// @param input The input data
void tmp_test(Node* node, const Input& input) {

  if (product(input.decoupage) != node->numberOfNodes()) {

    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'tmp_test(Node*, const Input&)' : Number of procs in each dimension ("
	     << input.decoupage << ") does not match total number of procs ("
	     << node->numberOfNodes() << "). STOP." << std::endl;

    exit(-1);

  }

}


/// @brief TBB setup
///
/// @param [in,out] input The input data
void tmp_threads(Input& input) {

  int nthreads = input.maxNumberOfThreadsPerNode;
  int nsimuthreads = input.maxNbSimulationWorkers;
  int nthreads_default = default_num_threads();

  if (nthreads<1 || nthreads>nthreads_default)
      nthreads = nthreads_default;

  if (nsimuthreads<1 || nsimuthreads>nthreads_default)
      nsimuthreads = nthreads_default;

  input.maxNumberOfThreadsPerNode = nthreads;
  input.maxNbSimulationWorkers = nsimuthreads;

  Global::maxNumberOfThreads = nthreads;

}
