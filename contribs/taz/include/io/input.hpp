/// @file 
/// @brief Definition of the input class

#ifndef __INPUT_HPP_INCLUDED
#define __INPUT_HPP_INCLUDED


#include <iostream>
#include <iomanip>
#include <string>


#include "io/io.hpp"
#include "io/StampV3LegacyIOStructures.hpp"

#include "utils/array/array.hpp"
#include "utils/vec3/vec3.hpp"

#include "omp.h"

class InputOutputManager;

/// @brief Store the whole input data
class Input {

public:

  const char comment = '#'; ///< Comment delimiter in the input file
  const char assign  = '='; ///< Assignment symbol in the input file

  const char string_delimiter = '\''; ///< String delimiter (not used)
  
  const char array_delimiter_left = '['; ///< Array left delimiter in the input file
  const char array_delimiter_right = ']'; ///< Array right delimiter in the input file
  const char array_separator = ','; ///< Array values separator in the input file
  
  const std::string env_flag     = "XSP_ENV"; ///< Flag indicating that a value must be read from an environment variable
  
  const std::string env_recover  = "XSP_RECOVERING"; ///< Name of the environment variable for the recovering
  const std::string env_symm      = "XSP_SYMMETRIZE"; ///< Name of the environment variable for the symmetrization
  const std::string env_sym      = "XSP_SYMETRIZE"; ///< Name of the environment variable for the symetrization (symmetrization with french spelling mistake)
  const std::string env_nthreads = "XSP_MAX_NUM_THREADS"; ///< Name of the environment variable for the maximum number of threads per node

  const std::string env_np_x     = "XSP_NPX"; ///< Name of the environment variable for the number of MPI processes in the x dimension
  const std::string env_np_y     = "XSP_NPY"; ///< Name of the environment variable for the number of MPI processes in the y dimension
  const std::string env_np_z     = "XSP_NPZ"; ///< Name of the environment variable for the number of MPI processes in the z dimension

  const std::string env_dcmp     = "XSP_DECOMPOSITION"; ///< Name of the environment variable for the type of decomposition
  const std::string env_lb_mthd  = "XSP_LB_METHOD"; ///< Name of the environment variable for the load balancing method

  const std::string env_log_rate  = "XSP_LOG_RATE"; ///< Name of the environment variable for the log writing rate
  const std::string env_lb_rate   = "XSP_LB_RATE"; ///< Name of the environment variable for the load balancing rate
  const std::string env_out_rate  = "XSP_OUT_RATE"; ///< Name of the environment variable for the output writing rate
  const std::string env_out_type  = "XSP_OUT_TYPE"; ///< Name of the environment variable for the output writing type
  const std::string env_out_dir  = "XSP_OUT_DIR"; ///< Name of the environment variable for the output writing directory
  const std::string env_dump_rate = "XSP_DUMP_RATE"; ///< Name of the environment variable for the dump writing rate
  const std::string env_dump_type = "XSP_DUMP_TYPE"; ///< Name of the environment variable for the dump writing type
  const std::string env_dump_dir = "XSP_DUMP_DIR"; ///< Name of the environment variable for the dump writing directory
  const std::string env_initialization  = "XSP_INITIALIZATION"; ///< Name of the environment variable for the initialization type
  const std::string env_init_step  = "XSP_INIT_STEP"; ///< Name of the environment variable for the initial step
  const std::string env_particles_file  = "XSP_PARTICLES_FILE"; ///< Name of the environment variable for the particles initialization file
  const std::string env_n_steps  = "XSP_N_STEPS"; ///< Name of the environment variable for the number of steps

  const std::string env_nc_x     = "XSP_NCX"; ///< Name of the environment variable for the number of lattice cells in the x dimension
  const std::string env_nc_y     = "XSP_NCY"; ///< Name of the environment variable for the number of lattice cells in the y dimension
  const std::string env_nc_z     = "XSP_NCZ"; ///< Name of the environment variable for the number of lattice cells in the z dimension

  /// @brief Enumeration of the types of node
  enum NodeType {
  	FLAT, ///< Default
  	SOTL ///< With GPU
  };
  /// @brief Enumeration of the possible types of simulation
  enum SimulationMode {
  	SINGLE_SPEC_ATOM, ///< Only one class of particles, atoms
  	SINGLE_SPEC_MESO, ///< Only one class of particles, mesoparticles
  	SINGLE_SPEC_SMOOTH, ///< Only one class of particles, smooth particles
  	SINGLE_SPEC_STIFF_MOLEC, ///< Only one class of particles, stiff molecule
  	ATOM_IN_MOL ///< Grid + molecule with only atoms as point particles
  };

  /// @brief Enumeration of the types of decomposition
  enum DecompositionType {
  	RECTILINEAR, ///< Rectilinear decomposition
  	ANY ///< Nonspecific decomposition
  };
  /// @brief Enumeration of the boundary conditions
  enum BoundaryCondition {
  	FREE, ///< Free (no position constraint for the particles)
  	PERIODIC, ///< Periodic (min bound = max bound)
  	WALL , ///< Wall force at the boundary to keep the particles inside
  	FREE_WALL, ///< Free at min bound and wall at max bound
  	WALL_FREE ///< Wall at min bound and free at max bound
  };
  /// @brief Shortcut for boundary conditions in the three dimension
  typedef vec3<BoundaryCondition> BoundaryConditions;

  /// @brief Enumeration of the integration schemes
  enum IntegrationScheme {
  	NEWTON_VERLET_LEAPFROG, ///< Leapfrog integration scheme
  	NEWTON_VERLET_VELOCITY, ///< Verlet Velocity integration scheme
  	LANGEVIN,  ///< Langevin Splitting integration scheme
	DPD, ///< Splitting scheme for DPD
	DPDE_SER, ///< SER scheme for DPDE
	SDPD_VV_SER, ///< SER scheme for SDPD
	SDPD_LANGEVIN_SER ///< SER+Langevin scheme for SDPD
  };
  /// @brief Enumeration of the load balancing methods
  enum BalancingMethod {
  	NONE, ///< No load balancing
  	BLOCK, ///< Block load balancing (designed for testing, cannot handle the double weights of ExaStamp)
  	RANDOM, ///< Random load balancing (for real; do not use)
  	COORD_BSC, ///< Recursive coordinate bisection (cut system in half according to the weights until the number of process is reached)
  	INERT_BSC, ///< Recursive inertial bisection (same as Recursive Coordinate Bisection but without axis constraint)
  	SPC_FILL_CRV, ///< Hilbert space filling curve
  	PHG, ///< Parallel hypergraph partitioning (balancing method that consider edges)
	METIS,
	SCOTCH
  };
  /// @brief Enumeration of the lattice types
  enum LatticeType {
  	FCC, ///< Face-centered cubic lattice
  	BCC, ///< Body centered cubic lattice
  	SC, ///< Simple cubic lattice
  	DIAM100 ///<
  };
  /// @brief Ennumeration of the force fields
  enum ForceField {
  	UFF, ///< General force field with parameters for the full periodic table, open source, no Coulomb interaction
  	COMPASS, ///< Condensed phase optimized force field, not fully available, simple charges
  	AMBER ///< Force field widely used for proteins and DNA, open source
  };
  /// @brief Enumeration of the methods to get the partial charges of atoms in molecules
  enum ChargeMethod {
  	NO_CHARGES, ///< All partial charges to zero
  	FROM_COMPASS, ///< Sum of bonds contributions
  	OPENBABEL, ///< OpenBabel default charges
  	NOT_SPECIFIED ///< Not explicitly defined
  };

  /// @brief Enumeration of the synchronous policy
  enum SynchronousPolicy {
    SYNCHRONOUS_WITHOUT_ARENAS,  ///< synchronous (same size for the arenas but simulation blocked during analytics)
    SYNCHRONOUS_WITH_ARENAS, ///< synchronous (different sizes for the arenas and simulation blocked during analytics)
    ASYNCHRONOUS_WITHOUT_ARENAS, ///< asynchronous (same size for the arenas but simulation unblocked during analytics)
    ASYNCHRONOUS_WITH_ARENAS ///< asynchronous (different sizes for the arenas and simulation unblocked during analytics)
  };

  /// @brief Enumeration of methods of constructing neighbour lists
  enum NeighMethod {
    VERLET_LIST,   ///< verlet list
    BLOCK_VERLET
  };

  Input(InputOutputManager* );

  /// @brief Destructor
  ///
	///
  ~Input() {
    if (legacyHeader!=nullptr) delete legacyHeader;
  }

  void convertToType(const std::string& field, const std::string& value);

  void convertToStampUnits();

  InputOutputManager* inputOutputManager; ///< Input/Output manager

  int version; ///< Version of the code to use

  NodeType nodeType; ///< Type of the nodes

  SimulationMode mode; ///< Type of simulation
  IntegrationScheme scheme; ///< Integration scheme
  double frictionLangevin; ///< Friction parameter for a Langevin scheme
  double temperatureThermostat; ///< Target temperature for thermostats
  int seed; ///< Random number generator seed
  
  bool recovering; ///< Indicates if there will be recovering in the force computation
  bool symmetrization; ///< Indicates if the force computation will be symmetrized
  int maxNumberOfThreadsPerNode; ///< Maximum number of threads per node
  bool bindThreads; ///< Indicates if threads must be bound

  // System extension data
  DecompositionType decomposition; ///< Type of decomposition
  vec3<int> decoupage; ///< Distribution of the MPI processes in the three dimensions
  BoundaryConditions boundaryConditions; ///< Boundary conditions on the system
  vec3<double> origin; ///< Origin of the system
  vec3<double> extension; ///< Extension of the system
  bool expandFreeLimits; ///< Indicates to check free boundaries and stop the simulation if a particle is about to cross one; not mandatory but if you don't some particles might mysteriously disappear
  double expandFactor; ///< Expansion factor to apply to the extension in dump in case where expandFreeLimits is activated

  // Time data
  int numberOfSteps; ///< Number of steps in the simulation
  double delta; ///< Size of a step
  double stopWallsTime; ///< Time when walls should be stopped
  int logRate; ///< Log writing rate (write log each logRate steps)
  int outputRate; ///< Output writing rate (write output each outputRate steps)
  OutputType outputType; ///< Type of output file
  bool outputCells; ///< Indicate if the cell output must be printed
  int dumpRate; ///< Dump writing rate (write dump each dumpRate steps)
  DumpType dumpType; ///< Type of dump file
  std::string dir; ///< Work directory
  std::string dumpDir; ///< Dump file directory
  std::string outputDir; ///< Output file directory
  SynchronousPolicy syncPolicy; ///< Indicates if outputs should be performed synchronously or asynchronously
  Array<std::string> analyticsGraph; ///< List of the analytics asked by the user
  bool isTimers; ///< Indicates if user wants timers more detailed
  bool isTimersPotential; ///< Indicates if user wants timers more detailed
  bool isTimersNeighbours; ///< Indicates if user wants timers more detailed
  bool isTimersGhost; ///< Indicates if user wants timers more detailed
  bool isTimersRefine; ///< Indicates if user wants timers more detailed
  bool isTimersGeneric; ///< Indicates if user wants timers more detailed
  bool doTypes; ///< Indicates if user wants to output the types
  bool doVelocities; ///< Indicates if user wants to output the velocities
  bool doEints; ///< Indicates if user wants to output the internal energy
  int maxNbAnalyticsWorkers; ///< Maximum number of workers for the analytics arena
  int maxNbSimulationWorkers; ///< Maximum number of workers for the simulation arena

  //Verlet list method
  double rVerlet; ///< raduis of verlet list
  NeighMethod neighMethod; ///< Methods of constructing neighbour lists

  // Load balancing data
  BalancingMethod balancingMethod; ///< Load balancing method
  int balanceRate; ///< Load balancing rate (balance workload each balanceRate steps)
  bool logWorkload;  ///< Indicates if workload is printed in a special log file

  // Atom types data
  Array<std::string> atomNames; ///< Names of the atoms used in the simulation
  Array<int> atomAtomicNumbers; ///< Atomic numbers of the atoms used in the simulation
  Array<double> atomMasses; ///< Masses of the atoms used in the simulation
  Array<double> atomCharges; ///< Charges of the atoms used in the simulation

  // Ideal gas interactions data
  Array<std::string> IGtypeA; ///< First elements of the couples of atoms types that interact as ideal gases
  Array<std::string> IGtypeB; ///< Second elements of the couples of atoms types that interact as ideal gases

  // Lennard-Jones interactions data
  Array<std::string> LJtypeA; ///< First elements of the couples of atoms types that interact by a Lennard-Jones potential
  Array<std::string> LJtypeB; ///< Second elements of the couples of atoms types that interact by a Lennard-Jones potential
  Array<double> LJrcut; ///< Cutoff radius for each Lennard-Jones interaction
  Array<double> LJepsilon; ///< Depth of the potential well radius for each Lennard-Jones interaction
  Array<double> LJsigma; ///< Finite distance at which the inter-particle potential is zero radius for each Lennard-Jones interaction

  // Exponential-6 interactions data
  Array<std::string> Exp6typeA; ///< First elements of the couples of atoms types that interact by an Exp. 6 potential
  Array<std::string> Exp6typeB; ///< Second elements of the couples of atoms types that interact by an Exp. 6 potential
  Array<double> Exp6rcut; ///< Cutoff radius for each Exp. 6 interaction
  Array<double> Exp6foo; ///< Some unidentified and unused property for each Exp. 6 interaction

  // Sutton-Chen interactions data
  Array<std::string> SCtypeA; ///< First elements of the couples of atoms types that interact by a Sutton-Chen potential
  Array<std::string> SCtypeB; ///< Second elements of the couples of atoms types that interact by a Sutton-Chen potential
  Array<double> SCrcut; ///< Cutoff radius for each Sutton-Chen interaction
  Array<double> SCc; ///< Shape parameter c for each Sutton-Chen interaction
  Array<double> SCepsilon; ///< Energy scale for each Sutton-Chen interaction
  Array<double> SCa0; ///< Length scale for each Sutton-Chen interaction
  Array<double> SCn; ///< Shape parameter n for each Sutton-Chen interaction
  Array<double> SCm; ///< Shape parameter m for each Sutton-Chen interaction

  // EAM VNIITF interactions data
  Array<std::string> EamVNIITFtypeA; ///< First elements of the couples of atoms types that interact by an EAM VNIITF potential
  Array<std::string> EamVNIITFtypeB; ///< Second elements of the couples of atoms types that interact by an EAM VNIITF potential
  Array<double> EamVNIITFrcut; ///< Cutoff radius for each EAM VNIITF interaction
  Array<double> EamVNIITFrmax; ///< Maximal distance for density contribution for each EAM VNIITF interaction
  Array<double> EamVNIITFrmin; ///< Minimal distance for density contribution for each EAM VNIITF interaction
  Array<double> EamVNIITFrt0; ///< Characteristic radius for each EAM VNIITF interaction
  Array<double> EamVNIITFEcoh; ///< Cohesive energy for each EAM VNIITF interaction
  Array<double> EamVNIITFE0; ///< Energy parameter for each EAM VNIITF interaction
  Array<double> EamVNIITFbeta; ///< Attenuation rate for each EAM VNIITF interaction
  Array<double> EamVNIITFA; ///< Empirical parameter A for each EAM VNIITF interaction
  Array<double> EamVNIITFZ; ///< Number of neighbors in the reference structure for each EAM VNIITF interaction
  Array<double> EamVNIITFn; ///< Empirical parameter n for each EAM VNIITF interaction
  Array<double> EamVNIITFalpha; ///< Parameter alpha for each EAM VNIITF interaction
  Array<double> EamVNIITFD; ///< Parameter D for each EAM VNIITF interaction
  Array<double> EamVNIITFeta; ///< Parameter eta for each EAM VNIITF interaction
  Array<double> EamVNIITFmu; ///< Parameter mu for each EAM VNIITF interaction


 // MEAM
  Array<std::string> MeamTypeA; ///< First elements of the couples of atoms types that interact by an MEAM potential
  Array<std::string> MeamTypeB; ///< Second elements of the couples of atoms types that interact by an MEAM potential

  Array<double> MeamRcut; ///< Cutoff radius for each MEAM interaction
  Array<double> MeamRmax; ///< Maximal distance for density contribution for each MEAM interaction
  Array<double> MeamRmin; ///< Minimal distance for density contribution for each MEAM interaction
  Array<double> MeamEcoh; ///< Cohesive energy for each MEAM interaction
  Array<double> MeamE0; ///< energy for each MEAM interaction
  Array<double> MeamA; ///< Empirical parameter A for each MEAM interaction
  Array<double> MeamR0; ///< Parameter r0 for each MEAM interaction
  Array<double> MeamAlpha; ///< Parameter alpha for each MEAM interaction
  Array<double> MeamDelta; ///< Parameter delta for each MEAM interaction
  Array<double> MeamBeta0; ///< Attenuation rate for each MEAM interaction
  Array<double> MeamBeta1; ///< Attenuation rate for each MEAM interaction
  Array<double> MeamBeta2; ///< Attenuation rate for each MEAM interaction
  Array<double> MeamBeta3; ///< Attenuation rate for each MEAM interaction
  Array<double> MeamT0; ///< Coefficient t0 for each MEAM interaction
  Array<double> MeamT1; ///< Coefficient t1 for each MEAM interaction
  Array<double> MeamT2; ///< Coefficient t2 for each MEAM interaction
  Array<double> MeamT3; ///< Coefficient t3 for each MEAM interaction
  Array<double> MeamS0; ///< Coefficient s0 for each MEAM interaction
  Array<double> MeamS1; ///< Coefficient s1 for each MEAM interaction
  Array<double> MeamS2; ///< Coefficient s2 for each MEAM interaction
  Array<double> MeamS3; ///< Coefficient s3 for each MEAM interaction
  Array<double> MeamCmin; ///< Boundary inf used by screening function for each MEAM interaction
  Array<double> MeamCmax; ///< Boundary sup used by screening function for each MEAM interaction
  Array<double> MeamZ; ///< Number of neighbors in the reference structure for each MEAM interaction
  Array<double> MeamRc; ///< distance f(rc)=0
  Array<double> MeamRp; ///< interval [rc-rp,rc] where we apply a polynome such that P(rc)=0

  // Gaussian interactions data
  Array<std::string> GaussTypeA; ///< First elements of the couples of atoms types that interact by a Gaussian potential
  Array<std::string> GaussTypeB; ///< Second elements of the couples of atoms types that interact by a Gaussian potential
  Array<double> GaussRcut; ///< Cutoff radius for each Gaussian interaction
  Array<double> GaussEpsilon; ///< Force amplitude for each Gaussian interaction
  Array<double> GaussRatt; ///< Attractive radius for each Gaussian interaction
  Array<double> GaussRrep; ///< Repulsive radius for each Gaussian interaction
  Array<double> GaussRatio; ///< Attractive/Repulsive ratio for each Gaussian interaction
  
  // Mesoparticles
  Array<std::string> mesoNames; ///< Names of the mesoparticles used in the simulation
  Array<double> mesoMasses; ///< Masses of the mesoparticles used in the simulation
  Array<double> mesoUnitMasses; ///< Unitary masses (=mass/size) of the mesoparticles used in the simulation
  Array<double> mesoGammaPara; ///< Parallel friction coefficient of the mesoparticles used in the simulation
  Array<double> mesoGammaOrtho; ///< Orthogonal parallel coefficient of the mesoparticles used in the simulation

  // Smmoth particles
  Array<std::string> smoothNames; ///< Names of the smooth particles used in the simulation
  Array<double> smoothMasses; ///< Masses of the smooth particles used in the simulation
  Array<double> smoothUnitMasses; ///< Unitary masses (=mass/size) of the smooth particles used in the simulation
  Array<double> smoothBulkViscosity; ///< Bulk viscosities of the smooth particles used in the simulation
  Array<double> smoothShearViscosity; ///< Shear viscosities of the smooth particles used in the simulation
  Array<double> smoothSmoothingLength; ///< Smoothing lengths of the smooth particles used in the simulation
  Array<std::string> smoothKernel; ///< Smoothing kernels of the smooth particles used in the simulation

  // Smmoth wall particles
  Array<std::string> wallNames; ///< Names of the wall particles used in the simulation
  Array<double> wallMasses; ///< Masses of the wall particles used in the simulation
  Array<double> wallUnitMasses; ///< Unitary masses (=mass/size) of the wall particles used in the simulation
  Array<double> wallSmoothingLength; ///< Walling lengths of the wall particles used in the simulation
  Array<std::string> wallKernel; ///< Smoothing kernels of the wall particles used in the simulation
  Array<double> wallVelocityX; ///< Wall velocity in the x-direction
  Array<double> wallVelocityY; ///< Wall velocity in the y-direction
  Array<double> wallVelocityZ; ///< Wall velocity in the z-direction
  
  // DPD interactions (FD)
  Array<std::string> DPDtypeA; ///< First elements of the couples of particles types that interact by DPD FD forces
  Array<std::string> DPDtypeB; ///< Second elements of the couples of particles types that interact by DPD FD forces

  // DPDE interactions (FD)
  Array<std::string> DPDEtypeA; ///< First elements of the couples of particles types that interact by DPDE FD forces
  Array<std::string> DPDEtypeB; ///< Second elements of the couples of particles types that interact by DPDE FD forces

  // SDPD interactions (P+FD)
  Array<std::string> SDPDtypeA; ///< First elements of the couples of particles types that interact by SDPD forces
  Array<std::string> SDPDtypeB; ///< Second elements of the couples of particles types that interact by SDPD forces

  // SDPD_WALL interactions (LJr + P)
  Array<std::string> SDPDWalltypeReal; ///< Type of the real smooth particle interacting with the wall
  Array<std::string> SDPDWalltypeVirtual; ///< Type of the virtual wall particles interacting with real smooth particles

  // DPDE EOS
  Array<std::string> DPDEtype; ///< Particle types using the DPDE equation of state
  Array<double> DPDEcv; ///< Cv for the DPDE equations of state

  // Ideal Gas EOS
  Array<std::string> IGtype; ///< Particle types using the ideal gas equation of state

  // Mie-Gruneisen EOS
  Array<std::string> MGtype; ///< Particle types using the Mie-Gruneisen equation of state
  Array<double> MGgamma0; ///< \f[ \Gamma_0 \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGgammaInf; ///< \f[ \Gamma_{\infty} \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGtheta0; ///< \f[ \theta_0 \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGq; ///< \f[ q \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGrho0; ///< \f[ \rho_0 \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGks; ///< \f[ K_S \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGns; ///< \f[ N_S \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGrhos; ///< \f[ \rho_S \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGur; ///< \f[ u_r \f] parameter for the Mie-Gruneisen equation of state
  Array<double> MGcvr; ///< \f[ Cv_r \f] parameter for the Mie-Gruneisen equation of state

  // HZ EOS
  Array<std::string> HZtype; ///< Particle types using the HZ equation of state
  Array<double> HZgamma0; ///< Gruneisen parameter for the HZ equation of state
  Array<double> HZrho0; ///< Reference density for the HZ equation of state
  Array<double> HZc0; ///< \f[ c_0 \f] parameter for the HZ equation of state
  Array<double> HZcv; ///< Heat capacity for the HZ equation of state
  Array<double> HZs; ///< \f[ s \f] parameter for the HZ equation of state

  // JWL EOS
  Array<std::string> JWLtype; ///< Particle types using the JWL equation of state
  Array<double> JWLgamma0; ///< Gruneisen parameter for the JWL equation of state
  Array<double> JWLrho0; ///< Reference density for the JWL equation of state
  Array<double> JWLe0; ///< Reference energy for the JWL equation of state
  Array<double> JWLdcj; ///< Detonation velocity
  Array<double> JWLpcj; ///< Pressure at the CJ point
  Array<double> JWLtcj; ///< Temperature at the CJ point
  Array<double> JWLcv; ///< Heat capacity for the JWL equation of state
  Array<double> JWLa; ///< \f[ a \f] parameter for the JWL equation of state
  Array<double> JWLb; ///< \f[ b \f] parameter for the JWL equation of state
  Array<double> JWLr1; ///< \f[ R_1 \f] parameter for the JWL equation of state
  Array<double> JWLr2; ///< \f[ R_2 \f] parameter for the JWL equation of state

  // Reactive EOS
  Array<std::string> Reactivetype; ///< Particle types using a reactive equation of state
  Array<std::string> ReactiveEOS0; ///< Type of the first equation of state
  Array<std::string> ReactiveEOS1; ///< Type of the second equation of state
  
  // Chemistry
  // Second order reaction
  Array<std::string> SOtype; ///< Particle types reacting with a second order kinetics
  Array<double> SOzab; ///< Arrhenius prefactor for direct second order reaction (A->B)
  Array<double> SOzba; ///< Arrhenius prefactor for reverse second order reaction (A<-B)
  Array<double> SOeab; ///< Activation energy for direct second order reaction (A->B)
  Array<double> SOeba; ///< Activation energy for reverse second order reaction (A<-B)
  
  // Particles initialization data
  InitType m_initType; ///< Type of particles initialization
  std::string initFile; ///< Initialization file (in the case of an initialization from a dump file)
  std::string initDir; ///< Initialization file directory (case of an initialization from a dump file)
  int initStep; ///< Initial step
  LatticeType latticeType; ///< Type of lattice (case of an initialization from a lattice)
  double latticeParameter; ///< Space parameter of the lattice (case of an initialization from a lattice)
  Array<std::string> latticeAtoms; ///< Lists of the atoms in a cell of the lattice (case of an initialization from a lattice)
  vec3<int> numberOfCells; ///< Number of cells in the lattice (case of an initialization from a lattice)
  double initialTemperature; ///< Initial temperature
  double initialTint; ///< Initial internal temperature

  double wallWidth; ///< Width of the wall
  vec3<std::string> wallLowerTypes; ///< Types of particle constituting the lower walls
  vec3<std::string> wallUpperTypes; ///< Types of particle constituting the upper walls
  
  LegacyHeaderIOStruct* legacyHeader; ///< Header of the initialization file (case of an initialization from a dump file)

  // Molecules simulations data
  ForceField m_forceField; 			///< Force field used to compute the interactions
  ChargeMethod m_chargeMethod;	///< Method used to get the partial charges
  std::string m_format;					///< Format of the molecule input file
  bool m_addHydrogens; 					///< Indicate if hydrogens must be added into the molecule input file
  int m_numConformer;						///< Conformer to read in the molecule input file
  vec3<double> m_margins;				///< Minimum margins between the molecule and the border of the system
  double m_maxBonds;						///< Maximum bond length


  // Adaptive Mesh Refinement
  int Dmax;
  int amrCriterion;

  // ctest
  bool checkValues;
  double checkEpsilon;
  int checkIteration;
  std::string checkFile;
};


/// @brief Base structure to extract some information from an Input structure in
/// order to initialize another class (see all specializations)
template <class T> class Configuration {

public:

	/// @brief Default constructor
  Configuration();
  /// @brief Destructor (nothing to do)
  ~Configuration();

  /// @brief Constructor from the input
	/// @param [in] input Input data
  Configuration(const Input& input);

};


void readFromFile(const std::string& filename, std::string& datastream, const char comment);


class CommManager;


void readFromString(std::string& datastream, Input* input);

#endif // __INPUT_HPP_INCLUDED
