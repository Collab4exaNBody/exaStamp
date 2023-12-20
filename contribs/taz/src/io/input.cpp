/// @file 
/// @brief Tools to read a input file and get the corresponding Input structure



#include <fstream>

#include "io/input.hpp"
#include "io/particleInput.hpp"

#include "utils/stampUnits.hpp"
#include "utils/stringUtils.hpp"
#include "omp.h"


template <> Input::NodeType from_string(const std::string& name,const std::string& str);
template <> Input::SimulationMode from_string(const std::string& name,const std::string& str);
template <> Input::DecompositionType from_string(const std::string& name,const std::string& str);
template <> Input::BoundaryCondition from_string(const std::string& name,const std::string& str);
template <> Input::IntegrationScheme from_string(const std::string& name,const std::string& str);
template <> Input::BalancingMethod from_string(const std::string& name,const std::string& str);
template <> InitType from_string(const std::string& name,const std::string& str);
template <> Input::LatticeType from_string(const std::string& name,const std::string& str);
template <> DumpType from_string(const std::string& field,const std::string& str);
template <> OutputType from_string(const std::string& field,const std::string& str);
template <> Input::SynchronousPolicy from_string(const std::string& name,const std::string& str);
template <> Input::ForceField from_string(const std::string& field,const std::string& str);
template <> Input::ChargeMethod from_string(const std::string& field,const std::string& str);
template <> Input::NeighMethod from_string(const std::string& name,const std::string& str);

// AMR
template <> omp_sched_t from_string(const std::string& name,const std::string& str);

/// @brief Copy a file into a datastream
/// @param [in] filename Filename
/// @param [out] datastream String that will temporarily store the input data
/// @param [in] comment Comment delimiter
void readFromFile(const std::string& filename, std::string& datastream, const char comment) {

  // Open file
  std::ifstream in;
  in.open(filename.c_str());

  // Opening test
  if (!in.is_open()) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'readFromFile(const std::string&, std::string&, const char)' : unable to open file " << filename << ". STOP." 
	     << std::endl;
    exit(-1);
  }

  // Read file
  std::string line, field, value;
  while (in.good()) {

    // Get line, remove spaces and comments
    getline(in, line);
    trim_spaces(line);
    trim_comment(line, comment);

    // Add endline character and add to the datastream
    if (!line.empty()) {
      line += "\n";
      datastream.append(line);
    }
  }

  // Close file
  in.close();

}


/// @brief Read program input from a datastream
/// @param [in,out] datastream String that temporarily store the input data
/// @param [in,out] input Input data structure
void readFromString(std::string& datastream, Input* input) {

	// Full env scan warning
  if (getenv("XSP_FORCE_SCAN_ENV")){
    std::cerr << "    env XSP_FORCE_SCAN_ENV" << std::endl;
    std::cerr << "    You need to define each field in input_file, if you want replace value per the value of environment variable." << std::endl;
  }
  std::string line, field, value;
  // While there are lines in the datastream
  while (!datastream.empty()) {
    // Get line, remove remaining spaces and comments
    auto position = datastream.find_first_of("\n");
    line = datastream.substr(0, position);
    trim_spaces(line);
    trim_comment(line, input->comment); // This was done previously and therefore seems useless

    if (!line.empty()) {
    	// Extract the field name and value from the line
      extract_content(line, field, input->assign, value);
      // Assign the value to the field
      input->convertToType(field, value);
    }
    // Delete read line
    datastream.erase(0, position+1);
  }

  // Set the origin and extension from the header in the case of an initialization from a Stamp dump file
  if (input->m_initType==InitType::INIT_STAMP_LEGACY_DUMP) {
    // TODO : a deleger au IOManager

    // Read the header from the initialization file
    input->legacyHeader = new LegacyHeaderIOStruct();
    readLegacyHeader(input->legacyHeader, input->initFile, input->initStep, input->inputOutputManager);

    // Get the minimum coordinates as the origin and set the extension from the maximum coordinates
    input->origin.x = input->legacyHeader->xmin;
    input->origin.y = input->legacyHeader->ymin;
    input->origin.z = input->legacyHeader->zmin;
    vec3<double> max(input->legacyHeader->xmax, input->legacyHeader->ymax, input->legacyHeader->zmax);
    input->extension = max-input->origin;
  }
  // Set the origin and extension from the header in the case of an initialization from a Hercule dump file
  else if (input->m_initType==InitType::INIT_HERCULE_DUMP) {
#ifdef __use_lib_hercule
  	// Read the header from the initialization file
    InputOutputManager *ioMgr =  input->inputOutputManager;
    input->legacyHeader = new LegacyHeaderIOStruct();
    ioMgr->herculeDumpR_readHeader(input->legacyHeader, input->initFile, input->initStep, 0, input->inputOutputManager);

    // Get the minimum coordinates as the origin and set the extension from the maximum coordinates
    input->origin.x = input->legacyHeader->xmin;
    input->origin.y = input->legacyHeader->ymin;
    input->origin.z = input->legacyHeader->zmin;
    vec3<double> max(input->legacyHeader->xmax, input->legacyHeader->ymax, input->legacyHeader->zmax);
    input->extension = max-input->origin;
#else
    std::cerr << "Can't use Hercule input." << std::endl;

#endif
  }
  // Correct the number of lattice cells to be positive
  input->numberOfCells=auxMax(input->numberOfCells,vec3<int>(1));
  // If the extension is unsatisfactory, get it from the lattice size except if the
  // initialization is done from a stampv4 file, in which case it will be done latter
  // The extension is increased by 1e-15 to avoid rounding error (number of cells in the lattice initialization)
  if (input->m_initType!=InitType::INIT_STAMPV4 && (input->extension.x<=0. || input->extension.y<=0. || input->extension.z<=0.)) {
    input->extension = input->numberOfCells*input->latticeParameter+1e-15;
  }
  // Correct the margins to be positive
  input->m_margins=auxMax(input->m_margins,vec3<double>(2e-10));
#ifndef __use_lib_openbabel
  if(input->m_initType==InitType::INIT_STAMPV4) {
  	std::cerr << "Error in readFromString :" << std::endl;
  	std::cerr << "The stampv4 initialization cannot be used without the OpenBabel library (NICOLAS: WHY?)." << std::endl;
  	exit(-1);
  }
  if(input->m_chargeMethod==Input::ChargeMethod::OPENBABEL) {
  	std::cerr << "Error in readFromString :" << std::endl;
  	std::cerr << "Default charges from OpenBabel cannot be used without the OpenBabel library." << std::endl;
  	exit(-1);
  }
#endif


  // Convert the values into Stamp units
  input->convertToStampUnits();

}


/// @brief Constructor
///
/// Set default values to all the fields
/// @param [in] _inputOutputManager Input/Output manager
Input::Input(InputOutputManager* _inputOutputManager)
  : inputOutputManager(_inputOutputManager),
    version(0), 
    nodeType(from_string<Input::NodeType>("nodeType","default")),
    mode(from_string<Input::SimulationMode>("mode","default")), 
    scheme(from_string<Input::IntegrationScheme>("scheme","default")), 
    frictionLangevin(0.), temperatureThermostat(0.),
    seed(0),
    recovering(false), symmetrization(true),
    maxNumberOfThreadsPerNode(-1),
    bindThreads(false),
    decomposition(from_string<Input::DecompositionType>("decomposition","default")), 
    decoupage(1),
    boundaryConditions(PERIODIC),
    origin(0.), extension(-1.),
    expandFreeLimits(false), expandFactor(1.1),
    numberOfSteps(0), delta(1.0e-15), stopWallsTime(-1.), logRate(1), outputRate(-1), outputType(VTK), outputCells(false),
    dumpRate(-1), dumpType(LEGACY_DUMP),
    dir("."), dumpDir(), outputDir(),
    syncPolicy(from_string<Input::SynchronousPolicy>("syncPolicy","default")),
    analyticsGraph(0), 
    isTimers(false),isTimersPotential(false),isTimersNeighbours(false),isTimersGhost(false),isTimersRefine(false),isTimersGeneric(false),
    doTypes(true), doVelocities(true), doEints(false),
    maxNbAnalyticsWorkers(1),
    maxNbSimulationWorkers(-1),
    balancingMethod(NONE), balanceRate(-1),
    neighMethod(VERLET_LIST), rVerlet(0),
    logWorkload(false),
    atomNames(0), 
    atomAtomicNumbers(0), atomMasses(0), atomCharges(0.),
    IGtypeA(0), IGtypeB(0), 
    LJtypeA(0), LJtypeB(0),
    LJrcut(0), LJepsilon(0), LJsigma(0), 
    Exp6typeA(0), Exp6typeB(0),
    Exp6rcut(0), Exp6foo(0),
    SCtypeA(0), SCtypeB(0),
    SCc(0), SCepsilon(0), SCa0(0), SCn(0), SCm(0),
    EamVNIITFtypeA(0), EamVNIITFtypeB(0),
    EamVNIITFrcut(0), EamVNIITFrmax(0), EamVNIITFrmin(0), EamVNIITFrt0(0),
    EamVNIITFEcoh(0), EamVNIITFE0(0), EamVNIITFbeta(0), EamVNIITFA(0), EamVNIITFZ(0), EamVNIITFn(0),
    EamVNIITFalpha(0), EamVNIITFD(0), EamVNIITFeta(0), EamVNIITFmu(0),
    MeamTypeA(0), MeamTypeB(0),
    MeamRcut(0), MeamRmax(0), MeamRmin(0),
    MeamEcoh(0),MeamE0(0), MeamA(0), MeamR0(0),MeamAlpha(0),MeamDelta(0),
    MeamBeta0(0), MeamBeta1(0), MeamBeta2(0), MeamBeta3(0),
    MeamT0(0),MeamT1(0),MeamT2(0),MeamT3(0),
    MeamS0(0),MeamS1(0),MeamS2(0),MeamS3(0),
    MeamCmin(0),MeamCmax(0), MeamZ(0),
    MeamRc(0), MeamRp(0),
    GaussTypeA(0), GaussTypeB(0),
    GaussRcut(0), GaussEpsilon(0), GaussRatt(0), GaussRrep(0), GaussRatio(0),
    mesoNames(0), 
    mesoMasses(0), mesoUnitMasses(0.), mesoGammaPara(0.), mesoGammaOrtho(0.),
    smoothNames(0), 
    smoothMasses(0), smoothUnitMasses(0), smoothBulkViscosity(0), smoothShearViscosity(0), smoothSmoothingLength(0),smoothKernel(0),
    wallNames(0), 
    wallMasses(0), wallUnitMasses(0), wallSmoothingLength(0),wallKernel(0),wallVelocityX(0),wallVelocityY(0),wallVelocityZ(0),
    DPDtypeA(0.), DPDtypeB(0.),
    DPDEtypeA(0.), DPDEtypeB(0.),
    SDPDtypeA(0.), SDPDtypeB(0.),
    SDPDWalltypeReal(0.), SDPDWalltypeVirtual(0.),
    DPDEtype(0), DPDEcv(0),
    IGtype(0),
    MGtype(0),
    MGgamma0(0), MGgammaInf(0), MGtheta0(0), MGq(0), MGrho0(0),
    MGks(0), MGns(0), MGrhos(0), MGur(0), MGcvr(0),
    HZtype(0),
    HZgamma0(0), HZrho0(0), HZc0(0), HZcv(0), HZs(0),
    JWLtype(0),
    JWLgamma0(0), JWLrho0(0), JWLe0(0), JWLdcj(0), JWLpcj(0), JWLtcj(0),
    JWLcv(0), JWLa(0), JWLb(0), JWLr1(0), JWLr2(0),
    Reactivetype(0), ReactiveEOS0(0), ReactiveEOS1(0),
    SOtype(0),
    SOzab(0), SOzba(0), SOeab(0), SOeba(0),
    m_initType(from_string<InitType>("particlesInitialization","default")),
    initFile(""),initDir("."),initStep(-1),
    latticeType(from_string<Input::LatticeType>("latticeType","default")), latticeParameter(0.367e-09),
    latticeAtoms(4, ""), numberOfCells(0), initialTemperature(300.), initialTint(0.),
    wallWidth(0), wallLowerTypes(""), wallUpperTypes(""),
    legacyHeader(nullptr),
    m_forceField(from_string<Input::ForceField>("m_forceField","default")),
    m_chargeMethod(Input::NOT_SPECIFIED),
    m_format(""),
    m_addHydrogens(false),
    m_numConformer(1),
    m_margins(0.),
    m_maxBonds(1.05),
    Dmax(0),
    amrCriterion(0),
    checkValues(false),
    checkEpsilon(0),
    checkIteration(0),
    checkFile("default")
    {}


/// @brief Aborts program if a field name is not recognized
/// @param [in] field Field name
/// @param [in] _trace Indicates if a warning must be printed, default=true
void error_field_name(const std::string& field, bool _trace=true)
{
  if (_trace){
    std::cerr<< "Error: unrecognized field '" << field << "'" << std::endl;
  }
  MPI_Abort(MPI_COMM_WORLD,2);
}


/// @brief Aborts program if a value is not recognized
/// @param [in] field Field name
/// @param [in] val Value
/// @param [in] _trace Indicates if a warning must be printed, default=true
void error_field_val(const std::string& field,const std::string& val, bool _trace=true)
{
  if (_trace){
    std::cerr<< "Error: unrecognized value '" << val << "' of field '" << field << "'" << std::endl;
  }
  MPI_Abort(MPI_COMM_WORLD,2);
}


/// @brief Converts the string value to an adequate value and assigns this value to the field
/// @param [in] field Field name
/// @param [in] val String value
// It may be possible to do this in a  more elegant way, but I had to go fast ...
void Input::convertToType(const std::string& field, const std::string& val) {

	// Shortcuts for the delimiters
  const char& ld =Input:: array_delimiter_left;
  const char& rd = Input::array_delimiter_right;
  const char& s = Input::array_separator;

  std::string value = val;

  // Check if the setup from environment variables is activated
  bool force_scan_env = getenv("XSP_FORCE_SCAN_ENV") != NULL;

  // If setup from environment variable for this field or all field
  if (is_equal(value, env_flag) || force_scan_env) {

  	// Lambda function to get a value from an environment variable
  	// Parameters :	[in] 			_force_scan_env	Indicates if the scan originate from a full environment scan
  	//							[in,out]	str							Value to get
  	//							[in]			env							Name of the environment variable
    auto read_from_env = [&] (bool _force_scan_env, std::string& str, const std::string& env) -> void {
    	// Get the environment variable
      char* tmp = getenv(env.c_str());
      // If there is no environment variable but it's a full scan, it's ok, we'll just assume the value is already defined in the file
      if (_force_scan_env && tmp == NULL) return;
      // Put the environment variable into the value and print to check
      str = std::string(tmp);
      if (tmp != NULL && this->inputOutputManager->isMaster()){
      	std::cerr << "    env " << env << " = " << str << std::endl;
      }
    };

  	// Lambda function to get a 3D value from three environment variables
  	// Parameters :	[in] 			_force_scan_env	Indicates if the scan originate from a full environment scan
  	//							[in,out]	str							Value to get
  	//							[in]			envX						Name of the environment variable for the x component
  	//							[in]			envY						Name of the environment variable for the y component
  	//							[in]			envZ						Name of the environment variable for the z component
    auto read_from_env_3 = [&] (bool _force_scan_env,std::string& str, const std::string& envX, const std::string& envY, const std::string& envZ) -> void {
    	// Get the environment variables
      char* tmpx = getenv(envX.c_str());
      char* tmpy = getenv(envY.c_str());
      char* tmpz = getenv(envZ.c_str());
      // If there is no environment variable but it's a full scan, it's ok, we'll just assume the value is already defined in the file
      if (_force_scan_env && (tmpx == NULL || tmpy == NULL || tmpz == NULL )) return;
      // Put the environment variables into the value (as an array) and print to check
      str = ld + std::string(tmpx) + s + std::string(tmpy) + s + std::string(tmpz) + rd;
      if ((tmpx != NULL && tmpy != NULL && tmpz != NULL) && this->inputOutputManager->isMaster()){
      	std::cerr << "env [" << envX << "," << envY << "," << envZ << "] = " << s << std::endl;
      }
    };

    // Scan the right environment variable to get a new string value
    if      (is_equal(field, "recovering"          )) read_from_env  (force_scan_env,value, env_recover);
    else if (is_equal(field, "symmetrization"      )) read_from_env  (force_scan_env,value, env_symm);
    // The following line allows spelling mistake for french speaking people
    else if (is_equal(field, "symetrization"       )) read_from_env  (force_scan_env,value, env_sym);
    else if (is_equal(field, "decoupage"           )) read_from_env_3(force_scan_env,value, env_np_x, env_np_y, env_np_z);
    else if (is_equal(field, "max_threads_per_node")) read_from_env  (force_scan_env,value, env_nthreads);
    else if (is_equal(field, "decomposition"       )) read_from_env  (force_scan_env,value, env_dcmp);
    else if (is_equal(field, "dynamic_balancing"   )) read_from_env  (force_scan_env,value, env_lb_mthd);
    else if (is_equal(field, "log_rate"            )) read_from_env  (force_scan_env,value, env_log_rate);
    else if (is_equal(field, "balance_rate"        )) read_from_env  (force_scan_env,value, env_lb_rate);
    else if (is_equal(field, "output_type"         )) read_from_env  (force_scan_env,value, env_out_type);
    else if (is_equal(field, "output_rate"         )) read_from_env  (force_scan_env,value, env_out_rate);
    else if (is_equal(field, "output_dir"          )) read_from_env  (force_scan_env,value, env_out_dir);
    else if (is_equal(field, "dump_type"           )) read_from_env  (force_scan_env,value, env_dump_type);
    else if (is_equal(field, "dump_rate"           )) read_from_env  (force_scan_env,value, env_dump_rate);
    else if (is_equal(field, "dump_dir"            )) read_from_env  (force_scan_env,value, env_dump_dir);
    else if (is_equal(field, "initialization"      )) read_from_env  (force_scan_env,value, env_initialization);
    else if (is_equal(field, "init_step"           )) read_from_env  (force_scan_env,value, env_init_step);
    else if (is_equal(field, "particles_file"      )) read_from_env  (force_scan_env,value, env_particles_file);
    else if (is_equal(field, "number_of_steps"     )) read_from_env  (force_scan_env,value, env_n_steps);
    else if (is_equal(field, "n_cells"             )) read_from_env_3(force_scan_env,value, env_nc_x, env_nc_y, env_nc_z);
  }

  // Convert the string value to an adequate value depending on the field
  if      (is_equal(field, "version")) version = from_string<int>(field,value);

  else if (is_equal(field, "node_type")) nodeType = from_string<NodeType>(field,value);

  else if (is_equal(field, "mode"                  )) mode                  = from_string<SimulationMode>(field,value);
  else if (is_equal(field, "scheme"                )) scheme                = from_string<IntegrationScheme>(field,value);
  else if (is_equal(field, "friction_langevin"     )) frictionLangevin      = from_string<double>(field,value);
  else if (is_equal(field, "temperature_thermostat")) temperatureThermostat = from_string<double>(field,value);
  else if (is_equal(field, "temperature_langevin"  )) temperatureThermostat = from_string<double>(field,value); // for compatibility
  else if (is_equal(field, "temperature_dpd"       )) temperatureThermostat = from_string<double>(field,value); // for compatibility
  else if (is_equal(field, "seed"                  )) seed                  = from_string<int>(field,value);

  else if (is_equal(field, "recovering"          )) recovering                = from_string<bool>(field,value);
  else if (is_equal(field, "symmetrization"      )) symmetrization            = from_string<bool>(field,value);
  // The following line allows spelling mistake for french speaking people
  else if (is_equal(field, "symetrization"       )) symmetrization            = from_string<bool>(field,value);
  else if (is_equal(field, "decoupage"           )) decoupage                 = from_string_vec3<int>(field,value, ld, rd, s);
  else if (is_equal(field, "max_threads_per_node")) maxNumberOfThreadsPerNode = from_string<int>(field,value);
  else if (is_equal(field, "bind_threads"        )) bindThreads               = from_string<bool>(field,value);

  else if (is_equal(field, "dynamic_balancing")) balancingMethod = from_string<BalancingMethod>(field,value);
  else if (is_equal(field, "neighbours_method"   )) neighMethod               = from_string<NeighMethod>(field,value);
  else if (is_equal(field, "raduis_verlet"      )) rVerlet                   = from_string<double>(field,value);
  else if (is_equal(field, "balance_rate"     )) balanceRate     = from_string<int>(field,value);
  else if (is_equal(field, "log_workload"     )) logWorkload     = from_string<bool>(field,value);

  else if (is_equal(field, "boundary_conditions")) boundaryConditions = from_string_vec3<BoundaryCondition>(field,value, ld, rd, s);
  else if (is_equal(field, "decomposition"      )) decomposition      = from_string<DecompositionType>(field,value);

  else if (is_equal(field, "origin"             )) origin           = from_string_vec3<double>(field,value, ld, rd, s);
  else if (is_equal(field, "extension"          )) extension        = from_string_vec3<double>(field,value, ld, rd, s);
  else if (is_equal(field, "expand_free_limits" )) expandFreeLimits = from_string<bool>(field,value);
  else if (is_equal(field, "expand_factor"      )) expandFactor     = from_string<double>(field,value);

  else if (is_equal(field, "number_of_steps"           )) numberOfSteps          = from_string<int>(field,value);
  else if (is_equal(field, "delta"                     )) delta                  = from_string<double>(field,value);
  else if (is_equal(field, "wall_stop_time"            )) stopWallsTime          = from_string<double>(field,value);
  else if (is_equal(field, "log_rate"                  )) logRate                = from_string<int>(field,value);
  else if (is_equal(field, "output_rate"               )) outputRate             = from_string<int>(field,value);
  else if (is_equal(field, "output_type"               )) outputType             = from_string<OutputType>(field,value);
  else if (is_equal(field, "output_cells"              )) outputCells            = from_string<bool>(field,value);
  else if (is_equal(field, "output_dir"                )) outputDir              = value;
  else if (is_equal(field, "dump_rate"                 )) dumpRate               = from_string<int>(field,value);
  else if (is_equal(field, "dump_type"                 )) dumpType               = from_string<DumpType>(field,value);
  else if (is_equal(field, "dump_dir"                  )) dumpDir                = value;
  else if (is_equal(field, "dir"                       )) dir                    = value;
  else if (is_equal(field, "sync_policy"               )) syncPolicy             = from_string<SynchronousPolicy>(field,value);
  else if (is_equal(field, "analytics_graph"           )) analyticsGraph         = from_string_array<std::string>(field, value, ld, rd, s);
  else if (is_equal(field, "doTypes"                   )) doTypes                = from_string<bool>(field, value);
  else if (is_equal(field, "doVelocities"              )) doVelocities           = from_string<bool>(field, value);
  else if (is_equal(field, "doEints"                   )) doEints                = from_string<bool>(field, value);
  else if (is_equal(field, "max_threads_for_analytics" )) maxNbAnalyticsWorkers  = from_string<int>(field,value);
  else if (is_equal(field, "max_threads_for_simulation")) maxNbSimulationWorkers = from_string<int>(field,value);
  else if (is_equal(field, "input_dir"        	       )) initDir                = value;

  else if (is_equal(field, "Chrono"                    )) isTimers               = from_string<bool>(field,value);
  else if (is_equal(field, "Chrono_potential"          )) isTimersPotential      = from_string<bool>(field,value);
  else if (is_equal(field, "Chrono_neighbours"         )) isTimersNeighbours     = from_string<bool>(field,value);
  else if (is_equal(field, "Chrono_ghost"              )) isTimersGhost          = from_string<bool>(field,value);
  else if (is_equal(field, "Chrono_refine"             )) isTimersRefine         = from_string<bool>(field,value);
  else if (is_equal(field, "Chrono_generic"            )) isTimersGeneric        = from_string<bool>(field,value);
  
  else if (is_equal(field, "atom_names"  )) atomNames         = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "atom_z"      )) atomAtomicNumbers = from_string_array<int>(field,value, ld, rd, s);
  else if (is_equal(field, "atom_masses" )) atomMasses        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "atom_charges")) atomCharges       = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "ig_type_A")) IGtypeA = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "ig_type_B")) IGtypeB = from_string_array<std::string>(field,value, ld, rd, s);

  else if (is_equal(field, "lj_type_A" )) LJtypeA   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "lj_type_B" )) LJtypeB   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "lj_rcut"   )) LJrcut    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "lj_epsilon")) LJepsilon = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "lj_sigma"  )) LJsigma   = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "sutton-chen_type_A" )) SCtypeA   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "sutton-chen_type_B" )) SCtypeB   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "sutton-chen_rcut"   )) SCrcut    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "sutton-chen_c"      )) SCc       = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "sutton-chen_epsilon")) SCepsilon = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "sutton-chen_a0"     )) SCa0      = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "sutton-chen_n"      )) SCn       = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "sutton-chen_m"      )) SCm       = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "eam-vniitf_type_A" )) EamVNIITFtypeA = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_type_B" )) EamVNIITFtypeB = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_rcut"   )) EamVNIITFrcut  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_rmax"   )) EamVNIITFrmax  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_rmin"   )) EamVNIITFrmin  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_rt0"    )) EamVNIITFrt0   = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_ecoh"   )) EamVNIITFEcoh  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_e0"     )) EamVNIITFE0    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_beta"   )) EamVNIITFbeta  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_A"      )) EamVNIITFA     = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_Z"      )) EamVNIITFZ     = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_n"      )) EamVNIITFn     = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_alpha"  )) EamVNIITFalpha = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_D"      )) EamVNIITFD     = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_eta"    )) EamVNIITFeta   = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "eam-vniitf_mu"     )) EamVNIITFmu    = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "meam_type_A" )) MeamTypeA = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_type_B" )) MeamTypeB = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_rcut"   )) MeamRcut  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_rmax"   )) MeamRmax  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_rmin"   )) MeamRmin  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_ecoh"   )) MeamEcoh  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_e0"     )) MeamE0    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_A"      )) MeamA     = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_r0"     )) MeamR0   = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_alpha"  )) MeamAlpha = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_delta"  )) MeamDelta = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_beta0"  )) MeamBeta0 = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_beta1"  )) MeamBeta1 = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_beta2"  )) MeamBeta2 = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_beta3"  )) MeamBeta3 = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_t0"     )) MeamT0    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_t1"     )) MeamT1    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_t2"     )) MeamT2    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_t3"     )) MeamT3    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_s0"     )) MeamS0    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_s1"     )) MeamS1    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_s2"     )) MeamS2    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_s3"     )) MeamS3    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_Cmin"   )) MeamCmin  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_Cmax"   )) MeamCmax  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_Z"      )) MeamZ     = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_rc"     )) MeamRc    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meam_rp"     )) MeamRp    = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "gauss_type_A" )) GaussTypeA   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "gauss_type_B" )) GaussTypeB   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "gauss_rcut"   )) GaussRcut    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "gauss_epsilon")) GaussEpsilon = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "gauss_r_att"  )) GaussRatt    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "gauss_r_rep"  )) GaussRrep    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "gauss_ratio"  )) GaussRatio   = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "meso_names"      )) mesoNames         = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "meso_masses"     )) mesoMasses        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meso_sizes"      )) mesoUnitMasses    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meso_gamma_para" )) mesoGammaPara     = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "meso_gamma_ortho")) mesoGammaOrtho    = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "smooth_names"           )) smoothNames              = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "smooth_masses"          )) smoothMasses             = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "smooth_unitmasses"      )) smoothUnitMasses         = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "smooth_bulk_viscosity"  )) smoothBulkViscosity      = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "smooth_shear_viscosity" )) smoothShearViscosity     = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "smooth_smoothing_length")) smoothSmoothingLength    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "smooth_kernel"          )) smoothKernel             = from_string_array<std::string>(field,value, ld, rd, s);

  else if (is_equal(field, "wall_names"           )) wallNames              = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "wall_masses"          )) wallMasses             = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "wall_unitmasses"      )) wallUnitMasses         = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "wall_smoothing_length")) wallSmoothingLength    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "wall_kernel"          )) wallKernel             = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "wall_velocity_x"      )) wallVelocityX          = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "wall_velocity_y"      )) wallVelocityY          = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "wall_velocity_z"      )) wallVelocityZ          = from_string_array<double>(field,value, ld, rd, s);
  
  else if (is_equal(field, "dpd_typeA"     )) DPDtypeA     = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "dpd_typeB"     )) DPDtypeB     = from_string_array<std::string>(field,value, ld, rd, s);

  else if (is_equal(field, "dpde_typeA"    )) DPDEtypeA    = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "dpde_typeB"    )) DPDEtypeB    = from_string_array<std::string>(field,value, ld, rd, s);

  else if (is_equal(field, "sdpd_typeA"    )) SDPDtypeA    = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "sdpd_typeB"    )) SDPDtypeB    = from_string_array<std::string>(field,value, ld, rd, s);

  else if (is_equal(field, "sdpd_wall_type_real"       )) SDPDWalltypeReal    = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "sdpd_wall_type_virtual"    )) SDPDWalltypeVirtual = from_string_array<std::string>(field,value, ld, rd, s);

  else if (is_equal(field, "dpde_type")) DPDEtype = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "dpde_cv"  )) DPDEcv   = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "ig_type")) IGtype = from_string_array<std::string>(field,value, ld, rd, s);

  else if (is_equal(field, "mg_type" )) MGtype      = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_g0"   )) MGgamma0    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_ginf" )) MGgammaInf  = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_t0"   )) MGtheta0    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_q"    )) MGq         = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_rho0" )) MGrho0      = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_ks"   )) MGks        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_ns"   )) MGns        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_rhos" )) MGrhos      = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_ur"   )) MGur        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "mg_cvr"  )) MGcvr       = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "hz_type" )) HZtype      = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "hz_g0"   )) HZgamma0    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "hz_rho0" )) HZrho0      = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "hz_c0"   )) HZc0        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "hz_cv"   )) HZcv        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "hz_s"    )) HZs         = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "jwl_type" )) JWLtype      = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_g0"   )) JWLgamma0    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_rho0" )) JWLrho0      = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_e0"   )) JWLe0        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_dcj"  )) JWLdcj       = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_pcj"  )) JWLpcj       = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_tcj"  )) JWLtcj       = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_cv"   )) JWLcv        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_a"    )) JWLa         = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_b"    )) JWLb         = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_r1"   )) JWLr1        = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "jwl_r2"   )) JWLr2        = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "reactive_type" )) Reactivetype   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "reactive_eos0" )) ReactiveEOS0   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "reactive_eos1" )) ReactiveEOS1   = from_string_array<std::string>(field,value, ld, rd, s);
  
  else if (is_equal(field, "reaction_so_type"               )) SOtype   = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "reaction_so_forward_prefactor"  )) SOzab    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "reaction_so_backward_prefactor" )) SOzba    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "reaction_so_forward_activation" )) SOeab    = from_string_array<double>(field,value, ld, rd, s);
  else if (is_equal(field, "reaction_so_backward_activation")) SOeba    = from_string_array<double>(field,value, ld, rd, s);

  else if (is_equal(field, "initialization"    )) m_initType = from_string<InitType>(field,value);
  else if (is_equal(field, "init_step"         )) initStep   = from_string<int>(field,value);
  else if (is_equal(field, "particles_file"    )) initFile                = from_string<std::string>(field,value);
  else if (is_equal(field, "lattice_type"      )) latticeType             = from_string<LatticeType>(field,value);
  else if (is_equal(field, "lattice_parameter" )) latticeParameter        = from_string<double>(field,value);
  else if (is_equal(field, "lattice_atoms"     )) latticeAtoms            = from_string_array<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "n_cells"           )) numberOfCells           = from_string_vec3<int>(field,value, ld, rd, s);
  else if (is_equal(field, "init_temperature"  )) initialTemperature      = from_string<double>(field,value);

  else if (is_equal(field, "force_field"))			m_forceField		= from_string<ForceField>(field,value);
  else if (is_equal(field, "charge_method"))		m_chargeMethod	= from_string<ChargeMethod>(field,value);
  else if (is_equal(field, "molecule_format"))	m_format 				= from_string<std::string>(field,value);
  else if (is_equal(field, "add_hydrogens")) 		m_addHydrogens	= from_string<bool>(field,value);
  else if (is_equal(field, "nb_conformer")) 		m_numConformer	= from_string<int>(field,value);
  else if (is_equal(field, "margins")) 					m_margins       = from_string_vec3<double>(field,value, ld, rd, s);
  else if (is_equal(field, "max_bond")) 				m_maxBonds      = from_string<double>(field,value);
 
  else if (is_equal(field, "init_tint"         )) initialTint             = from_string<double>(field,value);

  else if (is_equal(field, "wall_width"        )) wallWidth               = from_string<double>(field,value);
  else if (is_equal(field, "wall_lower_types"  )) wallLowerTypes          = from_string_vec3<std::string>(field,value, ld, rd, s);
  else if (is_equal(field, "wall_upper_types"  )) wallUpperTypes          = from_string_vec3<std::string>(field,value, ld, rd, s);
  
  // Adaptive Mesh Refinement
  else if (is_equal(field, "amr_octree_depth"  )) Dmax          = from_string<int>(field,value);  
  else if (is_equal(field, "amr_criterion"  )) amrCriterion          = from_string<int>(field,value);  

  // Ctest
  else if (is_equal(field, "check_values"  )) checkValues          = from_string<bool>(field,value);  
  else if (is_equal(field, "check_epsilon"  )) checkEpsilon          = from_string<double>(field,value);  
  else if (is_equal(field, "check_iteration"  )) checkIteration          = from_string<int>(field,value);  
  else if (is_equal(field, "check_file"  )) checkFile          = from_string<std::string>(field,value);  

  // If the field was not recognized, abort
  else{
    error_field_name(field);
  }

}


/// @brief Convert input data into local units for a better precision
void Input::convertToStampUnits() {

  // Lambda function to convert a 3D vector
	// Parameters :	[in,out]	vec					Vector to convert
	//							[in]			currentUnit	Old unit
	//							[in]			newUnit			New unit
  auto conv_vec3 = [](vec3<double>& vec, const Unit& currentUnit, const Unit& newUnit) -> void {
    vec.x = convert(vec.x, currentUnit, newUnit);
    vec.y = convert(vec.y, currentUnit, newUnit);
    vec.z = convert(vec.z, currentUnit, newUnit);
  };

  // Lambda function to convert an array
	// Parameters :	[in,out]	array				Array to convert
	//							[in]			currentUnit	Old unit
	//							[in]			newUnit			New unit
  auto conv_array = [](Array<double>& array, const Unit& currentUnit, const Unit& newUnit) -> void {
    for (auto& elem : array) elem = convert(elem, currentUnit, newUnit);
  };

  // Convert all the numerical values of the input
  conv_vec3(origin,    SI_Units_base::meter, Stamp_Units::length);
  conv_vec3(extension, SI_Units_base::meter, Stamp_Units::length);
  conv_vec3(m_margins, SI_Units_base::meter, Stamp_Units::length);

  delta = convert(delta, SI_Units_base::second, Stamp_Units::time);
  stopWallsTime = convert(stopWallsTime, SI_Units_base::second, Stamp_Units::time);
  
  frictionLangevin      = convert(frictionLangevin,    SI_Units_base::second.inv(), Stamp_Units::time.inv());
  temperatureThermostat = convert(temperatureThermostat, SI_Units_base::kelvin,       Stamp_Units::temperature);
  
  conv_array(atomMasses,  SI_Units_base::kilogram, Stamp_Units::mass);
  conv_array(atomCharges, SI_Units_base::coulomb,  Stamp_Units::charge);

  conv_array(LJrcut,    SI_Units_base::meter, Stamp_Units::length);
  conv_array(LJepsilon, SI_Units_base::joule, Stamp_Units::energy);
  conv_array(LJsigma,   SI_Units_base::meter, Stamp_Units::length);

  conv_array(SCrcut,    SI_Units_base::meter, Stamp_Units::length);
  conv_array(SCepsilon, SI_Units_base::joule, Stamp_Units::energy);
  conv_array(SCa0,      SI_Units_base::meter, Stamp_Units::length);

  conv_array(EamVNIITFrcut,      SI_Units_base::meter, Stamp_Units::length);
  conv_array(EamVNIITFrmax,      SI_Units_base::meter, Stamp_Units::length);
  conv_array(EamVNIITFrmin,      SI_Units_base::meter, Stamp_Units::length);
  conv_array(EamVNIITFrt0,       SI_Units_base::meter, Stamp_Units::length);
  conv_array(EamVNIITFEcoh,      SI_Units_base::joule, Stamp_Units::energy);
  conv_array(EamVNIITFE0,        SI_Units_base::joule, Stamp_Units::energy);

  conv_array(MeamRcut,      SI_Units_base::meter, Stamp_Units::length);
  conv_array(MeamRmax,      SI_Units_base::meter, Stamp_Units::length);
  conv_array(MeamRmin,      SI_Units_base::meter, Stamp_Units::length);
  conv_array(MeamR0,        SI_Units_base::meter, Stamp_Units::length);
  conv_array(MeamEcoh,      SI_Units_base::joule, Stamp_Units::energy);
  conv_array(MeamE0,        SI_Units_base::joule, Stamp_Units::energy);
  conv_array(MeamRc,        SI_Units_base::meter, Stamp_Units::length);
  conv_array(MeamRp,        SI_Units_base::meter, Stamp_Units::length);

  conv_array(GaussRcut,    SI_Units_base::meter, Stamp_Units::length);
  conv_array(GaussEpsilon, SI_Units_base::kilogram*SI_Units_base::second.inv()*SI_Units_base::second.inv(), Stamp_Units::mass*Stamp_Units::time.inv()*Stamp_Units::time.inv());
  conv_array(GaussRatt,    SI_Units_base::meter, Stamp_Units::length);
  conv_array(GaussRrep,    SI_Units_base::meter, Stamp_Units::length);


  
  conv_array(mesoMasses,     SI_Units_base::kilogram, Stamp_Units::mass);
  for (uint i = 0; i< mesoMasses.size(); i++)
    mesoUnitMasses[i] = mesoMasses[i]/mesoUnitMasses[i];
  conv_array(mesoGammaPara,  SI_Units_base::second.inv(), Stamp_Units::time.inv());
  conv_array(mesoGammaOrtho, SI_Units_base::second.inv(), Stamp_Units::time.inv());

  conv_array(smoothMasses,          SI_Units_base::kilogram, Stamp_Units::mass);
  conv_array(smoothUnitMasses,      SI_Units_base::kilogram, Stamp_Units::mass);
  conv_array(smoothBulkViscosity,   SI_Units_base::kilogram*SI_Units_base::meter.inv()*SI_Units_base::second.inv(), 
	     Stamp_Units::mass*Stamp_Units::length.inv()*Stamp_Units::time.inv());
  conv_array(smoothShearViscosity,  SI_Units_base::kilogram*SI_Units_base::meter.inv()*SI_Units_base::second.inv(), 
	     Stamp_Units::mass*Stamp_Units::length.inv()*Stamp_Units::time.inv());
  conv_array(smoothSmoothingLength, SI_Units_base::meter, Stamp_Units::length);

  conv_array(wallMasses,          SI_Units_base::kilogram, Stamp_Units::mass);
  conv_array(wallUnitMasses,      SI_Units_base::kilogram, Stamp_Units::mass);
  conv_array(wallSmoothingLength, SI_Units_base::meter, Stamp_Units::length);
  conv_array(wallVelocityX,       SI_Units_base::meter*SI_Units_base::second.inv(), Stamp_Units::length*Stamp_Units::time.inv());
  conv_array(wallVelocityY,       SI_Units_base::meter*SI_Units_base::second.inv(), Stamp_Units::length*Stamp_Units::time.inv());
  conv_array(wallVelocityZ,       SI_Units_base::meter*SI_Units_base::second.inv(), Stamp_Units::length*Stamp_Units::time.inv());
  
  conv_array(DPDEcv, SI_Units_base::joule/SI_Units_base::kelvin, Stamp_Units::energy/Stamp_Units::temperature);

  conv_array(MGtheta0, SI_Units_base::kelvin, Stamp_Units::temperature);
  conv_array(MGrho0, SI_Units_base::kilogram/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::mass/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length);
  conv_array(MGrhos, SI_Units_base::kilogram/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::mass/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length);
  conv_array(MGks, SI_Units_base::joule/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::energy/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length);
  conv_array(MGcvr, SI_Units_base::joule/SI_Units_base::kelvin/SI_Units_base::kilogram, Stamp_Units::energy/Stamp_Units::temperature/Stamp_Units::mass);

  conv_array(HZrho0, SI_Units_base::kilogram/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::mass/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length);
  conv_array(HZc0, SI_Units_base::meter/SI_Units_base::second, Stamp_Units::length/Stamp_Units::time);
  conv_array(HZcv, SI_Units_base::joule/SI_Units_base::kelvin/SI_Units_base::kilogram, Stamp_Units::energy/Stamp_Units::temperature/Stamp_Units::mass);

  conv_array(JWLrho0, SI_Units_base::kilogram/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::mass/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length);
  conv_array(JWLe0, SI_Units_base::joule/SI_Units_base::kilogram, Stamp_Units::energy/Stamp_Units::mass);
  conv_array(JWLdcj, SI_Units_base::meter/SI_Units_base::second, Stamp_Units::length/Stamp_Units::time);
  conv_array(JWLpcj, SI_Units_base::joule/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::energy/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length);
  conv_array(JWLtcj, SI_Units_base::kelvin, Stamp_Units::temperature);
  conv_array(JWLcv, SI_Units_base::joule/SI_Units_base::kelvin/SI_Units_base::kilogram, Stamp_Units::energy/Stamp_Units::temperature/Stamp_Units::mass);
  conv_array(JWLa, SI_Units_base::joule/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::energy/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length);
  conv_array(JWLb, SI_Units_base::joule/SI_Units_base::meter/SI_Units_base::meter/SI_Units_base::meter, Stamp_Units::energy/Stamp_Units::length/Stamp_Units::length/Stamp_Units::length);

  conv_array(SOzab, SI_Units_base::second.inv(), Stamp_Units::time.inv());
  conv_array(SOzba, SI_Units_base::second.inv(), Stamp_Units::time.inv());

  conv_array(SOeab, SI_Units_base::joule, Stamp_Units::energy);
  conv_array(SOeba, SI_Units_base::joule, Stamp_Units::energy);
  
  latticeParameter   = convert(latticeParameter,   SI_Units_base::meter,  Stamp_Units::length);
  initialTemperature = convert(initialTemperature, SI_Units_base::kelvin, Stamp_Units::temperature);
  initialTint        = convert(initialTint, SI_Units_base::kelvin, Stamp_Units::temperature);
  m_maxBonds=convert(m_maxBonds, SI_Units_base::meter, Stamp_Units::length);

  wallWidth = convert(wallWidth, SI_Units_base::meter, Stamp_Units::length);
  
  if (legacyHeader!=nullptr) {

    legacyHeader->xmin = convert(legacyHeader->xmin, SI_Units_base::meter,  Stamp_Units::length);
    legacyHeader->xmax = convert(legacyHeader->xmax, SI_Units_base::meter,  Stamp_Units::length);
    legacyHeader->ymin = convert(legacyHeader->ymin, SI_Units_base::meter,  Stamp_Units::length);
    legacyHeader->ymax = convert(legacyHeader->ymax, SI_Units_base::meter,  Stamp_Units::length);
    legacyHeader->zmin = convert(legacyHeader->zmin, SI_Units_base::meter,  Stamp_Units::length);
    legacyHeader->zmax = convert(legacyHeader->zmax, SI_Units_base::meter,  Stamp_Units::length);
    
    legacyHeader->time    = convert(legacyHeader->time   , SI_Units_base::second, Stamp_Units::time);
    legacyHeader->CPUtime = convert(legacyHeader->CPUtime, SI_Units_base::second, Stamp_Units::time);
    
  }

  if(neighMethod==Input::VERLET_LIST && !( scheme == Input::NEWTON_VERLET_VELOCITY ||  scheme == Input::LANGEVIN ) )
  {
    std::cout << " Verlet lists are not available with an other scheme than verlet velocity or langevin! " << std::endl;
    exit(0);
  }
  rVerlet = convert(rVerlet,      SI_Units_base::meter, Stamp_Units::length);

}


/// @version Specialization for a Input::NodeType enumeration type
template <> Input::NodeType from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "flat")) return Input::FLAT;
  else if (is_equal(str, "sotl")) return Input::SOTL;

  else return Input::FLAT;
}


/// @version Specialization for a Input::SimulationMode enumeration type
template <> Input::SimulationMode from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "default"))                    return Input::SINGLE_SPEC_ATOM;
  if      (is_equal(str, "single_spec_atom"))           return Input::SINGLE_SPEC_ATOM;
  else if (is_equal(str, "single_spec_meso"))           return Input::SINGLE_SPEC_MESO;
  else if (is_equal(str, "single_spec_smooth"))         return Input::SINGLE_SPEC_SMOOTH;
  else if (is_equal(str, "single_spec_stiff_molecule")) return Input::SINGLE_SPEC_STIFF_MOLEC;
  else if (is_equal(str, "atom_in_mol"))								return Input::ATOM_IN_MOL;

  error_field_val(field,str);

  return Input::SimulationMode();
}


/// @version Specialization for a Input::DecompositionType enumeration type
template <> Input::DecompositionType from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "default")) return Input::RECTILINEAR;
  if      (is_equal(str, "rectilinear")) return Input::RECTILINEAR;
  if      (is_equal(str, "any"))         return Input::ANY;

  error_field_val(field,str);

  return Input::RECTILINEAR;
}


/// @version Specialization for a Input::BoundaryCondition enumeration type
template <> Input::BoundaryCondition from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "free"))     return Input::FREE;
  else if (is_equal(str, "periodic")) return Input::PERIODIC;
  else if (is_equal(str, "wall"))     return Input::WALL;
  else if (is_equal(str, "free_wall")) return Input::FREE_WALL;
  else if (is_equal(str, "wall_free")) return Input::WALL_FREE;

  error_field_val(field,str);

  return Input::FREE;
}


/// @version Specialization for a Input::IntegrationScheme enumeration type
template <> Input::IntegrationScheme from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "default"))                return Input::NEWTON_VERLET_VELOCITY;
  if      (is_equal(str, "newton_verlet_leapfrog")) return Input::NEWTON_VERLET_LEAPFROG;
  else if (is_equal(str, "newton_verlet_velocity")) return Input::NEWTON_VERLET_VELOCITY;
  else if (is_equal(str, "langevin_splitting"    )) return Input::LANGEVIN;
  else if (is_equal(str, "dpd_splitting"         )) return Input::DPD;
  else if (is_equal(str, "dpde_ser"              )) return Input::DPDE_SER;
  else if (is_equal(str, "sdpd_vv_ser"           )) return Input::SDPD_VV_SER;
  else if (is_equal(str, "sdpd_langevin_ser"     )) return Input::SDPD_LANGEVIN_SER;

  error_field_val(field,str);

  return Input::NEWTON_VERLET_VELOCITY;
}


/// @version Specialization for a Input::BalancingMethod enumeration type
template <> Input::BalancingMethod from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "none"  )) return Input::NONE;
  else if (is_equal(str, "block" )) return Input::BLOCK;
  else if (is_equal(str, "random")) return Input::RANDOM;
  else if (is_equal(str, "coordinate_bisection")) return Input::COORD_BSC;
  else if (is_equal(str, "inertial_bisection"  )) return Input::INERT_BSC;
  else if (is_equal(str, "space_filling_curve" )) return Input::SPC_FILL_CRV;
  else if (is_equal(str, "phg"               )) return Input::PHG;
  else if (is_equal(str, "metis"               )) { 
#ifdef __use_lib_metis
		return Input::METIS;
#else
		return Input::NONE;
#endif
	}
  else if (is_equal(str, "scotch"               )) {
#ifdef __use_lib_scotch //not implemented
		return Input::SCOTCH;
#else
		return Input::NONE;
#endif
	}
  error_field_val(field,str);

  return Input::NONE;

}


/// @version Specialization for a InitType enumeration type
template <> InitType from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "default")) 						return InitType::INIT_DEFAULT;
  else if (is_equal(str, "legacy")) 						return InitType::INIT_STAMP_LEGACY_DUMP;
  else if (is_equal(str, "hercule")) 						return InitType::INIT_HERCULE_DUMP;
  else if (is_equal(str, "stampv4"))						return InitType::INIT_STAMPV4;

  error_field_val(field,str);

  return InitType::INIT_DEFAULT;
}


/// @version Specialization for a DumpType enumeration type
template <> DumpType from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "default")) return LEGACY_DUMP;
  else if (is_equal(str, "legacy")) return LEGACY_DUMP;
  else if (is_equal(str, "hercule")) return HERCULE_DUMP;

  error_field_val(field,str);

  return LEGACY_DUMP;
}


/// @version Specialization for a OutputType enumeration type
template <> OutputType from_string(const std::string& field,const std::string& str) {
  if      (is_equal(str, "default"     )) return VTK;
  else if (is_equal(str, "vtk"         )) return VTK;
  else if (is_equal(str, "dat"         )) return DAT;
  else if (is_equal(str, "hercule"     )) return HERCULE_PARTICULE_DEP;

  error_field_val(field,str);

  return VTK;
}


/// @version Specialization for a Input::LatticeType enumeration type
template <> Input::LatticeType from_string(const std::string& field,const std::string& str) {

  if (is_equal(str, "default" )) return Input::FCC;
  if (is_equal(str, "sc" )) return Input::SC;
  if (is_equal(str, "bcc" )) return Input::BCC;
  if (is_equal(str, "fcc" )) return Input::FCC;
  if (is_equal(str, "diam100" )) return Input::DIAM100;

  error_field_val(field,str);

  return Input::FCC;
}


/// @version Specialization for a Input::SynchronousPolicy enumeration type
template <> Input::SynchronousPolicy from_string(const std::string& field,const std::string& str) {

  if      (is_equal(str, "synchronous_without")) return Input::SYNCHRONOUS_WITHOUT_ARENAS;
  else if (is_equal(str, "synchronous_with")) return Input::SYNCHRONOUS_WITH_ARENAS;
  else if (is_equal(str, "asynchronous_without")) return Input::ASYNCHRONOUS_WITHOUT_ARENAS;
  else if (is_equal(str, "asynchronous_with")) return Input::ASYNCHRONOUS_WITH_ARENAS;

  else return Input::SYNCHRONOUS_WITHOUT_ARENAS;
}

/// @version Specialization for a Input::ForceField enumeration type
template <> Input::ForceField from_string(const std::string& field,const std::string& str) {

  if (is_equal(str, "default" )) 	return Input::UFF;
  if (is_equal(str, "uff" )) 			return Input::UFF;
  if (is_equal(str, "compass" )) 	return Input::COMPASS;
  if (is_equal(str, "amber"))			return Input::AMBER;

  error_field_val(field,str);

  return Input::UFF;
}


/// @version Specialization for a Input::ChargeMethod enumeration type
template <> Input::ChargeMethod from_string(const std::string& field,const std::string& str) {

  if (is_equal(str, "default" )) 		return Input::NOT_SPECIFIED;
  if (is_equal(str, "none" )) 			return Input::NO_CHARGES;
  if (is_equal(str, "openbabel" )) 	return Input::OPENBABEL;
  if (is_equal(str, "compass"))			return Input::FROM_COMPASS;

  return Input::NOT_SPECIFIED;
}

/// @version Specialization for a Input::NeighboursMethod enumeration type
template <> Input::NeighMethod from_string(const std::string& field,const std::string& str) {

  if (is_equal(str, "verlet_list" )) 			return Input::VERLET_LIST;
  if (is_equal(str, "block_verlet" )) 			return Input::BLOCK_VERLET;
  return Input::VERLET_LIST;
}
