/// @file 
/// @brief Definition of the InputOutputManager class
///
/// Also include the configuration structure used to construct an InputOutputManager
// ECom = to be commented or deleted by Estelle

#ifndef __OUTPUT_HPP_INCLUDED
#define __OUTPUT_HPP_INCLUDED


#include <string>

#include "parallel/commManager.hpp"
#include "parallel/mympi.hpp"

#include "io/input.hpp"
#ifdef __use_lib_hercule
#include "io/ReaderHerculeDump.hpp"
#include "HIc.h"


HIC_USE;
#endif


class InputOutputManager;
class ParticleWriterHercule;
class ReaderHerculeDump;
class WriterHerculeDump;


/// @brief Temporary structure gathering all data to initialize the input/output manager
template <> struct Configuration<InputOutputManager> {

	/// @brief The default constructor
	///
	///
  Configuration() 
    : m_outputType(HERCULE_PARTICULE_DEP), m_dumpType(HERCULE_DUMP), m_initType(), 
      m_inputDir("."), m_outputDir("."), m_dumpDir("."), m_outputCells(false) {}

  /// @brief Destructor (nothing to do)
  ~Configuration() {}

  /// @brief A constructor from the program input data
  ///
	///
  Configuration(const Input& input)
    : m_outputType(input.outputType),
      m_dumpType(input.dumpType),
      m_initType(input.m_initType),
      m_inputDir(input.initDir),
      m_outputDir(input.outputDir),
      m_dumpDir(input.dumpDir),
      m_outputCells(input.outputCells),
      m_initStep(input.initStep),
      m_checkSimulationValues(input.checkValues),
      m_checkSimulationEpsilon(input.checkEpsilon),
      m_checkSimulationFile(input.checkFile),
      m_checkSimulationIteration(input.checkIteration)
 {
    if (m_outputDir == "") m_outputDir=input.dir;
    if (m_dumpDir == "") m_dumpDir=input.dir;
    if (m_inputDir == "") m_inputDir=input.dir;
  }

  OutputType m_outputType; ///< Type of output file
  DumpType m_dumpType; ///< Type of dump file
  InitType m_initType; ///< Type of particles initialization
  std::string m_inputDir; ///< Initialization file directory
  std::string m_outputDir; ///< Output file directory
  std::string m_dumpDir; ///< Dump file directory
  bool m_outputCells; ///< Flag for the cell output
  int m_initStep; ///< Initial step
  bool m_checkSimulationValues;
  double m_checkSimulationEpsilon;
  std::string m_checkSimulationFile;
  int m_checkSimulationIteration;
};


/// @brief Inuput/output manager : the tool that handle all inputs and outputs communication, except debug ones
class InputOutputManager {

public:

  /// @brief Constructor from a configuration and a communication manager
  /// @param [in,out] _commManager Communication manager
  /// @param [in] configuration Input configuration
  InputOutputManager(CommManager *_commManager,Configuration<InputOutputManager>& configuration)
    : m_isMaster(_commManager->getRank() == 0),
#ifdef __use_lib_hercule
      m_herculeDepW(NULL), m_herculeDumpW(NULL), m_herculeDumpR(NULL),
#endif
      commManager(_commManager),
      m_dumpType(configuration.m_dumpType),
      m_dumpDir(configuration.m_dumpDir),
      m_outputType(configuration.m_outputType),
      m_outputDir(configuration.m_outputDir),
      m_initType(configuration.m_initType),
      m_initStep(configuration.m_initStep),
      m_initDir(configuration.m_inputDir),
      m_checkSimulationValues(configuration.m_checkSimulationValues),
      m_checkSimulationEpsilon(configuration.m_checkSimulationEpsilon),
      m_checkSimulationFile(configuration.m_checkSimulationFile),
      m_checkSimulationIteration(configuration.m_checkSimulationIteration),
      m_outputCells(configuration.m_outputCells){}

  /// @brief Reinitialize the input/output manager from a configuration
  /// @param [in] configuration Input configuration
  void changeConf(Configuration<InputOutputManager>& configuration){
    m_dumpType=configuration.m_dumpType;
    m_dumpDir=configuration.m_dumpDir;
    m_outputType=configuration.m_outputType;
    m_outputDir=configuration.m_outputDir;
    m_initDir=configuration.m_inputDir;
    m_initType=configuration.m_initType;
    m_outputCells=configuration.m_outputCells;
    m_initStep=configuration.m_initStep;
    m_checkSimulationValues = configuration.m_checkSimulationValues;
    m_checkSimulationEpsilon = configuration.m_checkSimulationEpsilon;
    m_checkSimulationFile = configuration.m_checkSimulationFile;
    m_checkSimulationIteration = configuration.m_checkSimulationIteration;
  }

  /// @brief Destructor (nothing to do)
  ~InputOutputManager() {}

  /// @brief Accessor to the particles initialization type
  InitType initType()   { return m_initType; }
  /// @brief Accessor to the initialization file directory
  std::string initDir() {return m_initDir;}
  /// @brief Accessor to the initial step
  int initStep() {return m_initStep;}

  /// @brief Accessor to the output file type
  OutputType outputType() { return m_outputType; }
  /// @brief Accessor to the dump file type
  DumpType dumpType()   { return m_dumpType; }
  /// @brief Accessor to the output file directory
  std::string outputDir() {return m_outputDir;}
  /// @brief Accessor to the dump file directory
  std::string dumpDir() {return m_dumpDir;}
  /// @brief  Accessor to the cell output flag
  bool outputCells() {return m_outputCells;}

  /// @brief true if simulation results are to be compared with a reference file
  inline bool checkValues() const { return m_checkSimulationValues; }
  
  /// @brief method used to check results against a reference file
  inline double checkEpsilon() const { return m_checkSimulationEpsilon; }

  inline std::string checkFile() const { return m_checkSimulationFile; }
  inline int checkIteration() const { return m_checkSimulationIteration; }

  /// @brief Shortcut for an MPI file
  typedef MPI_File   FileId;
  /// @brief Shortcut for an MPI offset
  typedef MPI_Offset Offset;

  void openIOFile (const std::string& name, const std::string& rw, FileId& mpiIOFileId);
  void closeIOFile(FileId& mpiIOFileId);
  template <class T> static void read (Array<T>& out, FileId& mpiIOFileId, Offset& offSet);
  template <class T> static void write(Array<T>& in , FileId& mpiIOFileId, Offset& offSet);

  bool m_isMaster; ///< Indicates if the current node is the master node
  /// @brief Accessor to the master node flag
  bool isMaster(){return m_isMaster;}
  void load(const std::string& inputfile, std::string& datastream, const char comment);

#ifdef __use_lib_hercule
  ParticleWriterHercule *m_herculeDepW; ///< [ECom] Tool to write particles with Hercule
  WriterHerculeDump *m_herculeDumpW; ///< [ECom] Tool to write dump with Hercule
  ReaderHerculeDump *m_herculeDumpR; ///< [ECom] Tool to read dump with Hercule
  void herculeDumpR_openBase();
  void herculeDumpR_closeBase();
  int herculeDumpR_nbStep();
  std::vector<int> herculeDumpR_steps();
  int herculeDumpR_nbSSDom(uint _step=-1);
  void herculeDumpR_openStep(uint _step, int _numSDom=-1);
  void herculeDumpR_closeStep(int _numSDom);
  void herculeDumpR_readHeader(LegacyHeaderIOStruct* header, const std::string& filename, int _step, int _ssDom, InputOutputManager* ioManager);
  template <class DumpStruct> void herculeDumpR_readParticles(int _numSDom, DumpStruct& _particles);
  void herculeDumpW_openBase();
  void herculeDumpW_closeBase();
  void herculeDumpW_openStep(uint _step, int _numSDom=-1);
  void herculeDumpW_closeStep();
  int herculeDumpR_modeInitHercule();
//  void hercule_close(); Neither used nor implemented
  void initHerculeDumpR();
  void initHerculeDumpW();
  void initHerculeDepW();
#endif

private:

#ifdef __use_lib_hercule
  HIc_Api m_api; ///< ECom
  void initHercule();
#endif

  CommManager *commManager; ///< Communication manager
  DumpType m_dumpType; ///< Type of dump file
  std::string m_dumpDir; ///< Dump file directory
  OutputType m_outputType; ///< Type of output file
  std::string m_outputDir; ///< Output file directory
  InitType m_initType; ///< Type of particles initialization
  int m_initStep; ///< Initial step
  std::string m_initDir; ///< Initialization file directory
  bool m_outputCells; ///< Flag for the cell output
  bool m_checkSimulationValues;
  double m_checkSimulationEpsilon;
  std::string m_checkSimulationFile;
  int m_checkSimulationIteration;
};


/// @brief Parallel read an array in a file
/// @tparam T Class in the array
/// @param [out] out Array read
/// @param [in] mpiIOFileId File to read in
/// @param [in] offset Position to read at
template <class T>
void InputOutputManager::read(Array<T>& out, FileId& mpiIOFileId, Offset& offset) {
  MPI_Status status;
  MPI_File_read_at(mpiIOFileId, offset, out.data(), out.size(), MPI__Type_get<T>(), &status);
}


/// @brief Parallel write an array in a file
/// @tparam T Class in the array
/// @param [in] in Array to write
/// @param [in,out] mpiIOFileId File to write in
/// @param [in] offset Position to write at
template <class T>
void InputOutputManager::write(Array<T>& in, FileId& mpiIOFileId, Offset& offset){
  MPI_Status status;
  int cond = 0;
  int after;
  while(cond==0)
  {
    MPI_File_write_at(mpiIOFileId, offset, in.data(), in.size(), MPI__Type_get<T>(), &status);
    after=0;
    MPI_Get_count(&status, MPI__Type_get<T>(), &after);
    if(in.size()==after)
    {
      cond=1;
    }
  }
}

#ifdef __use_lib_hercule
/// @brief [ECom] InputOutputManager::herculeDumpR_readParticles
/// @tparam DumpStruct Dump structure (Hercule dump or Hercule dump for DPDE)
/// @param [in] _numSDom 
/// @param [in] _particles Particles data
template <class DumpStruct>
void InputOutputManager::herculeDumpR_readParticles(int _numSDom, DumpStruct& _particles)
{
	m_herculeDumpR->readParticles(_numSDom,_particles);
}
#endif

#endif // __OUTPUT_HPP_INCLUDED
