/// @file
/// @brief Implementation of the InputOutputManager class


#include <string>
#include <iostream>

#include "io/outputManager.hpp"
#ifdef __use_lib_hercule
#include "io/particleOutputHercule.hpp"
#include "io/WriterHerculeDump.hpp"
#include "io/ReaderHerculeDump.hpp"
#include "HIc.h"
#include "HIc_MPI.h"

#include "referenceMap.hpp"

HIC_USE;


#endif


#ifdef __use_lib_hercule


/// @brief [ECom] InputOutputManager::initHerculeDumpR
void InputOutputManager::initHerculeDumpR(){
	initHercule();
	if (!m_herculeDumpR){
		m_herculeDumpR = new ReaderHerculeDump(commManager, m_api,m_dumpDir);
	}
}


/// @brief [ECom] InputOutputManager::initHerculeDumpW
void InputOutputManager::initHerculeDumpW(){
	initHercule();
	if (!m_herculeDumpW){
		m_herculeDumpW = new WriterHerculeDump(commManager, m_api,m_dumpDir);
	}
}


/// @brief [ECom] InputOutputManager::initHerculeDepW
void InputOutputManager::initHerculeDepW(){
	initHercule();
	if (!m_herculeDepW){
	  m_herculeDepW = new ParticleWriterHercule(commManager, m_api, true, true,Global::reference.isDPD() || Global::reference.isSDPD(),
						    Global::reference.isReactive(),m_outputDir, "xstamp-particles");
	}
}


/// @brief [ECom] InputOutputManager::initHercule
void InputOutputManager::initHercule() {
	if (!m_api.isNull())
		return;
    if (getenv("CODE") == NULL){
    	putenv(strdup("CODE=exaStamp"));
    }
    if (getenv("GME_DEFAULT_MODE") == NULL){
    	putenv(strdup("GME_DEFAULT_MODE=gmev2"));
    }
	/* initialisation de l'Api */
	m_api = HIc_Api("myApi", "stdDif");

	HIc_Init_Standard_Services(m_api);
	HIc_Init_Standard_CEA_DAM_DIF(m_api, "");
	HIc_Init_Parallel_MPI_Services(m_api, commManager->getCommunicator());
}


#endif


/// @brief Open an MPI_IO file
/// @param [in] name Name of the file
/// @param [in] rw Opening mode ("r" for read or "w" for write)
/// @param [out] mpiIOFileId File
void InputOutputManager::openIOFile(const std::string& name, const std::string& rw, FileId& mpiIOFileId){
  if (rw == "r")
    MPI_File_open(commManager->getCommunicator(), (char*)name.c_str(), MPI_MODE_RDONLY,                 MPI_INFO_NULL, &mpiIOFileId);
  else
    MPI_File_open(commManager->getCommunicator(), (char*)name.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &mpiIOFileId);
}


/// @brief Close an MPI_IO file
/// @param [in,out] mpiIOFileId File
void InputOutputManager::closeIOFile(FileId& mpiIOFileId){
  MPI_File_close(&mpiIOFileId);
}


/// @brief Read and broadcast the input file
/// @param [in] inputfile Input file
/// @param [out] datastream String that will temporarily store the input data
/// @param [in] comment Comment delimiter
void InputOutputManager::load(const std::string& inputfile, std::string& datastream, const char comment) {

  Array<int> streamSize(1);
  Array<char> tmp;

  // Job for the master node only
  if (m_isMaster) {

    // =========================================================================
    std::cout<< "  Reading input file";
    // =========================================================================

    // Get the file into a string
    readFromFile(inputfile, datastream, comment);

    // =========================================================================
    std::cout<< "\r  Reading input file ............... ok" << std::endl;
    // =========================================================================

    // =========================================================================
    std::cout<< "  Broadcasting input stream";
    // =========================================================================

    // Copy the datastream into a temporary char array
    tmp = Array<char>(datastream.size());
    memcpy(tmp.data(), datastream.data(), datastream.size());

    streamSize[0] = datastream.size();

  }

  // Broadcast the datastream size
  commManager->broadcast(streamSize);

  // Job for the non-master nodes
  if (!m_isMaster)
  	// Allocate place to get the datastream
    tmp = Array<char>(streamSize[0]);

  // Broadcast the temporary char array
  commManager->broadcast(tmp);

  // Copy the temporary char array into the datastream
  datastream = std::string(tmp.data(), tmp.size());

  // ===========================================================================
  if (m_isMaster) std::cout<< "\r  Broadcasting input stream ........ ok " << std::endl;
  // ===========================================================================

}


#ifdef __use_lib_hercule


/// @brief [ECom] InputOutputManager::herculeDumpR_openBase
void InputOutputManager::herculeDumpR_openBase(){
	initHerculeDumpR();
	m_herculeDumpR->openBase();
}


/// @brief [ECom] InputOutputManager::herculeDumpR_closeBase
void InputOutputManager::herculeDumpR_closeBase(){
	m_herculeDumpR->closeBase();
}


/// @brief [ECom] InputOutputManager::herculeDumpR_nbStep
int InputOutputManager::herculeDumpR_nbStep(){
	return m_herculeDumpR->nbStep();
}


/// @brief [ECom] InputOutputManager::herculeDumpR_steps
vector<int> InputOutputManager::herculeDumpR_steps(){
	return m_herculeDumpR->steps();
}


/// @brief [ECom] InputOutputManager::herculeDumpR_nbSSDom
int InputOutputManager::herculeDumpR_nbSSDom(uint _step){
	return m_herculeDumpR->nbSSDom(_step);
}


/// @brief [ECom] InputOutputManager::herculeDumpR_openStep
void InputOutputManager::herculeDumpR_openStep(uint _step, int _numSDom){
	m_herculeDumpR->openStep(_step,_numSDom);
}


/// @brief [ECom] InputOutputManager::herculeDumpR_closeStep
void InputOutputManager::herculeDumpR_closeStep(int _numSDom){
	m_herculeDumpR->closeStep(_numSDom);
}


/// @brief [ECom] InputOutputManager::herculeDumpR_readHeader
void InputOutputManager::herculeDumpR_readHeader(LegacyHeaderIOStruct* header, const std::string& filename, int _step, int _ssDom, InputOutputManager* ioManager){

	Array<LegacyHeaderIOStruct> headers(1);
	//fprintf(stderr,"ReaderHerculeDump::herculeDumpR_readHeader\n");

	herculeDumpR_openBase();
	if (_step == -2){
		vector<int> steps = herculeDumpR_steps();
		if (steps.size()== 0){
			std::cerr << "Error no hercule dump base" << std::endl;
			abort();
		}
		_step=steps[steps.size()-1];
		if (isMaster()){
			std::cerr << "herculeDumpR_readHeader: Restart last step in hercule dump base: " << _step << std::endl;
		}
	}
	herculeDumpR_openStep(_step,_ssDom);
	m_herculeDumpR->readHeader(_ssDom,headers);
	*header=headers[0];

  herculeDumpR_closeStep(_ssDom);

}


/// @brief [ECom] InputOutputManager::herculeDumpR_modeInitHercule
int InputOutputManager::herculeDumpR_modeInitHercule(){
	return m_herculeDumpR->m_modeInitHercule;
}


/// @brief [ECom] InputOutputManager::herculeDumpW_openBase
void InputOutputManager::herculeDumpW_openBase(){
	m_herculeDumpW->openBase();
}


/// @brief [ECom] InputOutputManager::herculeDumpW_closeBase
void InputOutputManager::herculeDumpW_closeBase(){
	m_herculeDumpW->closeBase();
}


/// @brief [ECom] InputOutputManager::herculeDumpW_openStep
void InputOutputManager::herculeDumpW_openStep(uint _step, int _numSDom){
	m_herculeDumpW->openStep(_step,_numSDom);
}


/// @brief [ECom] InputOutputManager::herculeDumpW_closeStep
void InputOutputManager::herculeDumpW_closeStep(){
	m_herculeDumpW->closeStep();
}

#endif
