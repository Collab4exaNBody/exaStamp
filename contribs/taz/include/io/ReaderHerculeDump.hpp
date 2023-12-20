/// @file 
/// @brief Definition of a tool to read particles from dump files with Hercule
// ECom = to be commented or deleted by Estelle

#ifndef __PARTICLE_READER_HERCULE_DUMP_HPP_INCLUDED
#define __PARTICLE_READER_HERCULE_DUMP_HPP_INCLUDED


#include "HIc.h"

#include "io/StampV3LegacyIOStructures.hpp"

#include "utils/array/array.hpp"


HIC_USE;


class CommManager;


/// @brief Tool to read particles from a dump file with Hercule
class ReaderHerculeDump {

private:
	CommManager* comm; ///< Communication Manager

public:
	/// @brief Constructor
	/// @param [in] addrCommManager Pointer to the communication manager
	/// @param [in] _api ECom
	/// @param [in] _rep ECom
	ReaderHerculeDump(CommManager* addrCommManager,const HIc_Api& _api,
			const std::string& _rep) :
			comm(addrCommManager), m_rep(_rep), m_api(_api) {
		  if (getenv("XSP_MODE_INIT_HERCULE")){
		  	  m_modeInitHercule=atoi(getenv("XSP_MODE_INIT_HERCULE"));
		  }
	}
	/// @brief Destructor (nothing to do)
	~ReaderHerculeDump() {
	}

	void readHeader(int _numSDom, Array<LegacyHeaderIOStruct>& header);

	void readParticles(int _numSDom,HerculeParticleIODumpStruct& particles);
  	void readParticles(int _numSDom,HerculeDPDEParticleIODumpStruct& particles);

	std::string m_rep; ///< ECom
	int m_modeInitHercule=0; ///< ECom

	HIc_Api m_api; ///< ECom
	HIc_Base m_base; ///< ECom
	Array<float_8> m_steps; ///< Ecom
	Array<int_8> m_nbDoms; ///< ECom

	map<int,HIc_Ctx> m_ctx; ///< ECom
	void openBase();
	int nbStep(void);
	vector<int> steps();
	int nbSSDom(int _step=-1);
	void closeBase();
	void readMetaInfo();

	void openStep(uint _step, int _numSDom=-1);
	void closeStep(int _numSDom);
	bool isReader(int _numSDom);
	int masterIo(int _numSDom);
};

#endif // __PARTICLE_OUTPUT_HPP_INCLUDED
