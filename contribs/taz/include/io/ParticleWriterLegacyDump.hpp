/// @file 
/// @brief Definition of a tool to write legacy dump files

#ifndef __PARTICLE_OUTPUT_LEGACYDUMP_HPP_INCLUDED
#define __PARTICLE_OUTPUT_LEGACYDUMP_HPP_INCLUDED


#include "io/outputManager.hpp"


class CommManager;


/// @brief Tool to write the legacy dump file
class ParticleWriterLegacyDump {

private:
	std::string legacyDumpName; ///< Name of the file
	InputOutputManager::FileId mpiIOFileId; ///< File
	InputOutputManager::Offset localOffset; ///< Position where this node must write in the dump file
	CommManager* comm; ///< Communication manager

public:
	/// @brief Constructor
	/// @param [in] addrCommManager Pointer to the communication manager
	/// @param [in] _inputOutput Pointer to the input/output manager
	ParticleWriterLegacyDump(CommManager* addrCommManager,InputOutputManager* _inputOutput) :
			comm(addrCommManager),m_inputOutputManager(_inputOutput) {
	}
	/// @brief Destructor (nothing to do)
	~ParticleWriterLegacyDump() {
	}

	void setFileName(uint step);
	void open();
	void writeHeader(Array<LegacyHeaderIOStruct>& header);
	void writeParticles(Array<LegacyParticleIOStruct>& particlesArray);
  	void writeParticles(Array<LegacyDPDEParticleIOStruct>& particlesArray);
	void close();

	void write(uint step, Array<LegacyHeaderIOStruct>& header,
			Array<LegacyParticleIOStruct>& particlesArray);

  	void write(uint step, Array<LegacyHeaderIOStruct>& header,
			Array<LegacyDPDEParticleIOStruct>& particlesArray);

	InputOutputManager* m_inputOutputManager; ///< Input/output manager
};

#endif // __PARTICLE_OUTPUT_HPP_INCLUDED
