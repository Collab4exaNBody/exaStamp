/// @file 
/// @brief Definition of a tool to write dump files with Hercule
// ECom = to be commented or deleted by Estelle

#ifndef __PARTICLE_WRITER_HERCULE_DUMP_HPP_INCLUDED
#define __PARTICLE_WRITER_HERCULE_DUMP_HPP_INCLUDED


#include "HIc.h"

#include "io/StampV3LegacyIOStructures.hpp"

#include "utils/array/array.hpp"


HIC_USE;


class CommManager;


/// @brief Tool to write dump file with Hercule
class WriterHerculeDump {

private:
  CommManager* comm; ///< Communication manager

public:
  /// @brief Constructor
  /// @param [in] addrCommManager Pointer to the communication manager
  /// @param [in] _api [ECom]
  /// @param [in] _rep [ECom]
  WriterHerculeDump(CommManager* addrCommManager,const HIc_Api& _api,
		    const std::string& _rep) :
    comm(addrCommManager), m_rep(_rep), m_api(_api) {
  }
  /// @brief Destructor (nothing to do)
  ~WriterHerculeDump() {
  }

  template <class DumpStruct>
  void write(Array<LegacyHeaderIOStruct>& header,
	     DumpStruct& particles, int step,
	     double cpuTime, double physTime);

  void writeHeader(Array<LegacyHeaderIOStruct>& header);

  void writeParticles(HerculeParticleIODumpStruct& particles);
  void writeParticles(HerculeDPDEParticleIODumpStruct& particles);

  std::string m_rep; ///< ECom

  HIc_Api m_api; ///< ECom
  HIc_Base m_base; ///< ECom
  HIc_Ctx m_ctx; ///< ECom
  void openBase();
  void closeBase();
  void openStep(uint _step, int _numSSDom=-1);
  void closeStep();
};

/// @brief [ECom] WriterHerculeDump::write
template <class DumpStruct>
void WriterHerculeDump::write(Array<LegacyHeaderIOStruct>& header,
			      DumpStruct& particles,int step, double cpuTime, double physTime) {

  openBase();

  openStep(step);

  writeHeader(header);

  writeParticles(particles);

  closeStep();

}

#endif // __PARTICLE_OUTPUT_HPP_INCLUDED
