/// @file 
/// @brief Definition of a tool to print the particle output with Hercule
// ECom = to be commented or deleted by Estelle

#ifdef __use_lib_hercule
#ifndef __PARTICLE_OUTPUT_HERCULE_HPP_INCLUDED
#define __PARTICLE_OUTPUT_HERCULE_HPP_INCLUDED


#include "io/particleOutput.hpp"
#include "HIc.h"


HIC_USE;


/// @brief Tool to write the particles with Hercule
class ParticleWriterHercule : public ParticleWriter {
public:
  //ParticleWriterHercule() {}
  virtual ~ParticleWriterHercule();

  CommManager* comm; ///< Communication manager

public:

  /// @brief Constructor
  /// @param [in] addrCommManager Pointer to the communication manager
  /// @param [in] _api Ecom
  /// @param [in] writeT Indicates if the types are to be written
  /// @param [in] writeV Indicate if the velocities are to be written
  /// @param [in] writeE Indicate if the internal energies are to be written
  /// @param [in] writeP Indicate if the internal energies are to be written
  /// @param [in] rep Ecom
  /// @param [in] name Root name of the file where the particles will be written
  ParticleWriterHercule(CommManager* addrCommManager, const HIc_Api& _api, bool writeT, bool writeV, bool writeE, bool writeP, const std::string& rep, const std::string& name)
    : ParticleWriter(NULL, writeT, writeV, writeE, writeP, name, ".vtk"), comm(addrCommManager), m_rep(rep), m_api(_api) {
  }

  void write(int step,ParticleOutput* buff);

  void writeMetaInfo(HIc_Ctx& _ctx, int step);
  void writePositions(HIc_Ctx& _ctx,HIc_Obj&_obj_particleSet);
  void writeIndexes(HIc_Ctx& _ctx,HIc_Obj&_obj_particleSet);
  void writeTypes(HIc_Ctx& _ctx,HIc_Obj&_obj_particleSet);
  void writeVelocities(HIc_Ctx& _ctx,HIc_Obj&_obj_particleSet);
  void writeEints(HIc_Ctx& _ctx,HIc_Obj&_obj_particleSet);
  void writeProgresses(HIc_Ctx& _ctx,HIc_Obj&_obj_particleSet);

  std::string m_rep; ///< ECom

  HIc_Api m_api; ///< ECom
  HIc_Base m_base; ///< ECom
  HIc_Ctx m_ctx; ///< ECom
  void openBase();
  void closeBase();
  void openStep(uint,ParticleOutput* buff);
  void closeStep();
 
};

#endif // __PARTICLE_OUTPUT_HPP_INCLUDED
#endif
