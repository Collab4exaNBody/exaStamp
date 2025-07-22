
set(EAM_POTENTIAL_TEMPALE_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/eam_potential.cu
  ${CMAKE_CURRENT_SOURCE_DIR}/eam_potential_multimat.cu
  PARENT_SCOPE)

macro(AddEAMPotential EamPotDir PluginName)
  list(APPEND ${PluginName}_LINK_LIBRARIES exaStampEAM exaStampPotential exanbParticleNeighbors)
  list(APPEND ${PluginName}_EXTERNAL_SRCS ${EAM_POTENTIAL_TEMPALE_SRCS})
  onika_add_plugin(${PluginName} ${EamPotDir})
endmacro()

