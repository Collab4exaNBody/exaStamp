
###
# Pair potential template replicated for all simple pair potentials
###

option(EXASTAMP_ENABLE_PAIR_WEIGHTING "Enables support for weighted pair potential" ON)

if(EXASTAMP_ENABLE_PAIR_WEIGHTING)
  set(PAIR_POTENTIAL_TEMPALE_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_class.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_common_vOFF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_common_vON.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_vOFF_wOFF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_vOFF_wON.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_vON_wOFF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_vON_wON.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_impl_multiparam.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_impl_rigidmol.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_impl_simple.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_plot.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_singlemat_symetric.cpp 
    PARENT_SCOPE )
else()
  set(PAIR_POTENTIAL_TEMPALE_SRCS
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_class.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_common_vOFF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_common_vON.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_vOFF_wOFF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_force_op_vON_wOFF.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_impl_multiparam.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_impl_rigidmol.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_impl_simple.cu
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_plot.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_singlemat_symetric.cpp 
    PARENT_SCOPE )
endif()

macro(AddPairPotential PairPotDir PluginName)
  list(APPEND ${PluginName}_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
  list(APPEND ${PluginName}_COMPILE_DEFINITIONS EXASTAMP_ENABLE_PAIR_WEIGHTING=$<BOOL:${EXASTAMP_ENABLE_PAIR_WEIGHTING}>)
  list(APPEND ${PluginName}_EXTERNAL_SRCS ${PAIR_POTENTIAL_TEMPALE_SRCS})
  onika_add_plugin(${PluginName} ${PairPotDir})
endmacro()

