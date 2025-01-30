# =================================
# === micro-stamp core library  ===
# =================================

set(PAIR_POTENTIAL_TEMPALE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pair_potential_template)
file(GLOB PAIR_POTENTIAL_TEMPALE_SRCS ${PAIR_POTENTIAL_TEMPALE_DIR}/*.cpp ${PAIR_POTENTIAL_TEMPALE_DIR}/*.cu)
macro(AddPairPotential PairPotDir PluginName)
  list(APPEND ${PluginName}_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
  #  list(APPEND ${PluginName}_INCLUDE_DIRS ${PAIR_POTENTIAL_TEMPALE_DIR})
  list(APPEND ${PluginName}_EXTERNAL_SRCS ${PAIR_POTENTIAL_TEMPALE_SRCS})
  onika_add_plugin(${PluginName} ${PairPotDir})
endmacro()

add_subdirectory(reaction_field)
add_subdirectory(ljexp6rf)

# Pair potentials
add_subdirectory(pair_potentials)
add_subdirectory(multimat_pair_potential)

# EAM potentials
#add_subdirectory(eam)
#add_subdirectory(eam_potentials)

#add_subdirectory(LCHBOP)
#add_subdirectory(meam)
#add_subdirectory(coul_dsf)
#add_subdirectory(coul_wolf)
#add_subdirectory(ewald)
#add_subdirectory(snap)
#add_subdirectory(snaplegacy)
#add_subdirectory(snaplmp)
#add_subdirectory(range_potential)
#add_subdirectory(n2p2)

#GenerateAllPairPotentials()
#GenerateAllEAMPotentials()

