# ===================================
# === micro-stamp meam potential  ===
# ===================================
set(exaStampEwald_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
onika_add_plugin(exaStampEwald ${CMAKE_CURRENT_SOURCE_DIR})

if(XSTAMP_BUILD_exaStampEwald)
	set(ewald_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/include PARENT_SCOPE)
  AddPairPotential(${CMAKE_CURRENT_SOURCE_DIR}/ewald_pair_potential ewald)
endif()

