
set(XSTAMP_MAX_RIGID_MOLECULE_ATOMS "4" CACHE STRING "Maximum atoms per rigid molecule")
set(exaStampParticleSpecies_COMPILE_DEFINITIONS XSTAMP_MAX_RIGID_MOLECULE_ATOMS=${XSTAMP_MAX_RIGID_MOLECULE_ATOMS})
set(exaStampParticleSpecies_LINK_LIBRARIES exanbCore exanbGridCellParticles exanbParticleNeighbors)

onika_add_plugin(exaStampParticleSpecies ${CMAKE_CURRENT_SOURCE_DIR})

