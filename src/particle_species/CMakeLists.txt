
set(XSTAMP_MAX_RIGID_MOLECULE_ATOMS "16" CACHE STRING "Maximum atoms per rigid molecule")
set(exaStampParticleSpecies_COMPILE_DEFINITIONS -DXSTAMP_MAX_RIGID_MOLECULE_ATOMS=${XSTAMP_MAX_RIGID_MOLECULE_ATOMS})
set(exaStampParticleSpecies_LINK_LIBRARIES exanbCore exanbGridCellParticles exanbParticleNeighbors)

xstamp_add_plugin(exaStampParticleSpecies ${CMAKE_CURRENT_SOURCE_DIR})

