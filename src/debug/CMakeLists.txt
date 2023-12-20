# =================================
# === micro-stamp core library  ===
# =================================

set(exaStampDebug_LINK_LIBRARIES exanbCore exanbDefBox exanbParticleNeighbors exaStampParticleSpecies exanbCompute)
xstamp_add_plugin(exaStampDebug ${CMAKE_CURRENT_SOURCE_DIR})

