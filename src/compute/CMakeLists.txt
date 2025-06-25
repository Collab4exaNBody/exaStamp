# =================================
# === micro-stamp core library  ===
# =================================
set(exaStampCompute_LINK_LIBRARIES exanbCompute exaStampParticleSpecies exanbGridCellParticles)
onika_add_plugin(exaStampCompute ${CMAKE_CURRENT_SOURCE_DIR})

