# ====================================================================
# === micro-stamp pair potential compute operator library library  ===
# ====================================================================

set(exaStampTTM_LINK_LIBRARIES exaStampParticleSpecies exanbGridCellParticles exanbCompute)
onika_add_plugin(exaStampTTM ${CMAKE_CURRENT_SOURCE_DIR})

