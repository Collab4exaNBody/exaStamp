# =================================
# === micro-stamp Analytics library  ===
# =================================

set(exaStampAnalytics_LINK_LIBRARIES exaStampCompute exaStampMechanical exanbParticleNeighbors exanbMPI exaStampParticleSpecies exanbGridCellParticles)
onika_add_plugin(exaStampAnalytics ${CMAKE_CURRENT_SOURCE_DIR})

