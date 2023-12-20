# =================================
# === micro-stamp Analytics library  ===
# =================================

set(exaStampAnalytics_LINK_LIBRARIES exaStampMechanical exaStampCompute exanbParticleNeighbors exanbMPI exaStampParticleSpecies exanbGridCellParticles)
xstamp_add_plugin(exaStampAnalytics ${CMAKE_CURRENT_SOURCE_DIR})

