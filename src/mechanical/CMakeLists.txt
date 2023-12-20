# =================================
# === micro-stamp Analytics library  ===
# =================================

set(exaStampMechanical_LINK_LIBRARIES exaStampCompute exanbDefBox exanbIO exanbParticleNeighbors exaStampSnapLegacy exaStampParticleSpecies)
xstamp_add_plugin(exaStampMechanical ${CMAKE_CURRENT_SOURCE_DIR})

