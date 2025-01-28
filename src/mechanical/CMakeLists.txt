# =================================
# === micro-stamp Analytics library  ===
# =================================

set(exaStampMechanical_LINK_LIBRARIES exaStampCompute exanbDefBox exanbIO exanbParticleNeighbors exaStampSnapLegacy exaStampParticleSpecies)
onika_add_plugin(exaStampMechanical ${CMAKE_CURRENT_SOURCE_DIR})

