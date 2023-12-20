# =================================
# === micro-stamp core library  ===
# =================================

set(exaStampSliplink_LINK_LIBRARIES exanbMPI exaStampParticleSpecies)
xstamp_add_plugin(exaStampSliplink ${CMAKE_CURRENT_SOURCE_DIR})

