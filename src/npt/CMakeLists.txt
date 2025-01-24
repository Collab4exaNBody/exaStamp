# ====================================================================
# === micro-stamp pair potential compute operator library library  ===
# ====================================================================

set(exaStampNPT_LINK_LIBRARIES exaStampCompute exaStampParticleSpecies)
xstamp_add_plugin(exaStampNPT ${CMAKE_CURRENT_SOURCE_DIR})

