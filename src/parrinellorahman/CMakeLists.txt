# ====================================================================
# === micro-stamp pair potential compute operator library library  ===
# ====================================================================

set(exaStampParrinelloRahman_LINK_LIBRARIES exaStampCompute exaStampParticleSpecies)
xstamp_add_plugin(exaStampParrinelloRahman ${CMAKE_CURRENT_SOURCE_DIR})

