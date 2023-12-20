# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(exaStampSnap_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
xstamp_add_plugin(exaStampSnap ${CMAKE_CURRENT_SOURCE_DIR})

