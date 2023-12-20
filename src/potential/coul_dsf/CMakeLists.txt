# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(exaStampCoulDsfParticleCharge_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
xstamp_add_plugin(exaStampCoulDsfParticleCharge ${CMAKE_CURRENT_SOURCE_DIR})

