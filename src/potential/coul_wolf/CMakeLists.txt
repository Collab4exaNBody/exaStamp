# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(exaStampCoulWolfParticleCharge_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
xstamp_add_plugin(exaStampCoulWolfParticleCharge ${CMAKE_CURRENT_SOURCE_DIR})

