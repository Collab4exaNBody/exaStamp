# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(exaStampSnapLMP_LINK_LIBRARIES exaStampSnap exaStampPotential exanbParticleNeighbors)
set(XSTAMP_BUILD_exaStampSnapLMP_DEFAULT OFF) # disbaled, have to be ativated by hand
xstamp_add_plugin(exaStampSnapLMP ${CMAKE_CURRENT_SOURCE_DIR})

