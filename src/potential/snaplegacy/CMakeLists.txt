# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(exaStampSnapLegacy_COMPILE_DEFINITIONS -DLAMMPS=1 ${KMP_ALIGNED_ALLOCATOR_DEFINITIONS})
set(exaStampSnapLegacy_INCLUDE_DIRS ${SNAP_INCLUDE_DIRS})
set(exaStampSnapLegacy_LINK_LIBRARIES exaStampSnap exaStampPotential exanbParticleNeighbors ${SNAP_LIBRARIES})
xstamp_add_plugin(exaStampSnapLegacy ${CMAKE_CURRENT_SOURCE_DIR})

