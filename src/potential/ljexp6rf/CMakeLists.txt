# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(exaStampLJExp6RFParticleCharge_LINK_LIBRARIES exaStampCompute exp6 lennard_jones exaStampReactionFieldParticleCharge exaStampPotential exanbParticleNeighbors)
xstamp_add_plugin(exaStampLJExp6RFParticleCharge ${CMAKE_CURRENT_SOURCE_DIR})

