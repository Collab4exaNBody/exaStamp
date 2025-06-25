# ====================================================================
# === micro-stamp pair potential compute operator library library  ===
# ====================================================================

set(exaStampRigidMolecule_LINK_LIBRARIES exaStampCompute exaStampPotential exanbParticleNeighbors)
onika_add_plugin(exaStampRigidMolecule ${CMAKE_CURRENT_SOURCE_DIR})

