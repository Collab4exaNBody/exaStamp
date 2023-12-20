# =================================
# === micro-stamp core library  ===
# =================================


set(exaStampMolecule_LINK_LIBRARIES exanbCompute exaStampCompute exanbParticleNeighbors exaStampParticleSpecies lennard_jones exp6 exaStampReactionFieldParticleCharge exaStampPotential)
xstamp_add_plugin(exaStampMolecule ${CMAKE_CURRENT_SOURCE_DIR})

