# ===================================
# ===  multi-material potential   ===
# ===================================

set(exaStampMultimat_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
onika_add_plugin(exaStampMultimat ${CMAKE_CURRENT_SOURCE_DIR})

