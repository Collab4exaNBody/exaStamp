# =================================
# === micro-stamp mpi library  ===
# =================================

set(exaStampMPI_COMPILE_DEFINITIONS ${EXASTAMP_ZOLTAN_DEFINITIONS})
set(exaStampMPI_INCLUDE_DIRS ${EXASTAMP_ZOLTAN_INCLUDE_DIRS})
set(exaStampMPI_LINK_LIBRARIES exanbMPI exaStampCompute exaStampParticleSpecies exanbGridCellParticles ${EXASTAMP_ZOLTAN_LIBRARIES})
xstamp_add_plugin(exaStampMPI ${CMAKE_CURRENT_SOURCE_DIR})

