# =================================
# === micro-stamp core library  ===
# =================================
find_package(ZLIB REQUIRED)

set(exaStampIO_LINK_LIBRARIES exaStampMechanical exaStampCompute exanbMPI exanbIO exanbDefBox exaStampParticleSpecies exanbGridCellParticles exanbParticleNeighbors ${ZLIB_LIBRARIES} stdc++fs)
set(exaStampIO_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR} ${BASEN_INCLUDE_DIRS}) 
onika_add_plugin(exaStampIO ${CMAKE_CURRENT_SOURCE_DIR})

