# =================================
# === micro-stamp core library  ===
# =================================
find_package(ZLIB REQUIRED)

set(exaStampIO_LINK_LIBRARIES exaStampCompute exanbMPI exanbIO exanbDefBox exaStampParticleSpecies exanbGridCellParticles exanbParticleNeighbors exaStampMolecule ${ZLIB_LIBRARIES} stdc++fs)
set(exaStampIO_INCLUDE_DIRS ${ZLIB_INCLUDE_DIR} ${BASEN_INCLUDE_DIRS}) 
xstamp_add_plugin(exaStampIO ${CMAKE_CURRENT_SOURCE_DIR})

