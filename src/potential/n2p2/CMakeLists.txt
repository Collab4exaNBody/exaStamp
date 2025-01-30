# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(N2P2_ROOT_DIR /usr/local/n2p2 CACHE STRING "N2P2 installation directory")

# default is disabled, user must explicitly enable it
set(ONIKA_BUILD_exaStampN2P2_DEFAULT OFF)

if(NOT IS_DIRECTORY ${N2P2_ROOT_DIR})
	message(STATUS "N2P2 directory not found in ${N2P2_ROOT_DIR}, N2P2 is forced OFF")
  set(ONIKA_BUILD_exaStampN2P2_FORCED ON)
else()
  set(ONIKA_BUILD_exaStampN2P2_FORCED OFF)
endif()

# ${MPI_CXX_LIBRARIES}
set(microStampN2P2_INCLUDE_DIRS ${N2P2_ROOT_DIR}/include)
set(microStampN2P2_LINK_LIBRARIES exaStampPotential ${N2P2_ROOT_DIR}/lib/libnnp.so ${N2P2_ROOT_DIR}/lib/libnnpif.so)
onika_add_plugin(exaStampN2P2 ${CMAKE_CURRENT_SOURCE_DIR})

