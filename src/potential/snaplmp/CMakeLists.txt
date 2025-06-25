# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(LMP_SRC_DIR "${XSTAMP_THIRD_PARTY_TOOLS_ROOT}/lammps" CACHE STRING "LAMMPS source directory")

set(ONIKA_BUILD_exaStampSnapLMP_DEFAULT OFF) # disbaled, have to be ativated by hand
if(NOT IS_DIRECTORY ${LMP_SRC_DIR})
	message(STATUS "LAMMPS directory not found in ${LMP_SRC_DIR}, SnapLMP is forced OFF")
  set(ONIKA_BUILD_exaStampSnapLMP_FORCED ON)
else()
  set(ONIKA_BUILD_exaStampSnapLMP_FORCED OFF)
endif()

set(exaStampSnapLMP_EXTERNAL_SRCS ${LMP_SRC_DIR}/src/memory.cpp ${LMP_SRC_DIR}/src/math_special.cpp ${LMP_SRC_DIR}/src/ML-SNAP/sna.cpp)
set(exaStampSnapLMP_LINK_LIBRARIES exaStampSnap exaStampPotential exanbParticleNeighbors)
set(exaStampSnapLMP_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR} ${LMP_SRC_DIR}/src/ML-SNAP ${LMP_SRC_DIR}/src)
onika_add_plugin(exaStampSnapLMP ${CMAKE_CURRENT_SOURCE_DIR})

