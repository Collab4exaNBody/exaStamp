# ===================================
# === micro-stamp meam potential  ===
# ===================================

set(LMP_SRC_DIR "${XSTAMP_THIRD_PARTY_TOOLS_ROOT}/lammps" CACHE STRING "LAMMPS source directory")

set(exaStampSnapLMP_EXTERNAL_SRCS ${LMP_SRC_DIR}/src/memory.cpp ${LMP_SRC_DIR}/src/math_special.cpp ${LMP_SRC_DIR}/src/ML-SNAP/sna.cpp)
set(exaStampSnapLMP_EXTERNAL_HDRS ${LMP_SRC_DIR}/src/memory.h ${LMP_SRC_DIR}/src/math_const.h ${LMP_SRC_DIR}/src/math_special.h ${LMP_SRC_DIR}/src/ML-SNAP/sna.h)
set(exaStampSnapLMP_LINK_LIBRARIES exaStampSnap exaStampPotential exanbParticleNeighbors)

set(XSTAMP_BUILD_exaStampSnapLMP_DEFAULT OFF) # disbaled, have to be ativated by hand
onika_add_plugin(exaStampSnapLMP ${CMAKE_CURRENT_SOURCE_DIR})

