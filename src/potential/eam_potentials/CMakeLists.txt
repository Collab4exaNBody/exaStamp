# ===========================================
# === micro-stamp EAM potential library  ===
# ===========================================

#define XSTAMP_EAM_CHECK_PAIR_DISTANCE 1.e-6
option(XSTAMP_EAM_CHECK_PAIR_DISTANCE "ExaStamp EAM potentials check whenever particle pair distance is below defined thresold" OFF)
if(XSTAMP_EAM_CHECK_PAIR_DISTANCE)
  set(XSTAMP_EAM_CHECK_PAIR_DISTANCE_VALUE "1.e-6" CACHE STRING "ExaStamp EAM potential minimal pair distance")
  set(USTAMP_COMPILE_DEFINITIONS "${USTAMP_COMPILE_DEFINITIONS} -DXSTAMP_EAM_CHECK_PAIR_DISTANCE=${XSTAMP_EAM_CHECK_PAIR_DISTANCE_VALUE}")
endif()

# Pair potentials generation from template sources
file(GLOB USTAMP_POT_SUBDIRS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*)
foreach(dir ${USTAMP_POT_SUBDIRS})
  if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
    AddEAMPotential(${CMAKE_CURRENT_SOURCE_DIR}/${dir} ${dir})
  endif()
endforeach()


