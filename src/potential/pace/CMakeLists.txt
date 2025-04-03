# ====================================
# ===== exaStamp PACE potential  =====
# ====================================

set(PACE_ROOT_DIR /usr/local/pace CACHE STRING "PACE installation directory")

if(IS_DIRECTORY ${PACE_ROOT_DIR})
  message(STATUS "PACE directory found in ${PACE_ROOT_DIR}, PACE is forced ON")
  set(ONIKA_BUILD_exaStampPace_FORCED ON)
else()
  message(STATUS "PACE directory not found in ${PACE_ROOT_DIR}, PACE is forced OFF")
  set(ONIKA_BUILD_exaStampPace_FORCED OFF)
endif()

set(ONIKA_BUILD_exaStampPace ${ONIKA_BUILD_exaStampPace_FORCED} CACHE BOOL "Build exaStampPace")

if(ONIKA_BUILD_exaStampPace)
  include_directories(
    ${PACE_ROOT_DIR}/wigner-cpp/include
    ${PACE_ROOT_DIR}/cnpy
    ${PACE_ROOT_DIR}/ML-PACE
  )
  message(STATUS "Building cnpy")
#  add_library(cnpy STATIC ${PACE_ROOT_DIR}/cnpy/cnpy.cpp)
#  set_target_properties(cnpy PROPERTIES POSITION_INDEPENDENT_CODE ON)
  add_library(cnpy SHARED ${PACE_ROOT_DIR}/cnpy/cnpy.cpp)  
  target_include_directories(cnpy PUBLIC ${PACE_ROOT_DIR}/cnpy)

  file(GLOB exaStampPace_EXTERNAL_SRCS
    ${PACE_ROOT_DIR}/ML-PACE/ace/*.cpp
    ${PACE_ROOT_DIR}/ML-PACE/ace-evaluator/*.cpp
  )
  add_compile_definitions(EXTRA_C_PROJECTIONS)
  add_compile_definitions(COMPUTE_B_GRAD)
  
  set(exaStampPace_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors cnpy)
  onika_add_plugin(exaStampPace ${CMAKE_CURRENT_SOURCE_DIR})
endif()
