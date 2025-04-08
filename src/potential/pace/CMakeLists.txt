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
  set(WIGNER_PATH ${PACE_ROOT_DIR}/wigner-cpp)
  set(WIGNER_INCLUDE_PATH ${WIGNER_PATH}/include)
  
  add_library(cnpy SHARED ${PACE_ROOT_DIR}/cnpy/cnpy.cpp)  
  target_include_directories(cnpy PUBLIC ${PACE_ROOT_DIR}/cnpy)

  file(GLOB exaStampPace_EXTERNAL_SRCS
    ${PACE_ROOT_DIR}/ML-PACE/ace-evaluator/*.cpp
    ${PACE_ROOT_DIR}/ML-PACE/ace/*.cpp
  )
  add_library(pace STATIC ${exaStampPace_EXTERNAL_SRCS})
  target_include_directories(pace PUBLIC ${WIGNER_INCLUDE_PATH})
  target_include_directories(pace PUBLIC ${PACE_ROOT_DIR}/ML-PACE)
  target_compile_definitions(pace PUBLIC EXTRA_C_PROJECTIONS)
  target_compile_definitions(pace PUBLIC COMPUTE_B_GRAD)
  target_link_libraries(pace PUBLIC cnpy)

  set_target_properties(pace PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
  )
  
  set(exaStampPace_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors cnpy pace)
  onika_add_plugin(exaStampPace ${CMAKE_CURRENT_SOURCE_DIR})
endif()
