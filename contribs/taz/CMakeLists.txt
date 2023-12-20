# ==========================================
# === taz : an AMR exaStamp spinoff      ===
# ==========================================

# =======================
# === build & install ===
# =======================
set(TAZ_GHOST_EXCHANGE_BUFFER_SIZE 50000 CACHE STRING "Number of particle per buffer chunk during ghost exchange")

# compile options
set(TAZ_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_BINARY_DIR}/gen/include
  ${CMAKE_CURRENT_SOURCE_DIR}/deps/gsl
  ${YAML_CPP_INCLUDE_DIR}
  ${SOATL_INCLUDE_DIRS}
  ${EXASTAMP_ZOLTAN_INCLUDE_DIRS}
  ${MPI_CXX_INCLUDE_PATH}
  ${EXASTAMP_HWLOC_INCLUDE_DIRS}
  )
# message(STATUS "TAZ_INCLUDE_DIRS = ${TAZ_INCLUDE_DIRS}")

set(TAZ_DEFINITIONS
  ${EXASTAMP_HOST_HW_DEFINITIONS}
  ${EXASTAMP_TEST_DEFINITIONS}
  -DVERSION="1.0"
  -D__xsp_compute_rep="${CMAKE_INSTALL_PREFIX}"
  -DTAZ_GHOST_EXCHANGE_BUFFER_SIZE=${TAZ_GHOST_EXCHANGE_BUFFER_SIZE}
  -D__use_lib_omp=1
  -D__evi_ser -D__evi_no_fma
  ${EXASTAMP_ZOLTAN_DEFINITIONS}
  ${EXASTAMP_HWLOC_DEFINITIONS}
  ${CCC_USER_DEFINITIONS}
  )
set(TAZ_CXXFLAGS
  ${OpenMP_CXX_FLAGS}
  ${METIS_CXX_FLAGS}
  )
set(TAZ_LIBRARIES
  ${EXASTAMP_HOST_HW_LIBRARIES}
  ${EXASTAMP_TEST_LIBRARIES}
  ${EXASTAMP_ZOLTAN_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  ${OpenMP_CXX_LIBRARIES}
  ${EXASTAMP_HWLOC_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${CCC_USER_LIBRARIES}
  )
# message(STATUS "TAZ_LIBRARIES=${TAZ_LIBRARIES}")
# Source files
file(GLOB_RECURSE TAZ_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Taz binary 
add_executable(taz ${TAZ_SRCS})
target_compile_definitions(taz PRIVATE ${TAZ_DEFINITIONS})
target_compile_options(taz PRIVATE ${TAZ_CXXFLAGS})
target_include_directories(taz PRIVATE ${TAZ_INCLUDE_DIRS})
target_link_libraries(taz ${TAZ_LIBRARIES})

install(TARGETS taz DESTINATION bin)

# generate launch script
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/taz-run.sh "#!/bin/sh\n${PROJECT_SETUP_ENV_COMMANDS}\n${CMAKE_INSTALL_PREFIX}/bin/taz \$*\n")
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/taz-run.sh DESTINATION bin)

file(GLOB input_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_SOURCE_DIR}/data/*.xsp.in)
foreach(inputfile ${input_files})
  get_filename_component(basename ${inputfile} NAME_WE)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/${inputfile} ${CMAKE_CURRENT_BINARY_DIR}/${basename}.xsp)
endforeach()

# tests
# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/regression_singlemat_lj_compv2.xsp.in ${CMAKE_CURRENT_BINARY_DIR}/regression_singlemat_lj_compv2.xsp)
set(TazSingleMatTestCommand ${CMAKE_CURRENT_BINARY_DIR}/taz ${CMAKE_CURRENT_BINARY_DIR}/regression_singlemat_lj_compv2.xsp)
AddTestWithDebugTarget(taz_singlemat_lj_compv2_mt TazSingleMatTestCommand 1 -1)

# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/regression_singlemat_meam_compv2.xsp.in ${CMAKE_CURRENT_BINARY_DIR}/regression_singlemat_meam_compv2.xsp)
set(TazSingleMatTestCommand ${CMAKE_CURRENT_BINARY_DIR}/taz ${CMAKE_CURRENT_BINARY_DIR}/regression_singlemat_meam_compv2.xsp)
AddTestWithDebugTarget(taz_singlemat_meam_compv2_mt TazSingleMatTestCommand 1 -1)

# configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/wall_meam_compv2.xsp.in ${CMAKE_CURRENT_BINARY_DIR}/wall_meam_compv2.xsp)
set(TazSingleMatTestCommand ${CMAKE_CURRENT_BINARY_DIR}/taz ${CMAKE_CURRENT_BINARY_DIR}/wall_meam_compv2.xsp)
AddTestWithDebugTarget(taz_wall_meam_compv2_mt TazSingleMatTestCommand 1 -1)

