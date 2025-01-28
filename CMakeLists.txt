cmake_minimum_required(VERSION 3.26)

# ========================================
# === Compiler toolchain configuration ===
# ========================================
# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS NO)

# ===============
# === Project ===
# ===============
project(exaStamp VERSION 3.6 LANGUAGES C CXX)
set(EXASTAMP_VERSION "exaSPH-${CMAKE_PROJECT_VERSION}") 
message(STATUS "Project version : ${EXASPH_VERSION}")

# ===========================
# === CMake customization ===
# ===========================
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release RelWithDebInfo Debug)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ======================================
# === exaStamp package specializations ===
# ======================================
set(XNB_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/include)
set(XNB_COMPILE_DEFINITIONS XNB_DOMAIN_SPECIFIC_FIELDS_INCLUDE=\<exaStamp/fields.h\>)
set(XNB_APP_HAS_GRID_SPECIALIZATIONS ON)

# ======================================
# === exaNBody application framework ===
# ======================================
find_package(exaNBody)

# ====================================
# === CCC supercomputing facility  ===
# ====================================
option(EXASTAMP_CCC_USER "Use CCC job resource interface" OFF)
if(EXASTAMP_CCC_USER)
  set(CCC_USER_DEFINITIONS "-D__use_lib_ccc_user=1" CACHE STRING "ccc_user macro definitions") 
  set(CCC_USER_INCLUDE_DIRS "" CACHE STRING "ccc_user include dir") 
  set(CCC_USER_LIBRARIES "ccc_user" CACHE STRING "ccc_user lib") 
  set(XNB_APP_DEFINITIONS ${XNB_APP_DEFINITIONS} ${CCC_USER_DEFINITIONS})
  set(XNB_APP_INCLUDE_DIRS ${XNB_APP_INCLUDE_DIRS} ${CCC_USER_INCLUDE_DIRS})
  set(XNB_APP_LIBRARIES ${XNB_APP_LIBRARIES} ${CCC_USER_LIBRARIES})
endif()

# ==========================================
# === build application plugins and libs ===
# ==========================================
add_subdirectory(src)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/exaStamp DESTINATION include)


# =================================================
# === generate and install config and test data ===
# =================================================
install(DIRECTORY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_INSTALL_PREFIX})
if(ONIKA_REGRESSION_TEST_ENABLE)
  AddRegressionTestDir("${CMAKE_SOURCE_DIR}/data/regression")
endif()

# =======================================
# === contributed tools and libraries ===
# =======================================
option(EXASTAMP_BUILD_CONTRIBS "Build contrib tools" OFF)
if(EXASTAMP_BUILD_CONTRIBS)
  add_subdirectory(contribs)
endif()

# ===================================
# === Documentation               ===
# ===================================
option(EXASTAMP_BUILD_DOC "Build documentation" OFF)
if(EXASTAMP_BUILD_DOC)
  add_subdirectory(doc)
endif()

# ===========================================================
# === generate execution wrappers and environment scripts ===
# ===========================================================
GenerateEndUserExecEnv()


