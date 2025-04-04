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
set(EXASTAMP_VERSION ${CMAKE_PROJECT_VERSION}) 
message(STATUS "Project version : ${EXASTAMP_VERSION}")

# ===========================
# === CMake customization ===
# ===========================
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release RelWithDebInfo Debug)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ======================================
# === exaStamp package specializations ===
# ======================================
set(XNB_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/include)
set(XNB_COMPILE_DEFINITIONS XNB_DOMAIN_SPECIFIC_FIELDS_INCLUDE=\<exaStamp/fields.h\> EXASTAMP_VERSION=${EXASTAMP_VERSION})
set(XNB_APP_HAS_GRID_SPECIALIZATIONS ON)
option(EXASTAMP_ENABLE_MOLECULE "Enables support for rigid and flexible molecules" ON)
if(EXASTAMP_ENABLE_MOLECULE)
  set(XNB_COMPILE_DEFINITIONS ${XNB_COMPILE_DEFINITIONS} EXASTAMP_ENABLE_MOLECULE=1)
endif()
option(EXASTAMP_ENABLE_MECHANICAL "Enables support for mechanical analysis" ON)
if(EXASTAMP_ENABLE_MECHANICAL)
  set(XNB_COMPILE_DEFINITIONS ${XNB_COMPILE_DEFINITIONS} EXASTAMP_ENABLE_MECHANICAL=1)
endif()

# ======================================
# === exaNBody application framework ===
# ======================================
find_package(exaNBody)


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


