# licensed to the apache software foundation (asf) under one
# or more contributor license agreements.  see the notice file
# distributed with this work for additional information
# regarding copyright ownership.  the asf licenses this file
# to you under the apache license, version 2.0 (the
# "license"); you may not use this file except in compliance
# with the license.  you may obtain a copy of the license at
# 
#   http://www.apache.org/licenses/license-2.0
# 
# unless required by applicable law or agreed to in writing,
# software distributed under the license is distributed on an
# "as is" basis, without warranties or conditions of any
# kind, either express or implied.  see the license for the
# specific language governing permissions and limitations
# under the license.

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
project(exaStamp VERSION 3.7 LANGUAGES C CXX)
set(EXASTAMP_VERSION ${CMAKE_PROJECT_VERSION}) 
message(STATUS "Project version : ${EXASTAMP_VERSION}")

# ===========================
# === CMake customization ===
# ===========================
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release RelWithDebInfo Debug)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
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
#  AddRegressionTestDir("${CMAKE_SOURCE_DIR}/data/regression")
  AddRegressionTestDir("${CMAKE_SOURCE_DIR}/data/regression/eam_potentials")
  AddRegressionTestDir("${CMAKE_SOURCE_DIR}/data/regression/pair_potentials")
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


