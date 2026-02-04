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
# FindXSZoltan.cmake
#
#  ZOLTAN_FOUND - system has ZOLTAN
#  ZOLTAN_INCLUDE_DIRS - the ZOLTAN include directories
#  ZOLTAN_LIBRARIES - Link these to use ZOLTAN

if(NOT XSZOLTAN_FOUND)

  if(NOT ZOLTAN_DIR)
    #message(STATUS "ZOLTAN_ROOT='$ENV{ZOLTAN_ROOT}'")
    if(DEFINED ENV{ZOLTAN_ROOT})
      #message(STATUS "found ZOLTAN_ROOT env var, using as ZOLTAN_DIR default value")
      set(ZOLTAN_DIR_DEFAULT "$ENV{ZOLTAN_ROOT}")
    else()
      set(ZOLTAN_DIR_DEFAULT "/usr/local/zoltan")
    endif()
    set(ZOLTAN_DIR ${ZOLTAN_DIR_DEFAULT} CACHE PATH "ZOLTAN install path")
  else()
    if(NOT ZOLTAN_INCLUDE_DIRS)
        find_path(ZOLTAN_INCLUDE_DIRS zoltan.h PATHS ${ZOLTAN_DIR}/include)
    endif()
    if(NOT ZOLTAN_LIBRARIES)
      find_library(ZOLTAN_LIBRARY NAMES zoltan PATHS ${ZOLTAN_DIR}/lib)
      if(ZOLTAN_LIBRARY)
        set(ZOLTAN_LIBRARIES ${ZOLTAN_LIBRARY})
      endif()
    endif()
  endif()

  if(ZOLTAN_INCLUDE_DIRS AND ZOLTAN_LIBRARIES)
    set(XSZOLTAN_FOUND TRUE)
  endif()

 #message("ZOLTAN_INCLUDE_DIRS=${ZOLTAN_INCLUDE_DIRS}")
 #message("ZOLTAN_LIBRARIES=${ZOLTAN_LIBRARIES}")

endif()

