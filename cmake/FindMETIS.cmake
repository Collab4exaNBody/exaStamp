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

# FindMetis.cmake
#
#  METIS_FOUND - system has METIS
#  METIS_INCLUDE_DIRS - the METIS include directories
#  METIS_LIBRARIES - Link these to use METIS

if(NOT METIS_FOUND)

  if(NOT METIS_DIR)
    set(METIS_DIR "/usr/local/metis" CACHE PATH "METIS install path")
  else()
    if(NOT METIS_INCLUDE_DIRS)
        find_path(METIS_INCLUDE_DIRS parmetis.h PATHS ${METIS_DIR}/include)
    endif()
    if(NOT METIS_LIBRARIES)
      find_library(METIS_LIBRARY NAMES parmetis PATHS ${METIS_DIR}/lib)
      if(METIS_LIBRARY)
        set(METIS_LIBRARIES ${METIS_LIBRARY})
      endif()
    endif()
  endif()

  if(METIS_INCLUDE_DIRS AND METIS_LIBRARIES)
    set(METIS_FOUND TRUE)
  endif()

# message(STATUS "METIS_INCLUDE_DIRS=${METIS_INCLUDE_DIRS}")
# message(STATUS "METIS_LIBRARIES=${METIS_LIBRARIES}")

endif()
