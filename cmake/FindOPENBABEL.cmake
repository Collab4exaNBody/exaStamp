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

# FindOpenBabel.cmake
# Copyright (C) 2012-2013 NextMove Software
# Try to find Open Babel headers and libraries
# Defines:
#
#  OPENBABEL_FOUND - system has Open Babel
#  OPENBABEL_INCLUDE_DIRS - the Open Babel include directories
#  OPENBABEL_LIBRARIES - Link these to use Open Babel

if(NOT OPENBABEL_FOUND)

  if(NOT OPENBABEL_DIR)
    set(OPENBABEL_DIR "/usr/local/openbabel" CACHE PATH "OpenBabel install path")
  endif()

  #message("OPENBABEL_DIR=${OPENBABEL_DIR}")

  if(NOT OPENBABEL_INCLUDE_DIRS)
	find_path(OPENBABEL_INCLUDE_DIRS openbabel/babelconfig.h PATHS ${OPENBABEL_DIR}/include/openbabel-2.0)
  endif()

  if(NOT OPENBABEL_LIBRARIES)
	find_library(OPENBABEL_LIBRARIES NAMES openbabel PATHS ${OPENBABEL_DIR}/lib NO_DEFAULT_PATH)
  endif()

endif()

if(OPENBABEL_INCLUDE_DIRS AND OPENBABEL_LIBRARIES)
  set(OPENBABEL_FOUND TRUE)
  #message("OPENBABEL_INCLUDE_DIRS=${OPENBABEL_INCLUDE_DIRS}")
  #message("OPENBABEL_LIBRARIES=${OPENBABEL_LIBRARIES}")
endif()

