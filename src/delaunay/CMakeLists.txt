#
# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements. See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership. The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.
#
# exaStamp delaunay (Delaunay tessellation of the particle system) is available only if Geogram
# is found. Set GEOGRAM_INSTALL_PREFIX to Geogram's installation prefix to enable it.

set(GEOGRAM_INSTALL_PREFIX "" CACHE PATH "Geogram installation prefix (empty: delaunay plugin disabled unless Geogram is found on the default search paths)")
if(GEOGRAM_INSTALL_PREFIX)
  list(APPEND CMAKE_MODULE_PATH "${GEOGRAM_INSTALL_PREFIX}/lib/cmake/modules")
endif()
find_package(Geogram QUIET)

if(Geogram_FOUND)
  message(STATUS "Geogram found: ${GEOGRAM_INCLUDE_DIR}")
  set(exaStampDelaunay_LINK_LIBRARIES exanbCore Geogram::geogram)
  onika_add_plugin(exaStampDelaunay ${CMAKE_CURRENT_SOURCE_DIR})
  # libgeogram.so isn't on the default runtime search path (not a system package) --
  # bake its location in so the plugin finds it without needing LD_LIBRARY_PATH set.
  if(GEOGRAM_INSTALL_PREFIX AND TARGET ${exaStampDelaunay_PLUGIN_LIB})
    set_target_properties(${exaStampDelaunay_PLUGIN_LIB} PROPERTIES
      BUILD_RPATH "${GEOGRAM_INSTALL_PREFIX}/lib"
      INSTALL_RPATH "${GEOGRAM_INSTALL_PREFIX}/lib"
    )
  endif()
else()
  message(STATUS "Geogram not found (set GEOGRAM_INSTALL_PREFIX) -- exaStampDelaunay plugin disabled")
endif()
