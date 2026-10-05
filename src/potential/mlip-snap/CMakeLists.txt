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
# exaStamp snap is available only if contribs/md/snap has been built in exaNBody

option(EXASTAMP_MLIP_SNAP_BUILD "Build SNAP potential (snap + snaplegacy)" OFF)

if(EXASTAMP_MLIP_SNAP_BUILD AND XNB_BUILD_mdSnap)
#  message(STATUS "exaNBody has snap")
  set(exaStampMlipSnap_LINK_LIBRARIES mdSnap)
  onika_add_plugin(exaStampMlipSnap ${CMAKE_CURRENT_SOURCE_DIR})
endif()
