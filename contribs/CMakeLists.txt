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

# ============================================================
# === Additional applications used togeteher with ExaStamp ===
# ============================================================


add_subdirectory(saru)

if(XSTAMP_ENABLE_SLIPLINK)
  option(XSTAMP_BUILD_SLIPLINK_REFAPP "Build SLink (slip-link model sample implementation)" OFF)
  if(XSTAMP_BUILD_SLIPLINK_REFAPP)
      add_subdirectory(sliplink)
  endif()
endif()

# =============================
# === legacy AMR code (Taz) ===
# =============================
option(EXASTAMP_BUILD_TAZ "Build Taz (AMR-ExaStamp)" OFF)
if(EXASTAMP_BUILD_TAZ)
    add_subdirectory(taz)
endif()

