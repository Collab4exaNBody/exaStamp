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

