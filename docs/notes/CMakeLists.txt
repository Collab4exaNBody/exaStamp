# ==================================================
# === Documentation compilation and installation ===
# ==================================================

option(EXASTAMP_BUILD_LATEX_DOC "Build LateX documentation" OFF)

if(EXASTAMP_BUILD_LATEX_DOC)
  include(UseLATEX)
  add_subdirectory(reference_manual)
  add_subdirectory(developer_guide)
endif()

