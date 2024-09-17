# ===================================
# === micro-stamp Snap potential  ===
# ===================================

add_executable(MakeSnapAutoGen snap_gen.cxx)

set(exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR ${PROJECT_BINARY_DIR}/include)
file(MAKE_DIRECTORY ${exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR})

set(exaStampSnap_GENERATED_SOURCES
  ${exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR}/compute_ui_jmax2.hxx
  ${exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR}/compute_ui_jmax3.hxx
  ${exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR}/compute_ui_jmax4.hxx
  ${exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR}/ymap_jmax2.hxx
  ${exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR}/ymap_jmax3.hxx
  ${exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR}/ymap_jmax4.hxx)

add_custom_target( SnapAutoGen
  COMMAND MakeSnapAutoGen ${exaStampSnap_GENERATED_SOURCES_OUTPUT_DIR}
  DEPENDS MakeSnapAutoGen
  BYPRODUCTS ${exaStampSnap_GENERATED_SOURCES} )

set(exaStampSnap_LINK_LIBRARIES exaStampPotential exanbParticleNeighbors)
xstamp_add_plugin(exaStampSnap ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(exaStampSnapPlugin SnapAutoGen)

