# ===================================
# === micro-stamp meam potential  ===
# ===================================

#TODO: par la suite, changer la racine de n2p2 vers un lien plus "général"
set(XSTAMP_BUILD_exaStampN2P2_DEFAULT OFF) # disabled by default
set(N2P2_ROOT_DIR /usr/local/n2p2 CACHE STRING "N2P2 installation directory")

# ${MPI_CXX_LIBRARIES}
set(microStampN2P2_INCLUDE_DIRS ${N2P2_ROOT_DIR}/include)
set(microStampN2P2_LINK_LIBRARIES exaStampPotential ${N2P2_ROOT_DIR}/lib/libnnp.so ${N2P2_ROOT_DIR}/lib/libnnpif.so)
xstamp_add_plugin(exaStampN2P2 ${CMAKE_CURRENT_SOURCE_DIR})

