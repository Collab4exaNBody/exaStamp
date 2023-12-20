# ===============================================
# === slink: a slip link model implementation ===
# ===============================================

file(GLOB SLINK_SRCS *.cpp)
add_executable(slink ${SLINK_SRCS})
install(TARGETS slink DESTINATION bin)

