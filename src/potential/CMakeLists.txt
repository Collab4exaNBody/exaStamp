# =================================
# === micro-stamp core library  ===
# =================================

# compile options specific to pair potentials
option(XSTAMP_POTENTIAL_FJEXEC "Fork-Join specialized implementation for pair potentials" ON)
if(XSTAMP_POTENTIAL_FJEXEC)
  set(USTAMP_COMPILE_DEFINITIONS ${USTAMP_COMPILE_DEFINITIONS} -DXSTAMP_POTENTIAL_FJEXEC=1)
endif()

include(PairPotentialTemplate)
include(EAMPotentialTemplate)

# Pair potentials
add_subdirectory(pair_potentials)
add_subdirectory(multimat_pair_potential)

# EAM potentials
add_subdirectory(eam)
add_subdirectory(eam_potentials)

add_subdirectory(LCHBOP)
add_subdirectory(meam)
add_subdirectory(reaction_field)
add_subdirectory(coul_dsf)
add_subdirectory(coul_wolf)
add_subdirectory(ewald)
add_subdirectory(snap)
add_subdirectory(snaplegacy)
add_subdirectory(snaplmp)
add_subdirectory(range_potential)
add_subdirectory(n2p2)


GenerateAllPairPotentials()
GenerateAllEAMPotentials()

