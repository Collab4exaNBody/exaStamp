#pragma once

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <exanb/core/compact_grid_pair_weights.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <iostream>
#include <type_traits>

// this allows for parallel compilation of templated operator for different variant functor : singe pair pot, multi-param pair pot, and rigid molecule pair pot
#define USTAMP_POTENTIAL_KIND_PAIR 0
#define USTAMP_POTENTIAL_KIND_RIGIDMOL 1
#define USTAMP_POTENTIAL_KIND_MULTI_PAIR 2


