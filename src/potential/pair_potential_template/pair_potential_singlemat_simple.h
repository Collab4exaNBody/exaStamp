#pragma once

// simple pair pot
#define USTAMP_POTENTIAL_KIND USTAMP_POTENTIAL_KIND_PAIR
#define OPERATOR_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_compute_force)
#include "pair_potential_force_op_singlemat.h"
#define POTENTIAL_ADDITIONAL_FIELDS 

