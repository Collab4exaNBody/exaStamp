#pragma once

// multi-param pair pot
#define USTAMP_POTENTIAL_KIND USTAMP_POTENTIAL_KIND_MULTI_PAIR
#define USTAMP_POTENTIAL_MULTI_PARAM 1
#define OPERATOR_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_multi_force)
#include "pair_potential_force_op_multiparam.h"
#define POTENTIAL_ADDITIONAL_FIELDS ,field::_type

