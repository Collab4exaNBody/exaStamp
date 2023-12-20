#pragma once

// rigid molecule
#define USTAMP_POTENTIAL_KIND USTAMP_POTENTIAL_KIND_RIGIDMOL
#define USTAMP_POTENTIAL_RIGIDMOL 1
#define OPERATOR_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_rigidmol_force)
#include "rigidmol_force_op_singlemat.h"
#define POTENTIAL_ADDITIONAL_FIELDS ,field::_orient,field::_couple

