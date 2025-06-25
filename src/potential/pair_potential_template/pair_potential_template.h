#pragma once

#include "potential.h"                  // definitons from the potenital's potential.h

#ifndef USTAMP_POTENTIAL_NAME
#error USTAMP_POTENTIAL_NAME must be defined to the name of the potential, without quotes
#endif

#ifndef USTAMP_POTENTIAL_PARAMS
#error USTAMP_POTENTIAL_PARAMS must be defined to the name of a structure containing potential configuration
#endif

#ifndef USTAMP_POTENTIAL_COMPUTE
#error USTAMP_POTENTIAL_COMPUTE must be defined to the name of a function
#endif

#ifndef USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT
#define USTAMP_POTENTIAL_PAIR_PARAMS_EXTRACT(p) (p)
#endif

#ifndef USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING
#define USTAMP_POTENTIAL_HANDLE_FORCE_WEIGHTING 0
#endif

#include <onika/cpp_utils.h>

#define USTAMP_POTENTIAL_CLASS USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,Potential)
#define USTAMP_POTENTIAL_REGISTER_FUNCTION USTAMP_CONCAT(register_potential_factory_,USTAMP_POTENTIAL_NAME)
#define USTAMP_POTENTIAL_STRING USTAMP_STR(USTAMP_POTENTIAL_NAME)

#define USTAMP_POTENTIAL_FORCE_OP_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,ComputeForceOperator)
#define USTAMP_POTENTIAL_FORCE_VIRIAL_OP_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,ComputeForceVirialOperator)
#define USTAMP_POTENTIAL_FORCE_SYMETRIC_OP_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,ComputeForceSymetricOperator)
#define USTAMP_POTENTIAL_FORCE_SYMETRIC_VIRIAL_OP_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,ComputeForceSymetricVirialOperator)

