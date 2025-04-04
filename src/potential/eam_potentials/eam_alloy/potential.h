#pragma once

#include "eam_alloy.h"

#define USTAMP_POTENTIAL_NAME eam_alloy
#define USTAMP_POTENTIAL_PARMS EamAlloyParameters

#define USTAMP_POTENTIAL_ENABLE_CUDA 1
#define USTAMP_POTENTIAL_EAM_MM 1 // set if potential supports multimaterial and thus needs type numbers as additional arguments

#define USTAMP_POTENTIAL_EAM_MM_INIT_TYPES(p,nt,pe) (p).initialize_types_table( nt , pe )
#define USTAMP_POTENTIAL_EAM_MM_FORCE eam_alloy_mm_force
#define USTAMP_POTENTIAL_EAM_EMB eam_alloy_fEmbed
#define USTAMP_POTENTIAL_EAM_EMB_DERIV eam_alloy_fEmbed_deriv
#define USTAMP_POTENTIAL_EAM_RHO_NODERIV eam_alloy_rho_noderiv

