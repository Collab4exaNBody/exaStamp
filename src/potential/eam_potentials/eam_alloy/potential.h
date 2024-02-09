#pragma once

#include "eam_alloy.h"

#define USTAMP_POTENTIAL_NAME eam_alloy
#define USTAMP_POTENTIAL_PARMS EamAlloyParameters
#define USTAMP_POTENTIAL_EAM_PHI eam_alloy_phi
#define USTAMP_POTENTIAL_EAM_RHO eam_alloy_rho
#define USTAMP_POTENTIAL_EAM_EMB eam_alloy_fEmbed

#define USTAMP_POTENTIAL_ENABLE_CUDA 1
#define USTAMP_POTENTIAL_EAM_MM 1 // set if potential supports multimaterial and thus needs type numbers as additional arguements
#define USTAMP_POTENTIAL_EAM_MM_INIT_TYPES(p,nt,pe) (p).initialize_types_table( nt , pe )

