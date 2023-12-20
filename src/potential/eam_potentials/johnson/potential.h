#pragma once

#include "johnson.h"

#define USTAMP_POTENTIAL_NAME johnson
#define USTAMP_POTENTIAL_PARMS EamJohnsonParameters
#define USTAMP_POTENTIAL_RCUT 1
#define USTAMP_POTENTIAL_EAM_PHI eam_johnson_phi
#define USTAMP_POTENTIAL_EAM_RHO eam_johnson_rho
#define USTAMP_POTENTIAL_EAM_EMB eam_johnson_fEmbed
#define USTAMP_POTENTIAL_ENABLE_CUDA 1

