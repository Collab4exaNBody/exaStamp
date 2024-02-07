#pragma once

#include "eam_alloy.h"

#define USTAMP_POTENTIAL_NAME eam_alloy
#define USTAMP_POTENTIAL_PARMS EamAlloyParameters
#define USTAMP_POTENTIAL_RCUT 1
#define USTAMP_POTENTIAL_EAM_PHI eam_alloy_phi
#define USTAMP_POTENTIAL_EAM_PHI_MM eam_alloy_phi_mm
#define USTAMP_POTENTIAL_EAM_RHO eam_alloy_rho
#define USTAMP_POTENTIAL_EAM_RHO_MM eam_alloy_rho_mm
#define USTAMP_POTENTIAL_EAM_EMB eam_alloy_fEmbed
#define USTAMP_POTENTIAL_EAM_EMB_MM eam_alloy_fEmbed_mm

// if defined, tells that multimat needs only one parameter set, and not one per atom type pair,
// also changes how phi and rho functions are called. if this is defined, the same correctly ordered type_i and tye_j are passed
// to phi and rho functions instead of EAMSpecyPairInfo
#define USTAMP_POTENTIAL_EAM_MM_UNIQUE_PARAMETER_SET 1

