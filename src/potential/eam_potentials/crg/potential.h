#pragma once

// Cooper Rushton grimes potential
// M.W.D. Cooper, M.J.D Rushton, R.W. Grimes, J. Phys. Condens. Matter 26 (2014) 105401.

#include "crg.h"

#define USTAMP_POTENTIAL_NAME crg
#define USTAMP_POTENTIAL_PARMS EamCRGParameters
#define USTAMP_POTENTIAL_RCUT 1
#define USTAMP_POTENTIAL_EAM_PHI crg_phi
#define USTAMP_POTENTIAL_EAM_RHO crg_rho
#define USTAMP_POTENTIAL_EAM_PHI_MM crg_phi_mm
#define USTAMP_POTENTIAL_EAM_RHO_MM crg_rho_mm
#define USTAMP_POTENTIAL_EAM_EMB crg_fEmbed
//#define USTAMP_POTENTIAL_ENABLE_CUDA 1


