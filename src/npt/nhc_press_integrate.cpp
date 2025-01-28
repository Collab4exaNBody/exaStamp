#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_stream.h>
#include <exaStamp/npt/npt.h>
#include <exanb/core/domain.h>
#include <onika/physics/constants.h>
#include <exaStamp/compute/thermodynamic_state.h>

#include <iostream>
#include <string>

namespace exaStamp
{
  using namespace exanb;
  
  struct NHCPressIntegrateNode : public OperatorNode
  {
    
    ADD_SLOT( double                  , dt         , INPUT , REQUIRED );
    ADD_SLOT( long                    , timestep   , INPUT , REQUIRED );
    ADD_SLOT( long                    , simulation_end_iteration , INPUT , REQUIRED );
    ADD_SLOT( NPTContext              , npt_ctx    , INPUT_OUTPUT );
    ADD_SLOT( Domain                  , domain     , INPUT );    
    ADD_SLOT( ThermodynamicState      , thermodynamic_state , INPUT );

    inline void execute () override final
    {

      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      long natoms = sim_info.particle_count();
      
      int pdof;
      double expfac,factor_etap,kecurrent,lkt_press,nkt;
      double kt = npt_ctx->boltz * npt_ctx->t_target;
      
      if (npt_ctx->omega_mass_flag) {
        nkt = (natoms + 1) * kt;
        for (int i = 0; i < 3; i++)
          if (npt_ctx->p_flag[i])
            npt_ctx->omega_mass[i] = nkt/(npt_ctx->p_freq[i]*npt_ctx->p_freq[i]);
        if (npt_ctx->pstyle == "TRICLINIC") {
          for (int i = 3; i < 6; i++)
            if (npt_ctx->p_flag[i]) npt_ctx->omega_mass[i] = nkt/(npt_ctx->p_freq[i]*npt_ctx->p_freq[i]);
        }
      }

      if (npt_ctx->etap_mass_flag) {
      	if (npt_ctx->mpchain) {
      	  npt_ctx->etap_mass[0] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->p_freq_max*npt_ctx->p_freq_max);
      	  for (int ich = 1; ich < npt_ctx->mpchain; ich++) {
      	    npt_ctx->etap_mass[ich] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->p_freq_max*npt_ctx->p_freq_max);
          }
      	  for (int ich = 1; ich < npt_ctx->mpchain; ich++) {
      	    npt_ctx->etap_dotdot[ich] =
      	      (npt_ctx->etap_mass[ich-1]*npt_ctx->etap_dot[ich-1]*npt_ctx->etap_dot[ich-1] -
      	       npt_ctx->boltz * npt_ctx->t_target) / npt_ctx->etap_mass[ich];
          }
      	}
      }
      
      kecurrent = 0.0;
      pdof = 0;
      for (int i = 0; i < 3; i++)
      	if (npt_ctx->p_flag[i]) {
      	  kecurrent += npt_ctx->omega_mass[i]*npt_ctx->omega_dot[i]*npt_ctx->omega_dot[i];
      	  pdof++;
      	}
      
      if (npt_ctx->pstyle == "TRICLINIC") {
        for (int i = 3; i < 6; i++)
          if (npt_ctx->p_flag[i]) {
            kecurrent += npt_ctx->omega_mass[i]*npt_ctx->omega_dot[i]*npt_ctx->omega_dot[i];
            pdof++;
          }
      }
      
      if (npt_ctx->pstyle == "ISO") lkt_press = kt;
      else lkt_press = pdof * kt;
      
      npt_ctx->etap_dotdot[0] = (kecurrent - lkt_press)/npt_ctx->etap_mass[0];
      
      double ncfac = 1.0;

      for (int ich = npt_ctx->mpchain-1; ich > 0; ich--) {
        expfac = exp(-ncfac*npt_ctx->dt8*npt_ctx->etap_dot[ich+1]);
        npt_ctx->etap_dot[ich] *= expfac;
        npt_ctx->etap_dot[ich] += npt_ctx->etap_dotdot[ich] * ncfac*npt_ctx->dt4;
        npt_ctx->etap_dot[ich] *= npt_ctx->pdrag_factor;
        npt_ctx->etap_dot[ich] *= expfac;
      }

      expfac = exp(-ncfac*npt_ctx->dt8*npt_ctx->etap_dot[1]);
      npt_ctx->etap_dot[0] *= expfac;
      npt_ctx->etap_dot[0] += npt_ctx->etap_dotdot[0] * ncfac*npt_ctx->dt4;
      npt_ctx->etap_dot[0] *= npt_ctx->pdrag_factor;
      npt_ctx->etap_dot[0] *= expfac;
      
      for (int ich = 0; ich < npt_ctx->mpchain; ich++) {
        npt_ctx->etap[ich] += ncfac*npt_ctx->dthalf*npt_ctx->etap_dot[ich];
      }
      
      factor_etap = exp(-ncfac*npt_ctx->dthalf*npt_ctx->etap_dot[0]);
      for (int i = 0; i < 3; i++) {
        if (npt_ctx->p_flag[i]) npt_ctx->omega_dot[i] *= factor_etap;
      }

      if (npt_ctx->pstyle == "TRICLINIC") {
        for (int i = 3; i < 6; i++)
          if (npt_ctx->p_flag[i]) npt_ctx->omega_dot[i] *= factor_etap;
      }
      
      kecurrent = 0.0;
      for (int i = 0; i < 3; i++)
        if (npt_ctx->p_flag[i]) kecurrent += npt_ctx->omega_mass[i]*npt_ctx->omega_dot[i]*npt_ctx->omega_dot[i];

      if (npt_ctx->pstyle == "TRICLINIC") {
        for (int i = 3; i < 6; i++)
          if (npt_ctx->p_flag[i]) kecurrent += npt_ctx->omega_mass[i]*npt_ctx->omega_dot[i]*npt_ctx->omega_dot[i];
      }
      
      npt_ctx->etap_dotdot[0] = (kecurrent - lkt_press)/npt_ctx->etap_mass[0];
      npt_ctx->etap_dot[0] *= expfac;
      npt_ctx->etap_dot[0] += npt_ctx->etap_dotdot[0] * ncfac*npt_ctx->dt4;
      npt_ctx->etap_dot[0] *= expfac;
      
      for (int ich = 1; ich < npt_ctx->mpchain; ich++) {
        expfac = exp(-ncfac*npt_ctx->dt8*npt_ctx->etap_dot[ich+1]);
        npt_ctx->etap_dot[ich] *= expfac;
        npt_ctx->etap_dotdot[ich] = (npt_ctx->etap_mass[ich-1]*npt_ctx->etap_dot[ich-1]*npt_ctx->etap_dot[ich-1] - npt_ctx->boltz*npt_ctx->t_target) / npt_ctx->etap_mass[ich];
        npt_ctx->etap_dot[ich] += npt_ctx->etap_dotdot[ich] * ncfac*npt_ctx->dt4;
        npt_ctx->etap_dot[ich] *= expfac;
        ldbg << ich << std::endl;
      }
  
    }
  };
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(nhc_press_integrate)
  {
   OperatorNodeFactory::instance()->register_factory( "nhc_press_integrate", make_compatible_operator< NHCPressIntegrateNode > );
  }

}
