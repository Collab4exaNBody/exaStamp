#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/basic_types_stream.h>
#include <exaStamp/npt/npt.h>
#include <exanb/core/domain.h>
#include <exanb/core/physics_constants.h>
#include <exaStamp/compute/thermodynamic_state.h>

#include <iostream>
#include <string>

namespace exaStamp
{
  using namespace exanb;
  enum{ISO,ANISO,TRICLINIC};
  
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

      int ich,i,pdof;
      double expfac,factor_etap,kecurrent;
      double kt = npt_ctx->boltz * npt_ctx->t_target;
      ldbg << std::fixed;
      ldbg << std::setprecision(10);

      ldbg << "kt = " << kt << std::endl;      
      double lkt_press;
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      long natoms = sim_info.particle_count();
      
      if (npt_ctx->omega_mass_flag) {
        double nkt = (natoms + 1) * kt;
        for (i = 0; i < 6; i++)
          if (npt_ctx->p_flag[i])
            npt_ctx->omega_mass[i] = nkt/(npt_ctx->p_freq[i]*npt_ctx->p_freq[i]);
      }
      for (i = 0; i < 6; i++)
        ldbg << "omega_mass[="<<i<<"] = " << npt_ctx->omega_mass[i] << std::endl;      

      if (npt_ctx->etap_mass_flag) {
      	if (npt_ctx->mpchain) {
      	  npt_ctx->etap_mass[0] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->p_freq_max*npt_ctx->p_freq_max);
      	  for (ich = 1; ich < npt_ctx->mpchain; ich++) {
      	    npt_ctx->etap_mass[ich] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->p_freq_max*npt_ctx->p_freq_max);
            ldbg << "etap_mass["<<ich<<"] = " << npt_ctx->etap_mass[ich] << std::endl;      
          }
      	  for (ich = 1; ich < npt_ctx->mpchain; ich++) {
      	    npt_ctx->etap_dotdot[ich] =
      	      (npt_ctx->etap_mass[ich-1]*npt_ctx->etap_dot[ich-1]*npt_ctx->etap_dot[ich-1] -
      	       npt_ctx->boltz * npt_ctx->t_target) / npt_ctx->etap_mass[ich];
            ldbg << "etap_dotdot["<<ich<<"] = " << npt_ctx->etap_dotdot[ich] << std::endl;
          }
      	}
      }
      ldbg << "etap_mass_flag passed" << std::endl;
      ldbg << "etap_mass[0] = " << npt_ctx->etap_mass[0] << std::endl;      
      
      kecurrent = 0.0;
      pdof = 0;
      for (i = 0; i < 6; i++)
      	if (npt_ctx->p_flag[i]) {
      	  kecurrent += npt_ctx->omega_mass[i]*npt_ctx->omega_dot[i]*npt_ctx->omega_dot[i];
      	  pdof++;
      	}
      ldbg << "kecurrent = " << kecurrent << std::endl;
      
      //      if (npt_ctx->pstyle == ISO) lkt_press = kt;
      lkt_press = pdof * kt;
      npt_ctx->etap_dotdot[0] = (kecurrent - lkt_press)/npt_ctx->etap_mass[0];

      ldbg << "lkt_press = " << lkt_press << std::endl;
      ldbg << "etap_dotdot[0] = " << npt_ctx->etap_dotdot[0] << std::endl;      
      
      // NO SUBCYCLE FOR NOW IN NPT
      // MIGHT ADD LATER
      //     double ncfac = 1.0/npt_ctx->nc_pchain;
      // for (int iloop = 0; iloop < npt_ctx->nc_pchain; iloop++) {
      double ncfac = 1.0;

      for (ich = npt_ctx->mpchain-1; ich > 0; ich--) {
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

      ldbg << "expfac = " << expfac << std::endl;
      ldbg << "etap_dot[0] = " << npt_ctx->etap_dot[0] << std::endl;
      ldbg << "etap_dot[1] = " << npt_ctx->etap_dot[1] << std::endl;
      
      for (ich = 0; ich < npt_ctx->mpchain; ich++) {
        npt_ctx->etap[ich] += ncfac*npt_ctx->dthalf*npt_ctx->etap_dot[ich];
        ldbg << "etap["<<ich<<"] = " << npt_ctx->etap[ich] << std::endl;
      }
      
      factor_etap = exp(-ncfac*npt_ctx->dthalf*npt_ctx->etap_dot[0]);
      ldbg << "factor_etap = " << factor_etap << std::endl;
      for (i = 0; i < 6; i++) {
        if (npt_ctx->p_flag[i]) npt_ctx->omega_dot[i] *= factor_etap;
        ldbg << "omega_dot["<<i<<"] = " << npt_ctx->omega_dot[i] << std::endl;
      }
      
      kecurrent = 0.0;
      for (i = 0; i < 6; i++)
        if (npt_ctx->p_flag[i]) kecurrent += npt_ctx->omega_mass[i]*npt_ctx->omega_dot[i]*npt_ctx->omega_dot[i];
      ldbg << "kecurrent = " << kecurrent << std::endl;
      npt_ctx->etap_dotdot[0] = (kecurrent - lkt_press)/npt_ctx->etap_mass[0];

      npt_ctx->etap_dot[0] *= expfac;
      npt_ctx->etap_dot[0] += npt_ctx->etap_dotdot[0] * ncfac*npt_ctx->dt4;
      npt_ctx->etap_dot[0] *= expfac;

      ldbg << "etap_dotdot[0] = " << npt_ctx->etap_dotdot[0] << std::endl;
      ldbg << "etap_dot[0] = " << npt_ctx->etap_dot[0] << std::endl;
      
      for (ich = 1; ich < npt_ctx->mpchain; ich++) {
        expfac = exp(-ncfac*npt_ctx->dt8*npt_ctx->etap_dot[ich+1]);
        npt_ctx->etap_dot[ich] *= expfac;
        npt_ctx->etap_dotdot[ich] = (npt_ctx->etap_mass[ich-1]*npt_ctx->etap_dot[ich-1]*npt_ctx->etap_dot[ich-1] - npt_ctx->boltz*npt_ctx->t_target) / npt_ctx->etap_mass[ich];
        npt_ctx->etap_dot[ich] += npt_ctx->etap_dotdot[ich] * ncfac*npt_ctx->dt4;
        npt_ctx->etap_dot[ich] *= expfac;
        ldbg << ich << std::endl;
      }
      // SUBCYCLE LOOP ENDS
  //     }

      //      std::abort();
  
    }
  };
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "nhc_press_integrate", make_compatible_operator< NHCPressIntegrateNode > );
  }

}
