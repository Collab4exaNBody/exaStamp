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
  
  struct SetupNPTNode : public OperatorNode
  {
    
    ADD_SLOT( double                  , dt         , INPUT , REQUIRED );
    ADD_SLOT( long                    , timestep   , INPUT , REQUIRED );
    ADD_SLOT( long                    , simulation_end_iteration , INPUT , REQUIRED );
    ADD_SLOT( NPTContext              , npt_ctx    , INPUT_OUTPUT );
    ADD_SLOT( Domain                  , domain     , INPUT );    
    ADD_SLOT( ThermodynamicState      , thermodynamic_state , INPUT );
    
    inline void execute () override final
    {
      static const double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann;
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      long natoms = sim_info.particle_count();
      npt_ctx->t_current = sim_info.temperature_scal() / natoms * conv_temperature;
      npt_ctx->tdof = 3 * natoms - 3;
      
      if ( npt_ctx->tstat_flag) {
        npt_ctx->update_target_T_KE(*timestep, *(simulation_end_iteration));
        npt_ctx->eta_mass[0] = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->t_freq*npt_ctx->t_freq);
        for (int ich = 1; ich < npt_ctx->mtchain; ich++) npt_ctx->eta_mass[ich] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->t_freq*npt_ctx->t_freq);
        for (int ich = 1; ich < npt_ctx->mtchain; ich++) npt_ctx->eta_dotdot[ich] = (npt_ctx->eta_mass[ich-1]*npt_ctx->eta_dot[ich-1]*npt_ctx->eta_dot[ich-1] - npt_ctx->boltz * npt_ctx->t_target) / npt_ctx->eta_mass[ich];
      }

      npt_ctx->t0 = 0.;
      if ( npt_ctx->pstat_flag ) {
        npt_ctx->p_hydro = 0.0;
        npt_ctx->update_target_P(*timestep, *(simulation_end_iteration));
        if (npt_ctx->deviatoric_flag) npt_ctx->update_sigma();
        
        double kt = npt_ctx->boltz * npt_ctx->t_target;
        double nkt = (natoms + 1) * kt;
        
        for (int i = 0; i < 3; i++)
          if (npt_ctx->p_flag[i])
            npt_ctx->omega_mass[i] = nkt/(npt_ctx->p_freq[i]*npt_ctx->p_freq[i]);
        for (int i = 3; i < 6; i++)
          if (npt_ctx->p_flag[i]) npt_ctx->omega_mass[i] = nkt/(npt_ctx->p_freq[i]*npt_ctx->p_freq[i]);
          
        if (npt_ctx->mpchain) {
          npt_ctx->etap_mass[0] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->p_freq_max*npt_ctx->p_freq_max);

          for (int ich = 1; ich < npt_ctx->mpchain; ich++)
            npt_ctx->etap_mass[ich] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->p_freq_max*npt_ctx->p_freq_max);
          for (int ich = 1; ich < npt_ctx->mpchain; ich++)
            npt_ctx->etap_dotdot[ich] = (npt_ctx->etap_mass[ich-1]*npt_ctx->etap_dot[ich-1]*npt_ctx->etap_dot[ich-1] - npt_ctx->boltz * npt_ctx->t_target) / npt_ctx->etap_mass[ich];
        }
      }
      
    }
  };
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "setup_npt", make_compatible_operator< SetupNPTNode > );
  }

}
