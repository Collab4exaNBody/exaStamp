#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_stream.h>
#include <exaStamp/npt/npt.h>
#include <exanb/core/domain.h>
#include <exanb/core/physics_constants.h>
#include <exaStamp/compute/thermodynamic_state.h>

#include <iostream>
#include <iomanip>
#include <string>

namespace exaStamp
{
  using namespace exanb;
  
  struct NHCTempIntegrateNode : public OperatorNode
  {
    ADD_SLOT( double                  , dt         , INPUT , REQUIRED );
    ADD_SLOT( long                    , timestep   , INPUT , REQUIRED );
    ADD_SLOT( long                    , simulation_end_iteration , INPUT , REQUIRED );    
    ADD_SLOT( NPTContext              , npt_ctx    , INPUT_OUTPUT );
    ADD_SLOT( Vec3d                   , vscale     , OUTPUT );
    ADD_SLOT( ThermodynamicState      , thermodynamic_state , INPUT );

    inline void execute () override final
    {

      static constexpr double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann;
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      long natoms = sim_info.particle_count();
      npt_ctx->t_current = sim_info.temperature_scal() / natoms * conv_temperature;
      double kecurrent = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_current;
      double expfac;

      npt_ctx->update_target_T_KE(*timestep, *(simulation_end_iteration));
      npt_ctx->eta_mass[0] = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->t_freq*npt_ctx->t_freq);
      
      for (int ich = 1; ich < npt_ctx->mtchain; ich++)
        npt_ctx->eta_mass[ich] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->t_freq*npt_ctx->t_freq);
      if (npt_ctx->eta_mass[0] > 0.0)
      	npt_ctx->eta_dotdot[0] = (kecurrent - npt_ctx->ke_target)/npt_ctx->eta_mass[0];
      else npt_ctx->eta_dotdot[0] = 0.0;
      
      double ncfac = 1.0;
      
      for (int ich = npt_ctx->mtchain-1; ich > 0; ich--) {
        expfac = exp(-ncfac*npt_ctx->dt8*npt_ctx->eta_dot[ich+1]);
        npt_ctx->eta_dot[ich] *= expfac;
        npt_ctx->eta_dot[ich] += npt_ctx->eta_dotdot[ich] * ncfac*npt_ctx->dt4;
        npt_ctx->eta_dot[ich] *= npt_ctx->tdrag_factor;
        npt_ctx->eta_dot[ich] *= expfac;
      }

      expfac = exp(-ncfac*npt_ctx->dt8*npt_ctx->eta_dot[1]);
      npt_ctx->eta_dot[0] *= expfac;
      npt_ctx->eta_dot[0] += npt_ctx->eta_dotdot[0] * ncfac*npt_ctx->dt4;
      npt_ctx->eta_dot[0] *= npt_ctx->tdrag_factor;
      npt_ctx->eta_dot[0] *= expfac;
      npt_ctx->factor_eta = exp(-ncfac*npt_ctx->dthalf*npt_ctx->eta_dot[0]);

      *vscale = npt_ctx->factor_eta * Vec3d{1.0, 1.0, 1.0};
      npt_ctx->t_current *= (npt_ctx->factor_eta*npt_ctx->factor_eta);
      kecurrent = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_current;
      
      if (npt_ctx->eta_mass[0] > 0.0)
        npt_ctx->eta_dotdot[0] = (kecurrent - npt_ctx->ke_target)/npt_ctx->eta_mass[0];
      else npt_ctx->eta_dotdot[0] = 0.0;
      
      for (int ich = 0; ich < npt_ctx->mtchain; ich++)
        npt_ctx->eta[ich] += ncfac*npt_ctx->dthalf*npt_ctx->eta_dot[ich];
      
      npt_ctx->eta_dot[0] *= expfac;
      npt_ctx->eta_dot[0] += npt_ctx->eta_dotdot[0] * ncfac*npt_ctx->dt4;
      npt_ctx->eta_dot[0] *= expfac;

      for (int ich = 1; ich < npt_ctx->mtchain; ich++) {
        expfac = exp(-ncfac*npt_ctx->dt8*npt_ctx->eta_dot[ich+1]);
        npt_ctx->eta_dot[ich] *= expfac;
        npt_ctx->eta_dotdot[ich] = (npt_ctx->eta_mass[ich-1]*npt_ctx->eta_dot[ich-1]*npt_ctx->eta_dot[ich-1] - npt_ctx->boltz * npt_ctx->t_target)/npt_ctx->eta_mass[ich];
        npt_ctx->eta_dot[ich] += npt_ctx->eta_dotdot[ich] * ncfac*npt_ctx->dt4;
        npt_ctx->eta_dot[ich] *= expfac;
      }

    }
  };
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(nhc_temp_integrate)
  {
   OperatorNodeFactory::instance()->register_factory( "nhc_temp_integrate", make_compatible_operator< NHCTempIntegrateNode > );
  }

}
