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
      int ich;
      double expfac;

      const ThermodynamicState& sim_info = *(this->thermodynamic_state);

      // MAY REPLACE npt_ctx->boltz by legacy_constant::boltzmann
      // THIS WORKS
      static constexpr double boltzloc = 1.3806504e-23;
      // THIS ALSO WORKS
      //      static constexpr double boltzloc = legacy_constant::boltzmann;
      //      npt_ctx->boltz = legacy_constant::boltzmann;
      
      static constexpr double conv_temperature = 1.e4 * legacy_constant::atomicMass / boltzloc;
      npt_ctx->t_current = sim_info.temperature_scal() / sim_info.particle_count() * conv_temperature;
      double kecurrent = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_current;

      long start_at = 0;
      double delta = *timestep - start_at;
      if (delta != 0.0) delta /= *(simulation_end_iteration) - start_at;
      npt_ctx->t_target = npt_ctx->t_start + delta * (npt_ctx->t_stop-npt_ctx->t_start);
      npt_ctx->ke_target = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_target;
      ldbg << std::fixed;
      ldbg << std::setprecision(10);

      ldbg << "Target temperature      = " << npt_ctx->t_target << std::endl;
      ldbg << "Current temperature     = " << npt_ctx->t_current << std::endl;
      ldbg << "Target kinetic energy   = " << npt_ctx->ke_target << std::endl;
      ldbg << "Current kinetic energy  = " << kecurrent << std::endl;    
      ldbg << "tdof                    = " << npt_ctx->tdof << std::endl;    
      
      npt_ctx->eta_mass[0] = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->t_freq*npt_ctx->t_freq);
      ldbg << "eta_mass[0]=" << npt_ctx->eta_mass[0] << std::endl;
      
      for (ich = 1; ich < npt_ctx->mtchain; ich++)
        npt_ctx->eta_mass[ich] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->t_freq*npt_ctx->t_freq);

      if (npt_ctx->eta_mass[0] > 0.0)
      	npt_ctx->eta_dotdot[0] = (kecurrent - npt_ctx->ke_target)/npt_ctx->eta_mass[0];
      else npt_ctx->eta_dotdot[0] = 0.0;
      
      ldbg << "eta_dotdot[0]=" << npt_ctx->eta_dotdot[0] << std::endl;

      double ncfac = 1.0;
      
      // NO SUBCYCLE FOR NOW IN NVT
      // MIGHT ADD LATER
      //      double ncfac = 1.0/npt_ctx->nc_tchain;
      //  for (int iloop = 0; iloop < npt_ctx->nc_tchain; iloop++) {
      
      for (ich = npt_ctx->mtchain-1; ich > 0; ich--) {
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

      ldbg << "eta_dot[0]" << npt_ctx->eta_dot[0] << std::endl;
      
      npt_ctx->factor_eta = exp(-ncfac*npt_ctx->dthalf*npt_ctx->eta_dot[0]);
      ldbg << "factor_eta = " << npt_ctx->factor_eta << std::endl;

      //      npt_ctx->factor_eta=1.0;
      *vscale = npt_ctx->factor_eta * Vec3d{1.0, 1.0, 1.0};
      npt_ctx->t_current *= (npt_ctx->factor_eta*npt_ctx->factor_eta);
      kecurrent = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_current;

      ldbg << "t_current = " << npt_ctx->t_current << std::endl;
      ldbg << "kecurrent = " << kecurrent << std::endl;    
      
      if (npt_ctx->eta_mass[0] > 0.0)
        npt_ctx->eta_dotdot[0] = (kecurrent - npt_ctx->ke_target)/npt_ctx->eta_mass[0];
      else npt_ctx->eta_dotdot[0] = 0.0;

      ldbg << "eta_dotdot[0]" << npt_ctx->eta_dotdot[0] << std::endl;
      
      for (ich = 0; ich < npt_ctx->mtchain; ich++)
        npt_ctx->eta[ich] += ncfac*npt_ctx->dthalf*npt_ctx->eta_dot[ich];
      
      ldbg << "eta[0]" << npt_ctx->eta[0] << std::endl;
      
      npt_ctx->eta_dot[0] *= expfac;
      npt_ctx->eta_dot[0] += npt_ctx->eta_dotdot[0] * ncfac*npt_ctx->dt4;
      npt_ctx->eta_dot[0] *= expfac;

      ldbg << "eta_dot[0]" << npt_ctx->eta_dot[0] << std::endl;
      
      ldbg << "END nhc_temp_int" << std::endl;

      for (ich = 1; ich < npt_ctx->mtchain; ich++) {
        expfac = exp(-ncfac*npt_ctx->dt8*npt_ctx->eta_dot[ich+1]);
        npt_ctx->eta_dot[ich] *= expfac;
        npt_ctx->eta_dotdot[ich] = (npt_ctx->eta_mass[ich-1]*npt_ctx->eta_dot[ich-1]*npt_ctx->eta_dot[ich-1] - npt_ctx->boltz * npt_ctx->t_target)/npt_ctx->eta_mass[ich];
        npt_ctx->eta_dot[ich] += npt_ctx->eta_dotdot[ich] * ncfac*npt_ctx->dt4;
        npt_ctx->eta_dot[ich] *= expfac;
      }
      //  }
    }
  };
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "nhc_temp_integrate", make_compatible_operator< NHCTempIntegrateNode > );
  }

}
