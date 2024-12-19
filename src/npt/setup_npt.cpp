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
      
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      long natoms = sim_info.particle_count();
      static constexpr double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann;
      static constexpr double conv_pressure = 1.e4 * legacy_constant::atomicMass * 1e30;      
      if ( npt_ctx->tstat_flag) {
        long start_at = 0;
        double delta = *timestep - start_at;
        if (delta != 0.0) delta /= *(simulation_end_iteration) - start_at;
        npt_ctx->t_target = npt_ctx->t_start + delta * (npt_ctx->t_stop-npt_ctx->t_start);
        npt_ctx->ke_target = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_target;

        npt_ctx->t_current = sim_info.temperature_scal() / natoms * conv_temperature;
        npt_ctx->tdof = 3 * natoms - 3;
        double kecurrent = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_current;
        npt_ctx->eta_mass[0] = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->t_freq*npt_ctx->t_freq);
        for (int ich = 1; ich < npt_ctx->mtchain; ich++) npt_ctx->eta_mass[ich] = npt_ctx->boltz * npt_ctx->t_target / (npt_ctx->t_freq*npt_ctx->t_freq);
        for (int ich = 1; ich < npt_ctx->mtchain; ich++) npt_ctx->eta_dotdot[ich] = (npt_ctx->eta_mass[ich-1]*npt_ctx->eta_dot[ich-1]*npt_ctx->eta_dot[ich-1] - npt_ctx->boltz * npt_ctx->t_target) / npt_ctx->eta_mass[ich];
      }

      npt_ctx->t0 = 0.;
      if ( npt_ctx->pstat_flag ) {
        long start_at=0;
        double delta = *timestep - start_at;
        if (delta != 0.0) delta /= *(simulation_end_iteration) - start_at;
        npt_ctx->p_hydro = 0.0;
        npt_ctx->p_target[0] = npt_ctx->m_config.m_Pstart.m11 + delta * (npt_ctx->m_config.m_Pend.m11-npt_ctx->m_config.m_Pstart.m11);
        npt_ctx->p_target[1] = npt_ctx->m_config.m_Pstart.m22 + delta * (npt_ctx->m_config.m_Pend.m22-npt_ctx->m_config.m_Pstart.m22);
        npt_ctx->p_target[2] = npt_ctx->m_config.m_Pstart.m33 + delta * (npt_ctx->m_config.m_Pend.m33-npt_ctx->m_config.m_Pstart.m33);
        npt_ctx->p_target[3] = npt_ctx->m_config.m_Pstart.m23 + delta * (npt_ctx->m_config.m_Pend.m23-npt_ctx->m_config.m_Pstart.m23);
        npt_ctx->p_target[4] = npt_ctx->m_config.m_Pstart.m13 + delta * (npt_ctx->m_config.m_Pend.m13-npt_ctx->m_config.m_Pstart.m13);
        npt_ctx->p_target[5] = npt_ctx->m_config.m_Pstart.m12 + delta * (npt_ctx->m_config.m_Pend.m12-npt_ctx->m_config.m_Pstart.m12);
        
        for (int i = 0; i < 3; i++) npt_ctx->p_hydro += npt_ctx->p_target[i];
        if (npt_ctx->pdim > 0) npt_ctx->p_hydro /= npt_ctx->pdim;
        
        if (npt_ctx->deviatoric_flag) {
          
          // generate upper-triangular half of
          // sigma = vol0*h0inv*(p_target-p_hydro)*h0inv^t
          // units of sigma are are PV/L^2 e.g. atm.A
          //
          // [ 0 5 4 ]   [ 0 5 4 ] [ 0 5 4 ] [ 0 - - ]
          // [ 5 1 3 ] = [ - 1 3 ] [ 5 1 3 ] [ 5 1 - ]
          // [ 4 3 2 ]   [ - - 2 ] [ 4 3 2 ] [ 4 3 2 ]
          
          npt_ctx->sigma[0] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*((npt_ctx->p_target[0]-npt_ctx->p_hydro)*npt_ctx->h0_inv[0] + npt_ctx->p_target[5]*npt_ctx->h0_inv[5]+npt_ctx->p_target[4]*npt_ctx->h0_inv[4]) + npt_ctx->h0_inv[5]*(npt_ctx->p_target[5]*npt_ctx->h0_inv[0] + (npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[5]+npt_ctx->p_target[3]*npt_ctx->h0_inv[4]) + npt_ctx->h0_inv[4]*(npt_ctx->p_target[4]*npt_ctx->h0_inv[0]+npt_ctx->p_target[3]*npt_ctx->h0_inv[5] + (npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[4]));
          npt_ctx->sigma[1] = npt_ctx->vol0*(npt_ctx->h0_inv[1]*((npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[1] + npt_ctx->p_target[3]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[3]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[1] + (npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[3]));
          npt_ctx->sigma[2] = npt_ctx->vol0*(npt_ctx->h0_inv[2]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
          npt_ctx->sigma[3] = npt_ctx->vol0*(npt_ctx->h0_inv[1]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[3]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
          npt_ctx->sigma[4] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*(npt_ctx->p_target[4]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[5]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[4]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
          npt_ctx->sigma[5] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*(npt_ctx->p_target[5]*npt_ctx->h0_inv[1]+npt_ctx->p_target[4]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[5]*((npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[1]+npt_ctx->p_target[3]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[4]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[1]+(npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[3]));
        }
        Mat3d ptensor_target { npt_ctx->p_target[0], npt_ctx->p_target[5], npt_ctx->p_target[4],
                               npt_ctx->p_target[5], npt_ctx->p_target[1], npt_ctx->p_target[3],
                               npt_ctx->p_target[4], npt_ctx->p_target[3], npt_ctx->p_target[2] };	    
        ldbg << "\t\tp_target = " << ptensor_target << " GPa" << std::endl;
      }

      // Get current pressure and send it to couple function (ADD LATER) for now : triclinic
      Mat3d ptensor_current = sim_info.full_stress_tensor() * conv_pressure;
      ldbg << "\t\tp_current = " << ptensor_current << " GPa" << std::endl;      
      npt_ctx->p_hydro = trace_matrix(ptensor_current) / npt_ctx->pdim;
      // ------------------------- sigma computation --------------------------- //
      npt_ctx->sigma[0] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*((npt_ctx->p_target[0]-npt_ctx->p_hydro)*npt_ctx->h0_inv[0] + npt_ctx->p_target[5]*npt_ctx->h0_inv[5]+npt_ctx->p_target[4]*npt_ctx->h0_inv[4]) + npt_ctx->h0_inv[5]*(npt_ctx->p_target[5]*npt_ctx->h0_inv[0] + (npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[5]+npt_ctx->p_target[3]*npt_ctx->h0_inv[4]) + npt_ctx->h0_inv[4]*(npt_ctx->p_target[4]*npt_ctx->h0_inv[0]+npt_ctx->p_target[3]*npt_ctx->h0_inv[5] + (npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[4]));
      npt_ctx->sigma[1] = npt_ctx->vol0*(npt_ctx->h0_inv[1]*((npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[1] + npt_ctx->p_target[3]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[3]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[1] + (npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[3]));
      npt_ctx->sigma[2] = npt_ctx->vol0*(npt_ctx->h0_inv[2]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
      npt_ctx->sigma[3] = npt_ctx->vol0*(npt_ctx->h0_inv[1]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[3]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
      npt_ctx->sigma[4] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*(npt_ctx->p_target[4]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[5]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[4]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
      npt_ctx->sigma[5] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*(npt_ctx->p_target[5]*npt_ctx->h0_inv[1]+npt_ctx->p_target[4]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[5]*((npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[1]+npt_ctx->p_target[3]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[4]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[1]+(npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[3]));
      // ------------------------- sigma computation --------------------------- //
      
      // masses and initial forces on barostat variables
      
      if (npt_ctx->pstat_flag) {
        double kt = npt_ctx->boltz * npt_ctx->t_target;
        double nkt = (natoms + 1) * kt;
        
        for (int i = 0; i < 3; i++)
          if (npt_ctx->p_flag[i])
            npt_ctx->omega_mass[i] = nkt/(npt_ctx->p_freq[i]*npt_ctx->p_freq[i]);
        for (int i = 3; i < 6; i++)
          if (npt_ctx->p_flag[i]) npt_ctx->omega_mass[i] = nkt/(npt_ctx->p_freq[i]*npt_ctx->p_freq[i]);
        
        // masses and initial forces on barostat thermostat variables
        
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
