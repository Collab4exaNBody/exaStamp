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
  
  struct NHOmegaDotNode : public OperatorNode
  {
    
    ADD_SLOT( double                  , dt         , INPUT , REQUIRED );
    ADD_SLOT( long                    , timestep   , INPUT , REQUIRED );
    ADD_SLOT( long                    , simulation_end_iteration , INPUT , REQUIRED );
    ADD_SLOT( NPTContext              , npt_ctx    , INPUT_OUTPUT );
    ADD_SLOT( Domain                  , domain     , INPUT );    
    ADD_SLOT( ThermodynamicState      , thermodynamic_state , INPUT );

/* ----------------------------------------------------------------------
   update omega_dot, omega
-----------------------------------------------------------------------*/

    inline void execute () override final
    {
      static constexpr double conv_pressure = 1.e4 * legacy_constant::atomicMass * 1e30;
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      double pascal_to_bar = 1.e-5;
      
      Mat3d ptensor_current = pascal_to_bar * sim_info.full_stress_tensor() * conv_pressure;
      
      //Mat3d ptensor_current = pascal_to_bar * sim_info.stress_tensor() * conv_pressure;      
      npt_ctx->p_current[0] = ptensor_current.m11;
      npt_ctx->p_current[1] = ptensor_current.m22;
      npt_ctx->p_current[2] = ptensor_current.m33;
      npt_ctx->p_current[3] = ptensor_current.m23;
      npt_ctx->p_current[4] = ptensor_current.m13;
      npt_ctx->p_current[5] = ptensor_current.m12;
      
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

      for (int i = 0; i < 3; i++)
        npt_ctx->p_hydro += npt_ctx->p_target[i];

      //      for (int i = 0; i < 3; i++) {
        //        npt_ctx->p_target[i] = npt_ctx->p_start[i] + delta * (npt_ctx->p_stop[i]-npt_ctx->p_start[i]);
      //      }
      npt_ctx->p_hydro /= npt_ctx->pdim;
        
      //      for (int i = 3; i < 6; i++)
      //        npt_ctx->p_target[i] = npt_ctx->p_start[i] + delta * (npt_ctx->p_stop[i]-npt_ctx->p_start[i]);      
      //      npt_ctx->p_hydro = trace_matrix(ptensor_current) / npt_ctx->pdim;
      //      npt_ctx->p_hydro = ( npt_ctx->p_target[0] + npt_ctx->p_target[1] + npt_ctx->p_target[2] )/ npt_ctx->pdim;
      ldbg << std::fixed;
      ldbg << std::setprecision(10);
      ldbg << "p_target0 = " << npt_ctx->p_target[0] << std::endl;
      ldbg << "p_target1 = " << npt_ctx->p_target[1] << std::endl;      
      ldbg << "p_target2 = " << npt_ctx->p_target[2] << std::endl;
      ldbg << "p_target3 = " << npt_ctx->p_target[3] << std::endl;
      ldbg << "p_target4 = " << npt_ctx->p_target[4] << std::endl;
      ldbg << "p_target5 = " << npt_ctx->p_target[5] << std::endl;      
      ldbg<< "p_hydro = " << npt_ctx->p_hydro << std::endl;
      // ------------------------- sigma computation --------------------------- //
      npt_ctx->sigma[0] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*((npt_ctx->p_target[0]-npt_ctx->p_hydro)*npt_ctx->h0_inv[0] + npt_ctx->p_target[5]*npt_ctx->h0_inv[5]+npt_ctx->p_target[4]*npt_ctx->h0_inv[4]) + npt_ctx->h0_inv[5]*(npt_ctx->p_target[5]*npt_ctx->h0_inv[0] + (npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[5]+npt_ctx->p_target[3]*npt_ctx->h0_inv[4]) + npt_ctx->h0_inv[4]*(npt_ctx->p_target[4]*npt_ctx->h0_inv[0]+npt_ctx->p_target[3]*npt_ctx->h0_inv[5] + (npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[4]));
      npt_ctx->sigma[1] = npt_ctx->vol0*(npt_ctx->h0_inv[1]*((npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[1] + npt_ctx->p_target[3]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[3]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[1] + (npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[3]));
      npt_ctx->sigma[2] = npt_ctx->vol0*(npt_ctx->h0_inv[2]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
      npt_ctx->sigma[3] = npt_ctx->vol0*(npt_ctx->h0_inv[1]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[3]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
      npt_ctx->sigma[4] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*(npt_ctx->p_target[4]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[5]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[2]) + npt_ctx->h0_inv[4]*((npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[2]));
      npt_ctx->sigma[5] = npt_ctx->vol0*(npt_ctx->h0_inv[0]*(npt_ctx->p_target[5]*npt_ctx->h0_inv[1]+npt_ctx->p_target[4]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[5]*((npt_ctx->p_target[1]-npt_ctx->p_hydro)*npt_ctx->h0_inv[1]+npt_ctx->p_target[3]*npt_ctx->h0_inv[3]) + npt_ctx->h0_inv[4]*(npt_ctx->p_target[3]*npt_ctx->h0_inv[1]+(npt_ctx->p_target[2]-npt_ctx->p_hydro)*npt_ctx->h0_inv[3]));
      // ------------------------- sigma computation --------------------------- //

      for (int i = 0; i < 6; i++) {
        ldbg << "sigma["<<i<<"] = " << npt_ctx->sigma[i] << std::endl;
      }

      //      ldbg << "\t\tp_current = " << ptensor_current << " GPa" << std::endl;      
      static constexpr double conv_temperature = 1.e4 * legacy_constant::atomicMass / legacy_constant::boltzmann;
      npt_ctx->t_current = sim_info.temperature_scal() / sim_info.particle_count() * conv_temperature;
      
      long natoms = sim_info.particle_count();
      double f_omega;
      double volume = sim_info.volume();
      ldbg << " Volume = " << volume << std::endl;
      // --------------- compute_deviatoric ------------- Move elsewhere later --------------- //
      //      Mat3d h = transpose(domain->xform() * diag_matrix( domain->extent() - domain->origin() ));
      Mat3d h = domain->xform() * diag_matrix( domain->extent() - domain->origin() );
      //      ldbg << "H matrix = " << h << std::endl;
      // Not sure about tilt values need to check w LAMMPS
      npt_ctx->fdev[0] =
        h.m11*(npt_ctx->sigma[0]*h.m11+npt_ctx->sigma[5]*h.m23+npt_ctx->sigma[4]*h.m13) +
        h.m23*(npt_ctx->sigma[5]*h.m11+npt_ctx->sigma[1]*h.m23+npt_ctx->sigma[3]*h.m13) +
        h.m13*(npt_ctx->sigma[4]*h.m11+npt_ctx->sigma[3]*h.m23+npt_ctx->sigma[2]*h.m13);
      npt_ctx->fdev[1] =
        h.m22*(              npt_ctx->sigma[1]*h.m22+npt_ctx->sigma[3]*h.m12) +
        h.m12*(              npt_ctx->sigma[3]*h.m22+npt_ctx->sigma[2]*h.m12);
      npt_ctx->fdev[2] =
        h.m33*(                            npt_ctx->sigma[2]*h.m33);
      npt_ctx->fdev[3] =
        h.m22*(                            npt_ctx->sigma[3]*h.m33) +
        h.m12*(                            npt_ctx->sigma[2]*h.m33);
      npt_ctx->fdev[4] =
        h.m11*(                            npt_ctx->sigma[4]*h.m33) +
        h.m23*(                            npt_ctx->sigma[3]*h.m33) +
        h.m13*(                            npt_ctx->sigma[2]*h.m33);
      npt_ctx->fdev[5] =
        h.m11*(              npt_ctx->sigma[5]*h.m22+npt_ctx->sigma[4]*h.m12) +
        h.m23*(              npt_ctx->sigma[1]*h.m22+npt_ctx->sigma[3]*h.m12) +
        h.m13*(              npt_ctx->sigma[3]*h.m22+npt_ctx->sigma[2]*h.m12);
      // --------------- compute_deviatoric ------------- Move elsewhere later --------------- //

      for (int i = 0; i < 6; i++) {
        ldbg << "fdev["<<i<<"] = " << npt_ctx->fdev[i] << std::endl;
      }
      
      // Martyna-Tobias-Klein correction, probably not needed
      npt_ctx->mtk_term1 = 0.0;
      // if (npt_ctx->mtk_flag) {
      // 	// if (npt_ctx->pstyle == ISO) {
      // 	//   npt_ctx->mtk_term1 = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_current;
      // 	//   npt_ctx->mtk_term1 /= npt_ctx->pdim * natoms;
      // 	// } else {
      //   Mat3d mvv_current = sim_info.ke_tensor();
      //   Vec3d mvv_vec = Vec3d{mvv_current.m11,mvv_current.m22,mvv_current.m33};
      //   //      	  for (int i = 0; i < 3; i++)
      //   //      	    if (npt_ctx->p_flag[i])
      //   npt_ctx->mtk_term1 += ( (mvv_vec.x+mvv_vec.y+mvv_vec.z) / npt_ctx->pdim * natoms);
      // }
      ldbg << "mtk_term1 = " << npt_ctx->mtk_term1 << std::endl;
      for (int i = 0; i < 6; i++) {
        ldbg << "p_current = " << npt_ctx->p_current[i] << std::endl;
      }
      for (int i = 0; i < 3; i++) {
        //      	if (npt_ctx->p_flag[i]) {
        f_omega = (npt_ctx->p_current[i]-npt_ctx->p_hydro)*volume /
          (npt_ctx->omega_mass[i] * npt_ctx->nktv2p) + npt_ctx->mtk_term1 / npt_ctx->omega_mass[i];
      if (npt_ctx->deviatoric_flag) f_omega -= npt_ctx->fdev[i]/(npt_ctx->omega_mass[i] * npt_ctx->nktv2p);
      npt_ctx->omega_dot[i] += f_omega*npt_ctx->dthalf;
      npt_ctx->omega_dot[i] *= npt_ctx->pdrag_factor;
      //      	}
      }
      ldbg << "f_omega = " << f_omega << std::endl;
      for (int i = 0; i < 3; i++) {
        ldbg << "omega_dot["<<i<<"] = " << npt_ctx->omega_dot[i] << std::endl;
      }      
      npt_ctx->mtk_term2 = 0.0;
      if (npt_ctx->mtk_flag) {
      	for (int i = 0; i < 3; i++)
      	  if (npt_ctx->p_flag[i])
      	    npt_ctx->mtk_term2 += npt_ctx->omega_dot[i];
      	if (npt_ctx->pdim > 0) npt_ctx->mtk_term2 /= npt_ctx->pdim * natoms;
      }
      ldbg << "mtk_term2 = " << npt_ctx->mtk_term2 << std::endl;

      //      if (npt_ctx->pstyle == TRICLINIC) {
      for (int i = 3; i < 6; i++) {
        if (npt_ctx->p_flag[i]) {
          f_omega = npt_ctx->p_current[i]*volume/(npt_ctx->omega_mass[i] * npt_ctx->nktv2p);
          if (npt_ctx->deviatoric_flag)
            f_omega -= npt_ctx->fdev[i]/(npt_ctx->omega_mass[i] * npt_ctx->nktv2p);
          npt_ctx->omega_dot[i] += f_omega*npt_ctx->dthalf;
          npt_ctx->omega_dot[i] *= npt_ctx->pdrag_factor;
        }
      }
        //      }
      ldbg << "f_omega = " << f_omega << std::endl;
      for (int i = 3; i < 6; i++) {
        ldbg << "omega_dot["<<i<<"] = " << npt_ctx->omega_dot[i] << std::endl;
      }      

      //      std::abort();
    }
  };
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "nh_omega_dot", make_compatible_operator< NHOmegaDotNode > );
  }

}
