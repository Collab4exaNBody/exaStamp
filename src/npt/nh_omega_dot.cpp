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
      
      const ThermodynamicState& sim_info = *(this->thermodynamic_state);
      long natoms = sim_info.particle_count();
      double volume = sim_info.volume();
      double f_omega;

      npt_ctx->p_hydro = 0.0;
      npt_ctx->update_target_P(*timestep, *(simulation_end_iteration));
      
      if (npt_ctx->deviatoric_flag) {
        npt_ctx->update_sigma();
        Mat3d h = domain->xform() * diag_matrix( domain->extent() - domain->origin() );
        npt_ctx->update_fdev(h);
      }
      
      // Martyna-Tobias-Klein correction, probably not needed
      npt_ctx->mtk_term1 = 0.0;
      // if (npt_ctx->mtk_flag) {
      // 	// if (npt_ctx->pstyle == "ISO") {
      // 	//   npt_ctx->mtk_term1 = npt_ctx->tdof * npt_ctx->boltz * npt_ctx->t_current;
      // 	//   npt_ctx->mtk_term1 /= npt_ctx->pdim * natoms;
      // 	// } else {
      //   Mat3d mvv_current = sim_info.ke_tensor();
      //   Vec3d mvv_vec = Vec3d{mvv_current.m11,mvv_current.m22,mvv_current.m33};
      //   //      	  for (int i = 0; i < 3; i++)
      //   //      	    if (npt_ctx->p_flag[i])
      //   npt_ctx->mtk_term1 += ( (mvv_vec.x+mvv_vec.y+mvv_vec.z) / npt_ctx->pdim * natoms);
      // }

      // ldbg << "######\n\n\n" << std::endl;
      // ldbg << "DEV FLAG  =" << npt_ctx->deviatoric_flag << std::endl;
      // ldbg << "VOLUME    =" << volume << std::endl;
      // ldbg << "mtk1    =" << npt_ctx->mtk_term1 << std::endl;
      // for (int i = 0; i < 3; i++) {
      //   ldbg << "fdev["<<i<<"] = " << npt_ctx->fdev[i] << std::endl;
      // }      
      
      for (int i = 0; i < 3; i++) {
        if (npt_ctx->p_flag[i]) {
          // ldbg << "p_current [" << i << "]=" << npt_ctx->p_current[i] << std::endl;
          // ldbg << "p_hydro =" << npt_ctx->p_hydro << std::endl;
          // ldbg << "omega_mass[" << i << "]=" << npt_ctx->omega_mass[i] << std::endl;          
          f_omega = (npt_ctx->p_current[i]-npt_ctx->p_hydro)*volume /
            (npt_ctx->omega_mass[i] * npt_ctx->nktv2p) + npt_ctx->mtk_term1 / npt_ctx->omega_mass[i];
          //          ldbg << "f_omega [" << i << "]=" << f_omega << std::endl;                
          if (npt_ctx->deviatoric_flag) 
            f_omega -= npt_ctx->fdev[i]/(npt_ctx->omega_mass[i] * npt_ctx->nktv2p);
          npt_ctx->omega_dot[i] += f_omega*npt_ctx->dthalf;
          npt_ctx->omega_dot[i] *= npt_ctx->pdrag_factor;
          //          ldbg << "f_omega [" << i << "]=" << f_omega << std::endl;      
        }
      }
      // for (int i = 0; i < 3; i++) {
      //   ldbg << "omega_dot["<<i<<"] = " << npt_ctx->omega_dot[i] << std::endl;
      // }      
      // for (int i = 0; i < 6; i++) {
      //   ldbg << "pflag["<<i<<"] = " << npt_ctx->p_flag[i] << std::endl;        
      // }      
      
      npt_ctx->mtk_term2 = 0.0;
      if (npt_ctx->mtk_flag) {
      	for (int i = 0; i < 3; i++)
      	  if (npt_ctx->p_flag[i])
      	    npt_ctx->mtk_term2 += npt_ctx->omega_dot[i];
      	if (npt_ctx->pdim > 0) npt_ctx->mtk_term2 /= npt_ctx->pdim * natoms;
      }
      
      if (npt_ctx->pstyle == "TRICLINIC") {
        for (int i = 3; i < 6; i++) {
          if (npt_ctx->p_flag[i]) {
            f_omega = npt_ctx->p_current[i]*volume/(npt_ctx->omega_mass[i] * npt_ctx->nktv2p);
            if (npt_ctx->deviatoric_flag)
              f_omega -= npt_ctx->fdev[i]/(npt_ctx->omega_mass[i] * npt_ctx->nktv2p);
            npt_ctx->omega_dot[i] += f_omega*npt_ctx->dthalf;
            npt_ctx->omega_dot[i] *= npt_ctx->pdrag_factor;
          }
        }
      }
      
    }
  };
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "nh_omega_dot", make_compatible_operator< NHOmegaDotNode > );
  }

}
