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

#include <iostream>
#include <string>

enum{NOBIAS,BIAS};
enum{NONE,XYZ,XY,YZ,XZ};
enum{ISO,ANISO,TRICLINIC};

namespace exaStamp
{
  using namespace exanb;

  struct InitNPTNode : public OperatorNode
  {
    static constexpr Mat3d all_one_matrix { 1.,1.,1., 1.,1.,1., 1.,1.,1. };

    ADD_SLOT( double , Tstart      , INPUT , REQUIRED);
    ADD_SLOT( double , Tend        , INPUT , OPTIONAL);
    ADD_SLOT( double , Tdamp       , INPUT , 0.1);
    ADD_SLOT( int    , tchain      , INPUT , 3);    
    ADD_SLOT( Mat3d  , Pstart      , INPUT , make_zero_matrix() );
    ADD_SLOT( Mat3d  , Pend        , INPUT , make_zero_matrix() );
    ADD_SLOT( double , Pdamp       , INPUT , 1.0);
    ADD_SLOT( int    , pchain      , INPUT , 3);
    ADD_SLOT( bool   , mtkcorr     , INPUT , false);
    ADD_SLOT( Domain , domain      , INPUT , REQUIRED);
    ADD_SLOT( double , dt          , INPUT , REQUIRED);
    ADD_SLOT( std::string , mode   , INPUT , "NVT");
    ADD_SLOT( NPTContext , npt_ctx , OUTPUT );
    ADD_SLOT( bool , pstat_flag    , OUTPUT );
    ADD_SLOT( bool , tstat_flag    , OUTPUT );
    ADD_SLOT( double , volume      , INPUT , REQUIRED );

    inline void execute () override final
    {
      double pascal_to_bar = 1.e-5;
      
      *npt_ctx = NPTContext { { *Tstart, *Tend, *Tdamp, (*Pstart)*pascal_to_bar, (*Pend)*pascal_to_bar, *Pdamp, *mode } };

      static const double conv_pressure = UnityConverterHelper::convert(1.e9, "kg.m^-1.s^-2");

      if ( *mode == "NVT" )
        {
          *pstat_flag = false;
          *tstat_flag = true;
          npt_ctx->tstat_flag = 1;
          npt_ctx->pstat_flag = 0;
          
          ldbg << "============== " << npt_ctx->m_config.m_mode <<" ==============" << std::endl;
          ldbg << "  Tstart     : " << npt_ctx->m_config.m_Tstart             << " K" << std::endl;
          ldbg << "  Tend       : " << npt_ctx->m_config.m_Tend               << " K" << std::endl;
          ldbg << "  Tdamp      : " << npt_ctx->m_config.m_Tdamp              << " ps"<< std::endl;
          ldbg << "=================================" <<std::endl<<std::endl;
          
        }
      else if ( *mode == "NPT" )
        {
          *pstat_flag = true;
          *tstat_flag = true;
          npt_ctx->tstat_flag = 1;
          npt_ctx->pstat_flag = 1;

          ldbg << "============== " << npt_ctx->m_config.m_mode <<" =============" << std::endl;
          ldbg << "  Tstart     : " << npt_ctx->m_config.m_Tstart             << " K" << std::endl;
          ldbg << "  Tend       : " << npt_ctx->m_config.m_Tend               << " K" << std::endl;
          ldbg << "  Tdamp      : " << npt_ctx->m_config.m_Tdamp              << " ps"<< std::endl;
          ldbg << "  Pstart     : " << npt_ctx->m_config.m_Pstart             << " Pa"<< std::endl;
          ldbg << "  Pend       : " << npt_ctx->m_config.m_Pend               << " Pa"<< std::endl;
          ldbg << "  Pdamp      : " << npt_ctx->m_config.m_Pdamp              << " ps"<< std::endl;
          ldbg << "================================" <<std::endl<<std::endl;
          
        }

      /* ------------------------------- NVT parameters ------------------------------- */
      npt_ctx->eta_mass_flag = 1;
      npt_ctx->t_start = npt_ctx->m_config.m_Tstart;
      npt_ctx->t_stop  = npt_ctx->m_config.m_Tend;
      if ( ! ( npt_ctx->m_config.m_Tdamp > 0.0 ) )
        {
          lerr << "Thermostat damping parameter must be > 0" << std::endl;
          std::abort();
        }
      npt_ctx->t_freq = 1./npt_ctx->m_config.m_Tdamp;
      npt_ctx->mtchain=*tchain;
      npt_ctx->eta = new double[npt_ctx->mtchain];
      npt_ctx->eta_dot = new double[npt_ctx->mtchain+1];
      npt_ctx->eta_dot[npt_ctx->mtchain] = 0.0;
      npt_ctx->eta_dotdot = new double[npt_ctx->mtchain];
      for (int ich = 0; ich < npt_ctx->mtchain; ich++) {
        npt_ctx->eta[ich] = npt_ctx->eta_dot[ich] = npt_ctx->eta_dotdot[ich] = 0.0;
      }
      npt_ctx->eta_mass = new double[npt_ctx->mtchain];
      int size_vector = 0;
      size_vector += 2*2*npt_ctx->mtchain;
      npt_ctx->tdrag_factor = 1.0;
      /* ------------------------------------------------------------------------------ */

      /* ------------------------------- NPT parameters ------------------------------- */
      npt_ctx->eta_mass_flag = 1;
      npt_ctx->etap_mass_flag = 0;
      
      npt_ctx->pcouple = NONE;
      npt_ctx->drag = 0.0;
      npt_ctx->allremap = 1;
      npt_ctx->dthalf = (*dt)/2.;
      npt_ctx->dt4 = (*dt)/4.;
      npt_ctx->dt8 = *dt/8.;
      npt_ctx->dto = npt_ctx->dthalf;
      npt_ctx->p_freq_max = 0.0;
      npt_ctx->boltz=8.617343e-5;
      npt_ctx->nktv2p=1.6021765e6;
      npt_ctx->mtk_flag = *mtkcorr;
      npt_ctx->deviatoric_flag = 1;
      npt_ctx->nreset_h0 = 0;

      Vec3d fpoint = 0.5 * domain->xform() * (domain->extent() - domain->origin());
      npt_ctx->fixedpoint[0] = fpoint.x;
      npt_ctx->fixedpoint[1] = fpoint.y;
      npt_ctx->fixedpoint[2] = fpoint.z;

      if ( ! ( npt_ctx->m_config.m_Pdamp > 0.0 ) )
        {
          lerr << "Barostat damping parameter must be > 0" << std::endl;
          std::abort();
        }
      for (int i = 0; i < 6; i++) {
        npt_ctx->p_freq[i] = 1./npt_ctx->m_config.m_Pdamp;
        npt_ctx->p_flag[i] = 1;
      }
      if (*pstat_flag) {
        npt_ctx->p_freq_max = std::max(npt_ctx->p_freq[0],npt_ctx->p_freq[1]);
        npt_ctx->p_freq_max = std::max(npt_ctx->p_freq_max,npt_ctx->p_freq[2]);
        npt_ctx->p_freq_max = std::max(npt_ctx->p_freq_max,npt_ctx->p_freq[3]);
        npt_ctx->p_freq_max = std::max(npt_ctx->p_freq_max,npt_ctx->p_freq[4]);
        npt_ctx->p_freq_max = std::max(npt_ctx->p_freq_max,npt_ctx->p_freq[5]);
        npt_ctx->pdrag_factor = 1.0 - ( (*dt) * npt_ctx->p_freq_max * npt_ctx->drag );
      }
      npt_ctx->mpchain=*pchain;

      if (npt_ctx->mpchain) {
        int ich;
        npt_ctx->etap = new double[npt_ctx->mpchain];
        npt_ctx->etap_dot = new double[npt_ctx->mpchain+1];
        npt_ctx->etap_dot[npt_ctx->mpchain] = 0.0;
        npt_ctx->etap_dotdot = new double[npt_ctx->mpchain];
        for (ich = 0; ich < npt_ctx->mpchain; ich++) {
          npt_ctx->etap[ich] = npt_ctx->etap_dot[ich] = npt_ctx->etap_dotdot[ich] = 0.0;
        }
        npt_ctx->etap_mass = new double[npt_ctx->mpchain];
        size_vector += 2*2*npt_ctx->mpchain;
      }
        
      if (npt_ctx->deviatoric_flag) size_vector += 1;	
      npt_ctx->vol0 = npt_ctx->t0 = 0.0;
      
      if (npt_ctx->pstat_flag) {
        npt_ctx->pdim = npt_ctx->p_flag[0] + npt_ctx->p_flag[1] + npt_ctx->p_flag[2];
        if (npt_ctx->vol0 == 0.0) {
          //	  if (dimension == 3) vol0 = domain->xprd * domain->yprd * domain->zprd;
          //	  else vol0 = domain->xprd * domain->yprd;
          npt_ctx->vol0 = *volume;
          Mat3d lattice = diag_matrix( domain->extent() - domain->origin() );
          //          Mat3d lot = transpose(domain->xform() * lattice);
          Mat3d lot = domain->xform() * lattice;
          Mat3d lot_inv = inverse( lot );
          // Might change h0_inv to a Mat3D later
          npt_ctx->h0_inv[0] = lot_inv.m11;
          npt_ctx->h0_inv[1] = lot_inv.m22;
          npt_ctx->h0_inv[2] = lot_inv.m33;
          npt_ctx->h0_inv[3] = lot_inv.m23;
          npt_ctx->h0_inv[4] = lot_inv.m13;
          npt_ctx->h0_inv[5] = lot_inv.m12;
        }

        double vnull = 0.0;
        ldbg << "Initial volume = " << npt_ctx->vol0 << std::endl;
        ldbg << "Initial reference cell = " << std::endl;
        ldbg << "| " << npt_ctx->h0_inv[0] << " | " << npt_ctx->h0_inv[5] << " | " << npt_ctx->h0_inv[4] << " | " << std::endl;
        ldbg << "| " << vnull << " | " << npt_ctx->h0_inv[1] << " | " << npt_ctx->h0_inv[3] << " | " << std::endl;
        ldbg << "| " << vnull << " | " << vnull << " | " << npt_ctx->h0_inv[2] << " | " << std::endl;
      }
            
      /* ------------------------------------------------------------------------------ */

      ldbg << "Initialization of " << *mode << " ensemble DONE" << std::endl;
    }
  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "init_npt", make_compatible_operator< InitNPTNode > );
  }

}
