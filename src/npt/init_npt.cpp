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
    
    ADD_SLOT( double , Pstart      , INPUT , OPTIONAL);
    ADD_SLOT( double , Pend        , INPUT , OPTIONAL);
    ADD_SLOT( double , Pdamp       , INPUT , OPTIONAL);

    ADD_SLOT( double , Pxxstart      , INPUT , OPTIONAL);
    ADD_SLOT( double , Pxxend        , INPUT , OPTIONAL);
    ADD_SLOT( double , Pxxdamp       , INPUT , OPTIONAL);

    ADD_SLOT( double , Pyystart      , INPUT , OPTIONAL);
    ADD_SLOT( double , Pyyend        , INPUT , OPTIONAL);
    ADD_SLOT( double , Pyydamp       , INPUT , OPTIONAL);

    ADD_SLOT( double , Pzzstart      , INPUT , OPTIONAL);
    ADD_SLOT( double , Pzzend        , INPUT , OPTIONAL);
    ADD_SLOT( double , Pzzdamp       , INPUT , OPTIONAL);

    ADD_SLOT( double , Pyzstart      , INPUT , OPTIONAL);
    ADD_SLOT( double , Pyzend        , INPUT , OPTIONAL);
    ADD_SLOT( double , Pyzdamp       , INPUT , OPTIONAL);

    ADD_SLOT( double , Pxzstart      , INPUT , OPTIONAL);
    ADD_SLOT( double , Pxzend        , INPUT , OPTIONAL);
    ADD_SLOT( double , Pxzdamp       , INPUT , OPTIONAL);

    ADD_SLOT( double , Pxystart      , INPUT , OPTIONAL);
    ADD_SLOT( double , Pxyend        , INPUT , OPTIONAL);
    ADD_SLOT( double , Pxydamp       , INPUT , OPTIONAL);
    
    ADD_SLOT( int    , pchain      , INPUT , 3);
    ADD_SLOT( std::string, Pmode   , INPUT , "OTHER");
    ADD_SLOT( std::string, Pcouple , INPUT , "NONE");
    ADD_SLOT( bool   , Scalexy     , INPUT, OPTIONAL);
    ADD_SLOT( bool   , Scalexz     , INPUT, OPTIONAL);
    ADD_SLOT( bool   , Scaleyz     , INPUT, OPTIONAL);
    ADD_SLOT( double , Pdrag       , INPUT , 0.);
    ADD_SLOT( bool   , mtkcorr     , INPUT , false);
    ADD_SLOT( std::string , algo   , INPUT , "NVT");
    
    ADD_SLOT( NPTContext , npt_ctx , OUTPUT );
    ADD_SLOT( bool , pstat_flag    , OUTPUT );
    ADD_SLOT( bool , tstat_flag    , OUTPUT );
    ADD_SLOT( Domain , domain      , INPUT , REQUIRED);
    ADD_SLOT( double , dt          , INPUT , REQUIRED);

    inline void execute () override final
    {
      double pascal_to_bar = 1.e-5;      
      *npt_ctx = NPTContext { *algo, *Pmode};

      npt_ctx->boltz = UnityConverterHelper::convert(1, "J/eV") * legacy_constant::boltzmann;
      npt_ctx->nktv2p = UnityConverterHelper::convert(1, "m^3") * legacy_constant::elementaryCharge * pascal_to_bar;      
      npt_ctx->dthalf = (*dt)/2.;
      npt_ctx->dt4 = (*dt)/4.;
      npt_ctx->dt8 = *dt/8.;
      npt_ctx->dto = npt_ctx->dthalf;
      npt_ctx->p_freq_max = 0.0;

      /* default values */
      npt_ctx->pcouple = *Pcouple;
      npt_ctx->drag = *Pdrag;
      npt_ctx->mtchain = *tchain;
      npt_ctx->mpchain = *pchain;
      npt_ctx->nc_tchain = npt_ctx->nc_pchain = 1;
      if (*mtkcorr) {
        npt_ctx->mtk_flag = 1;
      } else {
        npt_ctx->mtk_flag = 0;
      }
      npt_ctx->deviatoric_flag = 0;
      npt_ctx->eta_mass_flag = 1;
      npt_ctx->omega_mass_flag = 0;      
      npt_ctx->etap_mass_flag = 0;
      npt_ctx->start_it = 0;

      if (domain->periodic_boundary_y()) npt_ctx->scalexy = 1;
      if (domain->periodic_boundary_z()) {
        npt_ctx->scalexz = 1;
        npt_ctx->scaleyz = 1;
      }
      
      npt_ctx->scaleyz = npt_ctx->scalexz = npt_ctx->scalexy = 1;
      Mat3d H0 = domain->xform() * diag_matrix( domain->extent() - domain->origin() );
      Mat3d H0_inv = inverse( H0 );

      // set fixed-point to default = center of cell
      Vec3d fpoint = 0.5 * domain->xform() * (domain->extent() - domain->origin());
      npt_ctx->fixedpoint[0] = fpoint.x;
      npt_ctx->fixedpoint[1] = fpoint.y;
      npt_ctx->fixedpoint[2] = fpoint.z;

      npt_ctx->tstat_flag = 0;
      for (int i = 0; i < 6; i++) {
        npt_ctx->p_start[i] = npt_ctx->p_stop[i] = npt_ctx->p_period[i] = npt_ctx->p_target[i] = 0.0;
        npt_ctx->p_flag[i] = 0;
      }
      npt_ctx->t_freq = 0.0;
      npt_ctx->p_freq[0] = npt_ctx->p_freq[1] = npt_ctx->p_freq[2] = npt_ctx->p_freq[3] = npt_ctx->p_freq[4] = npt_ctx->p_freq[5] = 0.0;
      int size_vector = 0;

      if ( npt_ctx->algo == "NVT" )
        {
          *pstat_flag = false;
          *tstat_flag = true;
          npt_ctx->tstat_flag = 1;
          npt_ctx->pstat_flag = 0;          
        }
      else if ( npt_ctx->algo == "NPT" )
        {
          *pstat_flag = true;
          *tstat_flag = true;
          npt_ctx->tstat_flag = 1;
          npt_ctx->pstat_flag = 1;
        }
      
      /* ------------------------------- NVT parameters ------------------------------- */
      if (npt_ctx->tstat_flag) {
        npt_ctx->t_start = *Tstart;
        if ( Tend.has_value() )
          npt_ctx->t_stop  = *Tend;
        else
          npt_ctx->t_stop  = npt_ctx->t_start;
        if ( ( ! *Tdamp ) > 0.0 )
          {
            lerr << "Thermostat damping parameter must be > 0" << std::endl;
            std::abort();
          }
          npt_ctx->t_freq = 1./(*Tdamp);
        npt_ctx->mtchain=*tchain;
        npt_ctx->eta = new double[npt_ctx->mtchain];
        npt_ctx->eta_dot = new double[npt_ctx->mtchain+1];
        npt_ctx->eta_dot[npt_ctx->mtchain] = 0.0;
        npt_ctx->eta_dotdot = new double[npt_ctx->mtchain];
        for (int ich = 0; ich < npt_ctx->mtchain; ich++) {
          npt_ctx->eta[ich] = npt_ctx->eta_dot[ich] = npt_ctx->eta_dotdot[ich] = 0.0;
        }
        npt_ctx->eta_mass = new double[npt_ctx->mtchain];
        size_vector += 2*2*npt_ctx->mtchain;
        npt_ctx->tdrag_factor = 1.0;
      }
      /* ------------------------------------------------------------------------------ */
      
      /* ------------------------------- NPT parameters ------------------------------- */      
      if (npt_ctx->pstat_flag) {
        npt_ctx->pstyle = "ISO";
        npt_ctx->pmode = *Pmode;
        if (npt_ctx->pmode == "ISO") {
          npt_ctx->pstyle = "ISO";
          npt_ctx->pcouple = "XYZ";
          npt_ctx->p_start[0] = npt_ctx->p_start[1] = npt_ctx->p_start[2] = (*Pstart)*pascal_to_bar;
          npt_ctx->p_stop[0] = npt_ctx->p_stop[1] = npt_ctx->p_stop[2] = (*Pend)*pascal_to_bar;
          npt_ctx->p_flag[0] = npt_ctx->p_flag[1] = npt_ctx->p_flag[2] = 1;
          npt_ctx->p_flag[3] = npt_ctx->p_flag[4] = npt_ctx->p_flag[5] = 0;
          if ( ( ! *Pdamp )> 0.0 )
            {
              lerr << "Barostat damping parameter must be > 0" << std::endl;
              std::abort();
            }          
          npt_ctx->p_period[0] = npt_ctx->p_period[1] = npt_ctx->p_period[2] = *Pdamp;
        }
        else if (npt_ctx->pmode == "ANISO") {
          npt_ctx->pstyle = "ANISO";          
          npt_ctx->pcouple = "NONE";
          npt_ctx->p_start[0] = npt_ctx->p_start[1] = npt_ctx->p_start[2] = (*Pstart)*pascal_to_bar;
          npt_ctx->p_stop[0] = npt_ctx->p_stop[1] = npt_ctx->p_stop[2] = (*Pend)*pascal_to_bar;
          npt_ctx->p_flag[0] = npt_ctx->p_flag[1] = npt_ctx->p_flag[2] = 1;
          npt_ctx->p_flag[3] = npt_ctx->p_flag[4] = npt_ctx->p_flag[5] = 0;
          if ( ( ! *Pdamp )> 0.0 )
            {
              lerr << "Barostat damping parameter must be > 0" << std::endl;
              std::abort();
            }          
          npt_ctx->p_period[0] = npt_ctx->p_period[1] = npt_ctx->p_period[2] = *Pdamp;
        }
        else if (npt_ctx->pmode == "TRI") {
          npt_ctx->pstyle = "TRICLINIC";
          npt_ctx->pcouple = "NONE";
          npt_ctx->scalexy = npt_ctx->scalexz = npt_ctx->scaleyz = 0;
          npt_ctx->p_start[0] = npt_ctx->p_start[1] = npt_ctx->p_start[2] = (*Pstart)*pascal_to_bar;
          npt_ctx->p_stop[0] = npt_ctx->p_stop[1] = npt_ctx->p_stop[2] = (*Pend)*pascal_to_bar;
          npt_ctx->p_flag[0] = npt_ctx->p_flag[1] = npt_ctx->p_flag[2] = 1;
          npt_ctx->p_start[3] = npt_ctx->p_start[4] = npt_ctx->p_start[5] = 0.0;
          npt_ctx->p_stop[3] = npt_ctx->p_stop[4] = npt_ctx->p_stop[5] = 0.0;
          npt_ctx->p_flag[3] = npt_ctx->p_flag[4] = npt_ctx->p_flag[5] = 1;
          if ( ( ! *Pdamp )> 0.0 )
            {
              lerr << "Barostat damping parameter must be > 0" << std::endl;
              std::abort();
            }
          npt_ctx->p_period[0] = npt_ctx->p_period[1] = npt_ctx->p_period[2] = *Pdamp;
          npt_ctx->p_period[3] = npt_ctx->p_period[4] = npt_ctx->p_period[5] = *Pdamp;
        }
        else {
          if(Pxxstart.has_value()) {
            npt_ctx->p_start[0] = (*Pxxstart)*pascal_to_bar;
            if(Pxxend.has_value()) {          
              npt_ctx->p_stop[0] = (*Pxxend)*pascal_to_bar;
            } else {
              npt_ctx->p_stop[0] = npt_ctx->p_start[0];
            }
            if(Pxxdamp.has_value()) {
              npt_ctx->p_period[0] = *Pxxdamp;
            } else {
              lerr << "Barostat damping parameter Pxxdamp must be > 0" << std::endl;
              std::abort();
            }
          npt_ctx->p_flag[0] = 1;
          npt_ctx->deviatoric_flag = 1;
          }
          if(Pyystart.has_value()) {
            npt_ctx->p_start[1] = (*Pyystart)*pascal_to_bar;
            if(Pyyend.has_value()) {          
              npt_ctx->p_stop[1] = (*Pyyend)*pascal_to_bar;
            } else {
              npt_ctx->p_stop[1] = npt_ctx->p_start[1];
            }
            if(Pyydamp.has_value()) {
              npt_ctx->p_period[1] = *Pyydamp;
            } else {
              lerr << "Barostat damping parameter Pyydamp must be > 0" << std::endl;
              std::abort();
            }
          npt_ctx->p_flag[1] = 1;
          npt_ctx->deviatoric_flag = 1;
          }
          if(Pzzstart.has_value()) {
            npt_ctx->p_start[2] = (*Pzzstart)*pascal_to_bar;
            if(Pzzend.has_value()) {          
              npt_ctx->p_stop[2] = (*Pzzend)*pascal_to_bar;
            } else {
              npt_ctx->p_stop[2] = npt_ctx->p_start[2];
            }
            if(Pzzdamp.has_value()) {
              npt_ctx->p_period[2] = *Pzzdamp;
            } else {
              lerr << "Barostat damping parameter Pzzdamp must be > 0" << std::endl;
              std::abort();
            }
          npt_ctx->p_flag[2] = 1;
          npt_ctx->deviatoric_flag = 1;
          }
          if(Pyzstart.has_value()) {
            npt_ctx->p_start[3] = (*Pyzstart)*pascal_to_bar;
            if(Pyzend.has_value()) {          
              npt_ctx->p_stop[3] = (*Pyzend)*pascal_to_bar;
            } else {
              npt_ctx->p_stop[3] = npt_ctx->p_start[3];
            }
            if(Pyzdamp.has_value()) {
              npt_ctx->p_period[3] = *Pyzdamp;
            } else {
              lerr << "Barostat damping parameter Pyzdamp must be > 0" << std::endl;
              std::abort();
            }
          npt_ctx->p_flag[3] = 1;
          npt_ctx->deviatoric_flag = 1;
          npt_ctx->scaleyz = 0;
          }
          if(Pxzstart.has_value()) {
            npt_ctx->p_start[4] = (*Pxzstart)*pascal_to_bar;
            if(Pxzend.has_value()) {          
              npt_ctx->p_stop[4] = (*Pxzend)*pascal_to_bar;
            } else {
              npt_ctx->p_stop[4] = npt_ctx->p_start[4];
            }
            if(Pxzdamp.has_value()) {
              npt_ctx->p_period[4] = *Pxzdamp;
            } else {
              lerr << "Barostat damping parameter Pxzdamp must be > 0" << std::endl;
              std::abort();
            }
          npt_ctx->p_flag[4] = 1;
          npt_ctx->deviatoric_flag = 1;
          npt_ctx->scalexz = 0;          
          }
          if(Pxystart.has_value()) {
            npt_ctx->p_start[5] = (*Pxystart)*pascal_to_bar;
            if(Pxyend.has_value()) {          
              npt_ctx->p_stop[5] = (*Pxyend)*pascal_to_bar;
            } else {
              npt_ctx->p_stop[5] = npt_ctx->p_start[5];
            }
            if(Pxydamp.has_value()) {
              npt_ctx->p_period[5] = *Pxydamp;
            } else {
              lerr << "Barostat damping parameter Pxydamp must be > 0" << std::endl;
              std::abort();
            }
          npt_ctx->p_flag[5] = 1;
          npt_ctx->deviatoric_flag = 1;
          npt_ctx->scalexy = 0;
          }
          if(Pcouple.has_value()) {
            npt_ctx->pcouple = *Pcouple;
            lout << "COUPLE = " << npt_ctx->pcouple << std::endl;
          }
          if(Scalexy.has_value()) npt_ctx->scalexy = *Scalexy;
          if(Scalexz.has_value()) npt_ctx->scalexz = *Scalexz;
          if(Scaleyz.has_value()) npt_ctx->scaleyz = *Scaleyz;
        }

        if (npt_ctx->p_flag[3] || npt_ctx->p_flag[4] || npt_ctx->p_flag[5]) {
          npt_ctx->pstyle = "TRICLINIC";
        } else if (npt_ctx->pcouple == "XYZ") {
          npt_ctx->pstyle = "ISO";
        } else {
          npt_ctx->pstyle = "ANISO";
        }
        
        if (npt_ctx->p_flag[0]) npt_ctx->p_freq[0] = 1./npt_ctx->p_period[0];
        if (npt_ctx->p_flag[1]) npt_ctx->p_freq[1] = 1./npt_ctx->p_period[1];
        if (npt_ctx->p_flag[2]) npt_ctx->p_freq[2] = 1./npt_ctx->p_period[2];
        if (npt_ctx->p_flag[3]) npt_ctx->p_freq[3] = 1./npt_ctx->p_period[3];
        if (npt_ctx->p_flag[4]) npt_ctx->p_freq[4] = 1./npt_ctx->p_period[4];
        if (npt_ctx->p_flag[5]) npt_ctx->p_freq[5] = 1./npt_ctx->p_period[5];

        npt_ctx->omega[0] = npt_ctx->omega[1] = npt_ctx->omega[2] = 0.0;
        npt_ctx->omega_dot[0] = npt_ctx->omega_dot[1] = npt_ctx->omega_dot[2] = 0.0;
        npt_ctx->omega_mass[0] = npt_ctx->omega_mass[1] = npt_ctx->omega_mass[2] = 0.0;
        npt_ctx->omega[3] = npt_ctx->omega[4] = npt_ctx->omega[5] = 0.0;
        npt_ctx->omega_dot[3] = npt_ctx->omega_dot[4] = npt_ctx->omega_dot[5] = 0.0;
        npt_ctx->omega_mass[3] = npt_ctx->omega_mass[4] = npt_ctx->omega_mass[5] = 0.0;
        if (npt_ctx->pstyle == "ISO") size_vector += 2*2*1;
        else if (npt_ctx->pstyle == "ANISO") size_vector += 2*2*3;
        else if (npt_ctx->pstyle == "TRICLINIC") size_vector += 2*2*6;

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

        npt_ctx->p_freq_max = std::max(npt_ctx->p_freq[0],npt_ctx->p_freq[1]);
        npt_ctx->p_freq_max = std::max(npt_ctx->p_freq_max,npt_ctx->p_freq[2]);
        if (npt_ctx->pstyle == "TRICLINIC") {
          npt_ctx->p_freq_max = std::max(npt_ctx->p_freq_max,npt_ctx->p_freq[3]);
          npt_ctx->p_freq_max = std::max(npt_ctx->p_freq_max,npt_ctx->p_freq[4]);
          npt_ctx->p_freq_max = std::max(npt_ctx->p_freq_max,npt_ctx->p_freq[5]);
        }
        npt_ctx->pdrag_factor = 1.0 - ( (*dt) * npt_ctx->p_freq_max * npt_ctx->drag );

        npt_ctx->pdim = npt_ctx->p_flag[0] + npt_ctx->p_flag[1] + npt_ctx->p_flag[2];
        if (npt_ctx->vol0 == 0.0) {
          //	  if (dimension == 3) vol0 = domain->xprd * domain->yprd * domain->zprd;
          //	  else vol0 = domain->xprd * domain->yprd;
          //          npt_ctx->vol0 = *volume;
          npt_ctx->vol0 = determinant(H0);
          // Might change h0_inv to a Mat3D later
          npt_ctx->h0_inv[0] = H0_inv.m11;
          npt_ctx->h0_inv[1] = H0_inv.m22;
          npt_ctx->h0_inv[2] = H0_inv.m33;
          npt_ctx->h0_inv[3] = H0_inv.m23;
          npt_ctx->h0_inv[4] = H0_inv.m13;
          npt_ctx->h0_inv[5] = H0_inv.m12;
        }
      }      
      
    }
  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "init_npt", make_compatible_operator< InitNPTNode > );
  }

}
