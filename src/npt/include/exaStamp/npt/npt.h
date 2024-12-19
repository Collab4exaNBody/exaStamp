#pragma once

#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/core/domain.h>
#include <exanb/core/print_utils.h>

namespace exaStamp
{

  struct NPTConfig
  {
    using Mat3d = exanb::Mat3d;
    double m_Tstart;
    double m_Tend;    
    double m_Tdamp;
    Mat3d  m_Pstart { 0.,0.,0., 0.,0.,0., 0.,0.,0. };
    Mat3d  m_Pend { 0.,0.,0., 0.,0.,0., 0.,0.,0. };
    double m_Pdamp;
    std::string m_mode;
  };

  struct NPTContext
  {
    using Mat3d = exanb::Mat3d;
    NPTConfig m_config;

    int dimension, which;
    double dtv, dtf, dthalf, dt4, dt8, dto;
    double boltz, nktv2p, tdof;
    double vol0;    // reference volume
    double t0;      // reference temperature

    // For thermostat purposes
    double t_start, t_stop;
    double t_current, t_target, ke_target;
    double t_freq;
    int tstat_flag;

    // For barostat purposes
    int pstyle, pcouple, allremap;
    int p_flag[6];    // 1 if control P on this dim, 0 if not
    double p_start[6], p_stop[6];
    double p_freq[6], p_target[6];
    double omega[6], omega_dot[6];
    double omega_mass[6];
    double p_current[6];
    double drag, tdrag_factor;     // drag factor on particle thermostat
    double pdrag_factor;           // drag factor on barostat
    int pstat_flag;    // 1 if control P

    double p_temp;    // target temperature for barostat
    int p_temp_flag;

    double *eta, *eta_dot;    // chain thermostat for particles
    double *eta_dotdot;
    double *eta_mass;
    int mtchain;                 // length of chain
    int mtchain_default_flag;    // 1 = mtchain is default

    double *etap;    // chain thermostat for barostat
    double *etap_dot;
    double *etap_dotdot;
    double *etap_mass;
    int mpchain;    // length of chain

    bool mtk_flag;         // 0 if using Hoover barostat
    int pdim;             // number of barostatted dims
    double p_freq_max;    // maximum barostat frequency

    double p_hydro;    // hydrostatic target pressure

    int nc_tchain, nc_pchain;
    double factor_eta;
    double sigma[6];        // scaled target stress
    double fdev[6];         // deviatoric force on barostat
    int deviatoric_flag;    // 0 if target stress tensor is hydrostatic
    double h0_inv[6];       // h_inv of reference (zero strain) box
    int nreset_h0;          // interval for resetting h0

    double mtk_term1, mtk_term2;    // Martyna-Tobias-Klein corrections

    int eta_mass_flag;      // 1 if eta_mass updated, 0 if not.
    int omega_mass_flag;    // 1 if omega_mass updated, 0 if not.
    int etap_mass_flag;     // 1 if etap_mass updated, 0 if not.
    int dipole_flag;        // 1 if dipole is updated, 0 if not.
    int dlm_flag;           // 1 if using the DLM rotational integrator, 0 if not

    int scaleyz;     // 1 if yz scaled with lz
    int scalexz;     // 1 if xz scaled with lz
    int scalexy;     // 1 if xy scaled with ly
    int flipflag;    // 1 if box flips are invoked as needed

    int pre_exchange_flag;    // set if pre_exchange needed for box flips

    double fixedpoint[3];    // location of dilation fixed-point
    
    // evolutive parameters
    double m_gammaNVT = 0.0;
    double m_gammaNVTp = 0.0;
    Mat3d h, hp, hpp; // these matrices are to be initialized to 0,
    Mat3d G, Gp;      // this is done in the default constructor of Mat3d

    // derived quantities, updated in updateMembers
    Mat3d hi, ht, Gi, Giht, GiGp, hpt, hpthp;

    template<class StreamT>
    inline void print(StreamT& out)
    {
      out << exanb::default_stream_format;
      out << "Tstart        = " << this->m_config.m_Tstart << std::endl;
      out << "Tend          = " << this->m_config.m_Tend << std::endl;
      out << "Tdamp         = " << this->m_config.m_Tdamp << std::endl;
      out << "--------------- "  << std::endl;
      out << "Pxx_start     = " << this->m_config.m_Pstart.m11 << std::endl;
      out << "Pyy_start     = " << this->m_config.m_Pstart.m22 << std::endl;
      out << "Pzz_start     = " << this->m_config.m_Pstart.m33 << std::endl;      
      out << "Pxy_start     = " << this->m_config.m_Pstart.m12 << std::endl;
      out << "Pxz_start     = " << this->m_config.m_Pstart.m13 << std::endl;
      out << "Pyz_start     = " << this->m_config.m_Pstart.m23 << std::endl;      
      out << "Pdamp         = " << this->m_config.m_Pdamp << std::endl;
    }

    inline void updateMembers()
    {
      using namespace exanb;

      // h = xform * diag_matrix(domain.extent());
      hi = inverse(h);
      ht = transpose(h);
      hpt = transpose(hp);
      hpthp = hpt * hp;
      
      G = ht * h;
      Gi = inverse(G);
      Giht = Gi * ht;
      Gp = diag_matrix({0.,0.,0.}); // this is in stamp and yes it's strange
      GiGp = Gi * Gp;
    }

    inline void apply_mask()
    {
      using namespace exanb;

      // hp = comp_multiply( hp, m_config.m_hmask );
      // hpp = comp_multiply( hpp , m_config.m_hmask );
      // if( ! is_identity(m_config.m_hblend) )
      // {
      //   if( is_diagonal(hp) && is_diagonal(hpp) )
      //   {
      //     hp = diag_matrix( m_config.m_hblend * Vec3d{ hp.m11 , hp.m22 , hp.m33 } );
      //     hpp = diag_matrix( m_config.m_hblend * Vec3d{ hpp.m11 , hpp.m22 , hpp.m33 } );
      //   }
      //   else
      //   {
      //     lerr<<"Cannot apply Hp/Hpp diagonal blending matrix because Hp or Hpp is not diagonal"<<std::endl;
      //     lerr<<"Hp = "<< hp<<std::endl;
      //     lerr<<"Hpp = "<< hpp<<std::endl;
      //     lerr<<"Hblend = "<< m_config.m_hblend<<std::endl;
      //     std::abort();
      //   }
      // }
    }

  };

}

