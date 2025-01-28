#pragma once

#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/core/domain.h>
#include <exanb/core/print_utils.h>

namespace exaStamp
{

  struct ParrinelloRahmanConfig
  {
    using Mat3d = exanb::Mat3d;

    double m_Text;
    double m_masseNVT;
    double m_Pext;
    double m_masseB;    
    Mat3d m_hmask { 1.,1.,1., 1.,1.,1., 1.,1.,1. };
    Mat3d m_hblend { 1.,0.,0., 0.,1.,0., 0.,0.,1. };
  };

  struct ParrinelloRahmanContext
  {
    using Mat3d = exanb::Mat3d;

    ParrinelloRahmanConfig m_config;

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
      out << "Text        = " << this->m_config.m_Text << std::endl;
      out << "masseNVT    = " << this->m_config.m_masseNVT << std::endl;
      out << "Pext        = " << this->m_config.m_Pext << std::endl;
      out << "m_masseB    = " << this->m_config.m_masseB << std::endl;
      out << "m_gammaNVT  = " << this->m_gammaNVT << std::endl;
      out << "m_gammaNVTp = " << this->m_gammaNVTp << std::endl;
      out << "m_hmask     = " << this->m_config.m_hmask << std::endl;
      out << "m_hblend    = " << this->m_config.m_hblend << std::endl;
      out << "h           = " << this->h << std::endl;
      out << "hp          = " << this->hp << std::endl;
      out << "hpp         = " << this->hpp << std::endl;
      out << "hi          = " << this->hi << std::endl;
      out << "ht          = " << this->ht << std::endl;
      out << "hpt         = " << this->hpt << std::endl;
      out << "hpthp       = " << this->hpthp << std::endl;
      out << "G           = " << this->G << std::endl;
      out << "Gi          = " << this->Gi << std::endl;
      out << "Giht        = " << this->Giht << std::endl;
      out << "Gp          = " << this->Gp << std::endl;
      out << "GiGp        = " << this->GiGp << std::endl;
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

      hp = comp_multiply( hp, m_config.m_hmask );
      hpp = comp_multiply( hpp , m_config.m_hmask );
      if( ! is_identity(m_config.m_hblend) )
      {
        if( is_diagonal(hp) && is_diagonal(hpp) )
        {
          hp = diag_matrix( m_config.m_hblend * Vec3d{ hp.m11 , hp.m22 , hp.m33 } );
          hpp = diag_matrix( m_config.m_hblend * Vec3d{ hpp.m11 , hpp.m22 , hpp.m33 } );
        }
        else
        {
          lerr<<"Cannot apply Hp/Hpp diagonal blending matrix because Hp or Hpp is not diagonal"<<std::endl;
          lerr<<"Hp = "<< hp<<std::endl;
          lerr<<"Hpp = "<< hpp<<std::endl;
          lerr<<"Hblend = "<< m_config.m_hblend<<std::endl;
          std::abort();
        }
      }
    }

  };

}

