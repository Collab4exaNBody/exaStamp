#pragma once

#include <onika/memory/allocator.h>


namespace exaStamp
{

  using onika::memory::DEFAULT_ALIGNMENT;

  // Reaction Field Compute functor
  struct alignas(DEFAULT_ALIGNMENT) RangePotentialForceOp 
  {
    // potential parameters
    double r0min = 0.0;
    double r0max = 0.0;
    double r1min = 0.0;
    double r1max = 0.0;
    double r2min = 0.0;
    double r2max = 0.0;

    const double m_params;

    template<class ComputeBufferT, class CellParticlesT>
    inline void operator ()
      (        
      size_t n,
      const ComputeBufferT& tab,
      double& ep,
      double& fx,
      double& fy,
      double& fz,
      double charge,    // per particle (read only)
      CellParticlesT* unused
      ) const
    {
      FakeMat3d virial;
      this->operator() ( n,tab,ep, fx,fy,fz, charge, virial, unused );
    }

    template<class ComputeBufferT,class Mat3dT, class CellParticlesT>
    inline void operator ()
      (        
      size_t n,
      const ComputeBufferT& tab,
      double& ep,
      double& fx,
      double& fy,
      double& fz,
      double charge,    // per particle (read only)
      Mat3dT& virial,
      CellParticlesT*
      ) const
    {

      // energy and force contributions to the particle
      double _ep = 0.;
      double _fx = 0.;
      double _fy = 0.;
      double _fz = 0.;

      Mat3dT _vir; // default constructor defines all elements to 0
      // assert( _vir.m11==0 && _vir.m12==0 && _vir.m13==0 && _vir.m21==0 && _vir.m22==0 && _vir.m23==0 && _vir.m31==0 && _vir.m32==0 && _vir.m33==0);

      
      // si on parcours les atomes voisins par type
      const unsigned int n_type0 = tab.ext.number_of_nbh0();
      for(unsigned int i=0;i<n_type0;i++)
      {
        if( ! tab.ext.nbh_is_type0(i) ) { _Pragma("omp critical(dbg_mesg)") std::cerr<<"neighbor "<<i<<" is not type0 as expected"<<std::endl; }
        const int j = tab.ext.nbh0(i);
        const double r = sqrt( tab.d2[j] );
        if( ! ( r >= r0min && r <r0max ) ) { _Pragma("omp critical(dbg_mesg)") std::cerr<<"neighbor "<<i<<" is not in range ["<<r0min <<";"<<r0max <<"[ as expected"<<std::endl; }
      }
      const unsigned int n_type1 = tab.ext.number_of_nbh1();
      for(unsigned int i=0;i<n_type1;i++)
      {
        if( ! tab.ext.nbh_is_type1(i) ) { _Pragma("omp critical(dbg_mesg)") std::cerr<<"neighbor "<<i<<" is not type1 as expected"<<std::endl; }
        const int j = tab.ext.nbh1(i);
        const double r = sqrt( tab.d2[j] );
        if( ! ( r >= r1min && r <r1max ) ) { _Pragma("omp critical(dbg_mesg)") std::cerr<<"neighbor "<<i<<" is not in range ["<<r1min <<";"<<r1max <<"[ as expected"<<std::endl; }
      }
      const unsigned int n_type2 = tab.ext.number_of_nbh2();
      for(unsigned int i=0;i<n_type2;i++)
      {
        if( ! tab.ext.nbh_is_type2(i) ) { _Pragma("omp critical(dbg_mesg)") std::cerr<<"neighbor "<<i<<" is not type2 as expected"<<std::endl; }
        const int j = tab.ext.nbh2(i);
        const double r = sqrt( tab.d2[j] );
        if( ! ( r >= r2min && r <r2max ) ) { _Pragma("omp critical(dbg_mesg)") std::cerr<<"neighbor "<<i<<" is not in range ["<<r2min <<";"<<r2max <<"[ as expected"<<std::endl; }
      }

      // si on parcours tous les atomes voisins et qu'on change le calcul en fonction du type
      for(size_t i=0;i<n;i++)
      {
        const double r = std::sqrt(tab.d2[i]);
        const bool is_type0 = tab.ext.nbh_is_type0(i);
        const bool is_type1 = tab.ext.nbh_is_type1(i);
        const bool is_type2 = tab.ext.nbh_is_type2(i);

        double e=0.0, de=0.0;
        // calcul bindon de e,de faisant intervenir is_type0, is_type1, is_type2 mais aussi r et un paramètre bidon m_params
        e = cos( sin(r+is_type0)+is_type1 ) + is_type2 ;
        de = sin(e);
        
        // partie standard du calcul
        de /= r;
        const double drx = tab.drx[i];
        const double dry = tab.dry[i];
        const double drz = tab.drz[i];
        const double fe_x = de * drx;
        const double fe_y = de * dry;
        const double fe_z = de * drz;
        _fx += fe_x;
        _fy += fe_y;
        _fz += fe_z;
        _ep += .5 * e;
        _vir += tensor( Vec3d{fe_x,fe_y,fe_z}, Vec3d{drx,dry,drz} ) * -0.5;
      }
      ep += _ep;
      fx += _fx; 
      fy += _fy; 
      fz += _fz; 
      virial += _vir;
    }
  };

}


