#pragma once

#include <yaml-cpp/yaml.h>
#include <onika/physics/units.h>
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_yaml.h>
#include <cmath>

#include <onika/cuda/cuda.h>
#include <onika/cuda/ro_shallow_copy.h>

namespace exaStamp
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  namespace ewald_constants
  {
    static constexpr double epsilonZero = 0x1.337f14f782bb5p-21; //0x1.2c4a1a79b5aafp-11 * 0.001;
    static constexpr double fpe0 = 0x1.e303a50b0d1ap-18; // 0x1.d7b18f2ccacb6p-8 * 0.001 ;
  }

  struct EwaldRho
  {
    size_t nk = 0;
    onika::memory::CudaMMVector<Complexd> rho;
  };

  struct EwaldCoeffs
  {
    double Gx;
    double Gy;
    double Gz;
    double Gc;
  };
  
  struct alignas(DEFAULT_ALIGNMENT) EwaldParms
  {
    double alpha = 0.0;
    double radius = 0.0; 
    double epsilon = 0.0;
    ssize_t kmax = 0;
    ssize_t nk = 0;
    ssize_t nknz = 0;

    double gm = 0.0;
    double gm_sr = 0.0; // the gm used in short range computation
    double bt = 0.0;
    double bt_sr = 0.0; // the bt value used in short range computation
    
    double GnMax = 0.0;
    double volume = 0.0;
    Vec3d lm = { 0.0 , 0.0 , 0.0 };
    
    onika::memory::CudaMMVector<EwaldCoeffs> Gdata;
//    onika::memory::CudaMMVector<double> Gy;
//    onika::memory::CudaMMVector<double> Gz;
//    onika::memory::CudaMMVector<double> Gc;
  };

  struct ReadOnlyEwaldParms
  {
    double alpha = 0.0;
//    double radius = 0.0; 
//    double epsilon = 0.0;
//    ssize_t kmax = 0;
//    ssize_t nk = 0;
    ssize_t nknz = 0;

//    double gm = 0.0;
    double gm_sr = 0.0; // the gm used in short range computation
//    double bt = 0.0;
    double bt_sr = 0.0; // the bt value used in short range computation
    
    double GnMax = 0.0;
//    double volume = 0.0;
//    Vec3d lm = { 0.0 , 0.0 , 0.0 };
    
    const EwaldCoeffs* __restrict__ Gdata = nullptr;
    
    ReadOnlyEwaldParms() = default;
    ReadOnlyEwaldParms(const ReadOnlyEwaldParms&) = default;
    ReadOnlyEwaldParms(ReadOnlyEwaldParms&&) = default;
    ReadOnlyEwaldParms& operator = (const ReadOnlyEwaldParms&) = default;
    ReadOnlyEwaldParms& operator = (ReadOnlyEwaldParms&&) = default;
    
    inline ReadOnlyEwaldParms( const EwaldParms & p )
      : alpha( p.alpha )
//      , radius( p.radius )
//      , epsilon( p.epsilon )
//      , kmax( p.kmax )
//      , nk( p.nk )
      , nknz( p.nknz )

//      , gm( p.gm )
      , gm_sr( p.gm_sr )
//      , bt( p.bt )
      , bt_sr( p.bt_sr )
      
      , GnMax( p.GnMax )
//      , volume( p.volume )
//      , lm( p.lm )
      
      , Gdata( p.Gdata.data() )
//      , Gy( p.Gy.data() )
//      , Gz( p.Gz.data() )
//      , Gc( p.Gc.data() )
    {}
  };


  // core computation kernel for reaction field potential
  template<class EwaldParmsT>
  ONIKA_HOST_DEVICE_FUNC static inline void ewald_compute_energy(const EwaldParmsT& p, double c1, double c2, double r, double& e, double& de)
  {
    const double ar = p.alpha * r;
    const double cf = p.gm_sr * c1 * c2; // (1/(4.pi.epsilon0)qi.qj

    e = cf * erfc(ar) / r;
    de = - (cf * p.bt_sr * exp(-ar * ar) + e) / r;
  }

  inline void ewald_init_parameters(double alpha, double radius, double epsilon, unsigned int in_kmax, const Vec3d& domainSize, EwaldParms& p , ::exanb::LogStreamWrapper& ldbg )
  {
    // using ewald_constants::epsilonZero;
    using ewald_constants::fpe0;
    using ewald_constants::epsilonZero;

    p.alpha = alpha;
    p.radius = radius;
    p.epsilon = epsilon;
    p.kmax = in_kmax ;

    // ewald_init_parameters should be called again if volume changes
    p.volume = domainSize.x * domainSize.y * domainSize.z ;
    p.lm = Vec3d{ 2.*M_PI/domainSize.x , 2.*M_PI/domainSize.y , 2.*M_PI/domainSize.z };

    if( p.volume == 0.0 ) return;

    // compute alpha
    if(p.alpha <= 0.)
    {
        ldbg<<"alpha < 0 : calcul de alpha "<<std::endl<<std::flush;
        int test = 0;
        double alpha=0.;  
        double dx=5./p.radius;
        long i=0;
 
        long numit=0;
        long numitmax=10000000;

        int test1 = 0;
        while(test1 == 0 && numit < numitmax){ 
                i=0;
                while(test == 0){
                        i++;
                        alpha = i*dx;
                        double y = (std::erfc(alpha*p.radius));
                        if(y<p.epsilon){
                                ldbg<<"alpha="<<alpha<<" y="<<y<<" i="<<i<<" dx="<<dx<<" numit="<<numit<<std::endl;
                                test=1;
                        }
                }
                if(i > 1000){
                        test1++;
                }else{
                        dx = dx/2;
                        test = 0;
                        numit++;
                }      
        }
        p.alpha = alpha;
    }
    if(p.alpha == 0.)
    {
        lerr<<"######## ERROR in ewald_init_parameters : alpha ="<<p.alpha<<" - Decrease epsilon of Ewald method"<<std::endl<<std::flush;
        std::abort();
    }else{

        ldbg<<" alpha ="<<p.alpha<<std::endl<<std::flush;
    }

    // compute kmax
    if(p.kmax <= 0)
    {
        if(p.alpha == 0.)
        {
          lerr<<"######## ERROR in ewald_init_parameters : alpha ="<<p.alpha<<" in kmax calculation"<<std::endl<<std::flush;
          std::abort();
        }

        ldbg<<"kmax < 0 : calcul de kmax "<<std::endl;
        int test = 0;
        double gcible=0.;  
        double dx=1.e30;
        long i=0;

        long numit=0;
        long numitmax=10000000000;

        int test1 = 0;
        while(test1 == 0 && numit < numitmax){ 
                i=0;
                while(test == 0){
                        i++;
                        gcible = i*dx;
                        double y=std::exp(-gcible/(4.*p.alpha*p.alpha))/gcible;
                        if(y<p.epsilon){
                                ldbg<<"gcible="<<gcible<<" y="<<y<<" i="<<i<<" dx="<<dx<<std::endl;
                                test=1;
                        }
                }
                if(i > 1000){
                        test1++;
                }else{
                        dx = dx/2.;
                        test = 0;
                        numit++;
                }
        }       

        p.kmax = std::floor(std::sqrt(gcible)*domainSize.x/(2.*M_PI))+1;
    }

    if(p.kmax < 2)
    {
        lerr<<"######## ERROR in ewald_init_parameters : kmax ="<<p.kmax<<" - Decrease epsilon of Ewald method"<<std::endl;
        std::abort();
    }else{

        ldbg<<" kmax ="<<p.kmax<<std::endl<<std::flush;
    }
        


    double x = 2. * M_PI * p.kmax / domainSize.x;
    p.GnMax = x * x;
    p.nk = (2 * p.kmax + 1) * (2 * p.kmax + 1) * (2 * p.kmax + 1) - 1;
   
    p.Gdata.resize( p.nk );
//    p.Gy.assign( p.nk, 0. );
//    p.Gz.assign( p.nk, 0. );
//    p.Gc.resize( p.nk );

    p.bt = 2. * M_PI / fpe0 / p.volume;
    p.bt_sr = 2. * p.alpha / std::sqrt(M_PI);
    p.gm = 1. / (4. * p.alpha * p.alpha);
    p.gm_sr = 0.25 / M_PI / epsilonZero;


    ldbg << "Ewald parameters"          <<std::endl<<std::flush;
    ldbg << "   alpha="<<p.alpha        <<std::endl<<std::flush;
    ldbg << "   E short range(rc)="<<(std::erfc(p.alpha*p.radius))/p.radius<<std::endl<<std::flush;
    ldbg << "   radius="<<p.radius      <<std::endl<<std::flush;
    ldbg << "   epsilon="<<p.epsilon    <<std::endl<<std::flush;
    ldbg << "   kmax="<<p.kmax          <<std::endl<<std::flush;
    ldbg << "   nk="<<p.nk              <<std::endl<<std::flush;
    ldbg << "   bt="<<p.bt              <<std::endl<<std::flush;                //p.bt=1/(2.epsilon_0.V)
    ldbg << "   gm="<<p.gm              <<std::endl<<std::flush;
    ldbg << "   lm=" << p.lm.x <<','<<p.lm.y<<','<<p.lm.z  <<std::endl<<std::flush;   // lm.x = 2 Pi / Lx
    ldbg << "   volume="<< p.volume     <<std::endl<<std::flush;
    ldbg << "   fpe0="<< fpe0           <<std::endl<<std::flush;                // fpe0=1/p.gm_sr = (4piEpsilon0)
    ldbg << "   sr="<< p.gm_sr          <<std::endl<<std::flush;                // p.gm_sr = 1/(4piEpsilon0)
    ldbg << "   GnMax="<< p.GnMax       <<std::endl<<std::flush;



    // note : this triple-loop is parallelizable
    //const ssize_t kmin = - static_cast<ssize_t>(p.kmax);
    const ssize_t kmin = - static_cast<ssize_t>(p.kmax);
    const ssize_t kmax = static_cast<ssize_t>(p.kmax);
    size_t kk = 0;
    double gcmin = 1e30;
    double gcmax = 0.0;

    for (ssize_t kx=kmin; kx<=kmax; ++kx )
    {
      for (ssize_t ky=kmin; ky<=kmax; ++ky)
      {
        for (ssize_t kz=kmin; kz<=kmax; ++kz)
        {
	        const ssize_t nk = kx*kx + ky*ky + kz*kz;
	        if( nk > 0 ) // <=>  kx!=0 || ky!=0 || kz!=0
          {
           //if(kk<p.nk)
           //{
            const double Gx_kk = kx * p.lm.x;   // p.Gx[kk] = 2 pi kx / Lx
            const double Gy_kk = ky * p.lm.y;   // p.Gy[kk] = 2 pi ky / Ly
            const double Gz_kk = kz * p.lm.z;   // p.Gz[kk] = 2 pi kz / Lz
            const double Gn_kk = Gx_kk*Gx_kk + Gy_kk*Gy_kk + Gz_kk*Gz_kk; // Gn_kk = G^2
            if (Gn_kk <= p.GnMax)
            {
              assert( kk < static_cast<size_t>(p.nk) );

              double Gc_kk = std::exp(-p.gm * Gn_kk) / Gn_kk; //p.Gc[kk] = exp(-G^2/(4 alpha^2)/G^2

              if(Gc_kk<gcmin) gcmin=Gc_kk;
              if(Gc_kk>gcmax) gcmax=Gc_kk;

              ldbg<<"G^2="<<Gn_kk<<" exp(-G^2/(4 alpha^2)/G^2)="<<Gc_kk<<std::endl<<std::flush;

              Gc_kk *= p.bt; //p.Gc[kk] = 1/(2.epsilon_0.V) exp(-G^2/(4 alpha^2)/G^2

              p.Gdata[kk] = EwaldCoeffs{ Gx_kk , Gy_kk , Gz_kk , Gc_kk };

              ++kk;
            }
            //} // if(kk<p.nk)
          }
        }
      }
    }

    ldbg<<"   exp(-G^2/4alpha^2)/G^2 : minimum value :"<<gcmin<<std::endl<<std::flush;
    ldbg<<"                          : maximum value :"<<gcmax<<std::endl<<std::flush;

    // number of non zero values
    p.nknz = kk;    
    ldbg << "   number of k points="<< p.nknz       <<std::endl<<std::flush;

    // adjust coeffs size
    p.Gdata.resize( p.nknz );
    p.Gdata.shrink_to_fit();
  }

  inline void ewald_init_parameters(double alpha, double radius, double epsilon, unsigned int in_kmax, const Vec3d& domainSize, EwaldParms& p )
  {
    ewald_init_parameters(alpha,radius,epsilon,in_kmax,domainSize,p , ::exanb::ldbg<<"" );
  }

} // end of namespace exaStamp

// specialize ReadOnlyShallowCopyType so ReadOnlyEwaldParms is the read only type for EwaldParms
namespace onika { namespace cuda {
  template<> struct ReadOnlyShallowCopyType< exaStamp::EwaldParms > { using type = exaStamp::ReadOnlyEwaldParms; };
} }

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  using exaStamp::EwaldParms;
  
  using onika::physics::Quantity;
  using exanb::Vec3d;
  using exaStamp::ewald_init_parameters;

  template<> struct convert<EwaldParms>
  {
    static bool decode(const Node& node, EwaldParms& v)
    {
      if( !node.IsMap() ) { return false; }
      double alpha = 0.0;
      if( node["alpha"] )
      {
        alpha = node["alpha"].as<Quantity>().convert();
      }
      Vec3d domSize = node["size"].as<Vec3d>();
      ewald_init_parameters(
        alpha,
        node["radius"].as<Quantity>().convert(),
        node["epsilon"].as<Quantity>().convert(),
        node["kmax"].as<long>(), domSize, v);
      return true;
    }
  };
}

