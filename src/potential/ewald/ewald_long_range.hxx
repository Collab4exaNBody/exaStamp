#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <exaStamp/potential/ewald/ewald.h>
#include <exaStamp/particle_species/particle_specie.h>

#include <onika/cuda/cuda.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exanb/core/xform.h>

#include <mpi.h>

#ifdef EXASTAMP_EWALD_PER_PARTICLE_CHARGE
#define EwaldName1 ewald_long_range_pc
#define WORKING_FIELD_CHARGE field::_charge
#else
#define EwaldName1 ewald_long_range
#define WORKING_FIELD_CHARGE field::_type
#endif

#define WORKING_FIELDS field::_fx ,field::_fy ,field::_fz, WORKING_FIELD_CHARGE

#define EwaldName2 EwaldName1

#define EwaldPotentialOperatorName EwaldName2
#define EwaldPotentialStr USTAMP_STR(EwaldName2)

#define EWALD_REGISTER_INIT() _EWALD_REGISTER_INIT( EWALD_CONSTRUCTOR_FUNC_NAME )
#define EWALD_CONSTRUCTOR_FUNC_NAME USTAMP_CONCAT(EwaldPotentialOperatorName,_init)
#define _EWALD_REGISTER_INIT(name) CONSTRUCTOR_ATTRIB void MAKE_UNIQUE_NAME(name,_,__LINE__,ONIKA_CURRENT_PACKAGE_NAME) ()

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  template<class XFormT>
  struct EwaldLongRangeRhoComputeFunc
  {
    const XFormT xform;
    const ParticleSpecie* __restrict__ species = nullptr;
    ReadOnlyEwaldParms p;
    Complexd* __restrict__ ewald_rho = nullptr;
    
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double rx, double ry, double rz, double q ) const
    {
      const Vec3d r = xform.transformCoord( Vec3d{rx,ry,rz} );
      const unsigned int nk = p.nknz;
      for(unsigned int k=0;k<nk;k++)
      {
        const auto gdata = p.Gdata[k];
        const double ps = r.x * gdata.Gx + r.y * gdata.Gy + r.z * gdata.Gz;
        ONIKA_CU_BLOCK_ATOMIC_ADD( ewald_rho[k].r , q * cos(ps) );
        ONIKA_CU_BLOCK_ATOMIC_ADD( ewald_rho[k].i , q * sin(ps) );
      }
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double rx, double ry, double rz, uint8_t type ) const
    {
      this->operator() (rx,ry,rz, species[type].m_charge);
    }
  };

  template<class XFormT>
  struct EwaldLongRangeForceComputeFunc
  {
    const XFormT xform;
    const ParticleSpecie* __restrict__ species = nullptr;
    ReadOnlyEwaldParms p;
    const Complexd* __restrict__ ewald_rho = nullptr;
    
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double & fx, double & fy, double & fz, double rx, double ry, double rz, double q_in ) const
    {
      const double q = 2. * q_in;
      const Vec3d r = xform.transformCoord( Vec3d{rx,ry,rz} );
      const unsigned int nk = p.nknz;

      double lfx = 0.0;
      double lfy = 0.0;
      double lfz = 0.0;
      for(unsigned int k=0;k<nk;k++)
      {
        const auto gdata = p.Gdata[k];
        const double ps = r.x * gdata.Gx + r.y * gdata.Gy + r.z * gdata.Gz;
        const double al = q * gdata.Gc * ( ewald_rho[k].r * sin(ps) - ewald_rho[k].i * cos(ps) );
        lfx += al * gdata.Gx;
        lfy += al * gdata.Gy;
        lfz += al * gdata.Gz;
      }
      fx += lfx;
      fy += lfy;
      fz += lfz;
    }
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double & fx, double & fy, double & fz, double rx, double ry, double rz, uint8_t type ) const
    {
      this->operator() ( fx,fy,fz, rx,ry,rz, species[type].m_charge);
    }
  };

}

namespace exanb
{
  template<class XFormT> struct ComputeCellParticlesTraits< exaStamp::EwaldLongRangeRhoComputeFunc<XFormT> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };

  template<class XFormT> struct ComputeCellParticlesTraits< exaStamp::EwaldLongRangeForceComputeFunc<XFormT> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, WORKING_FIELDS >
    >
  class EwaldPotentialOperatorName : public OperatorNode
  {
//  private:
    // ========= I/O slots =======================
    ADD_SLOT( MPI_Comm   , mpi          , INPUT );
    ADD_SLOT( EwaldParms , ewald_config , INPUT , OPTIONAL ); // optional is here to allow compute_force block to be called before anything is setup
    ADD_SLOT( GridT      , grid         , INPUT_OUTPUT );
    ADD_SLOT( Domain     , domain       , INPUT , REQUIRED );
    ADD_SLOT( double     , rcut_max     , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( EwaldRho   , ewald_rho    , INPUT_OUTPUT , EwaldRho{} );
    

#   ifdef EXASTAMP_EWALD_PER_PARTICLE_CHARGE
    using ewald_rho_field_set_t = FieldSet<field::_rx,field::_ry,field::_rz,field::_charge>;
    using ewald_force_field_set_t = FieldSet<field::_fx,field::_fy,field::_fz,field::_rx,field::_ry,field::_rz,field::_charge>;

    inline double get_particle_charge(const double* charges, const uint8_t*, size_t j) { return charges[j]; }
    static inline constexpr ParticleSpecie* get_species() { return nullptr; }
#   else
    ADD_SLOT( ParticleSpecies       , species           , INPUT );

    using ewald_rho_field_set_t = FieldSet<field::_rx,field::_ry,field::_rz,field::_type>;
    using ewald_force_field_set_t = FieldSet<field::_fx,field::_fy,field::_fz,field::_rx,field::_ry,field::_rz,field::_type>;

    inline ParticleSpecie* get_species() { return species->data(); }
    inline double get_particle_charge(const double*, const uint8_t* types, size_t j) { return (*species)[ types[j] ].m_charge; }
#   endif

    static constexpr ewald_rho_field_set_t ewald_rho_field_set = {};
    static constexpr ewald_force_field_set_t ewald_force_field_set = {};

  public:

    // Operator execution
    inline void execute () override final
    {
      ldbg<<"------------------------------"<<std::endl<<std::flush;
      ldbg<<"Beginning of long range energy"<<std::endl<<std::flush;
      ldbg<<"------------------------------"<<std::endl<<std::flush;

      if( ! ewald_config.has_value() )
      {
        ldbg << "ewald_config not set, skip ewal_long_range computation" << std::endl;
        return ;
      }

      *rcut_max = std::max( *rcut_max , ewald_config->radius );
      
      if( grid->number_of_cells() == 0 ) return;

      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();

      // Recupération du nombre de point k
      const size_t nk = ewald_config->nknz;
      ewald_rho->nk = nk;
      ewald_rho->rho.resize(nk);
      
      // Recupération des Gx=2.pi*kx/Lx, ...
      const EwaldCoeffs * __restrict__ Gdata = ewald_config->Gdata.data();
      
      // Rho computation
      // Rho=qi * exp(iG.r)
      int max_nt = omp_get_max_threads();

      ldbg<<"      number of k points : "<<nk<<std::endl<<std::flush;
      ldbg<<"      number of threads : "<<max_nt<<std::endl<<std::flush;

      // size_t total_particles = 0;
      Mat3d xform = domain->xform();

      bool allow_cuda_exec = ( global_cuda_ctx() != nullptr );
      if( allow_cuda_exec ) allow_cuda_exec = global_cuda_ctx()->has_devices();
      if( allow_cuda_exec )
      {
        // const int streamIndex = 0;
	      //lout << "Ewald rho : Cuda version" << std::endl;
        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MEMSET( ewald_rho->rho.data(), 0, sizeof(Complexd)*nk, global_cuda_ctx()->getThreadStream(0) ) );
        /* if( domain->xform_is_identity() )
        {
          EwaldLongRangeRhoComputeFunc<NullXForm> func = { {} , get_species() , *ewald_config , ewald_rho->rho.data() };
          compute_cell_particles( *grid , false , func , ewald_rho_field_set , parallel_execution_context() );
        }
        else */
        {
          EwaldLongRangeRhoComputeFunc<LinearXForm> func = { {xform} , get_species() , *ewald_config , ewald_rho->rho.data() };
          compute_cell_particles( *grid , false , func , ewald_rho_field_set , parallel_execution_context() );
        }
      }
      else
      {
	      //lout << "Ewald rho : CPU version" << std::endl;
        onika::memory::CudaMMVector<Complexd> tmp_rho( nk * max_nt );

#       pragma omp parallel
        {
          const int tid = omp_get_thread_num();
          const int nt = omp_get_num_threads();
          assert( nt <= max_nt );
          assert( tid < max_nt );

          Complexd * __restrict__ local_rho = tmp_rho.data() + ( tid * nk );
          for(size_t k=0; k<nk; ++k) { local_rho[k] = Complexd{0.,0.}; }
          
          GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) nowait /*reduction(+:total_particles)*/ )
          {
            size_t i = grid_ijk_to_index( dims , loc + gl );
            size_t n = cells[i].size();
            // total_particles += n;
            
            const double* __restrict__ charges = cells[i].field_pointer_or_null(field::charge); ONIKA_ASSUME_ALIGNED(charges);
            const uint8_t* __restrict__ types = cells[i].field_pointer_or_null(field::type); ONIKA_ASSUME_ALIGNED(types);
            const double* __restrict__ rx = cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
            const double* __restrict__ ry = cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
            const double* __restrict__ rz = cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
            
            for(size_t j=0;j<n;j++)
            {
              const double q = get_particle_charge(charges,types,j);
              Vec3d r = { rx[j] , ry[j] , rz[j] };
              r = xform * r;

              for (size_t k=0; k<nk; ++k)
              {
                const auto gdata = Gdata[k];
                const double ps = r.x * gdata.Gx + r.y * gdata.Gy + r.z * gdata.Gz;
                local_rho[k] += q * Complexd{ cos(ps) , sin(ps) };
              }
            }
          }
          GRID_OMP_FOR_END

          size_t start = ( nk * tid ) / nt;
          size_t end = ( nk * (tid+1) ) / nt;

#         pragma omp barrier

          for(size_t i=start;i<end;i++)
          {
            ewald_rho->rho[i] = Complexd{ 0.0 , 0.0 };
          }
          for(int t=0;t<nt;t++)
          {
            local_rho = tmp_rho.data() + ( t * nk );
            for(size_t i=start;i<end;i++)
            {
              ewald_rho->rho[i] += local_rho[i]; // tmp_rho[i] = rho[i] totale du domaine =\sum_i  q_i * exp(i G.r)
            }
          }
        }
      
      }

      static_assert( sizeof(Complexd) == 2*sizeof(double) );
      MPI_Allreduce(MPI_IN_PLACE, (double*) ewald_rho->rho.data(),nk*2,MPI_DOUBLE,MPI_SUM,*mpi);

  // forces and energy
      if( allow_cuda_exec )
      {
	//lout << "Ewald force : Cuda version" << std::endl;
        /* if( domain->xform_is_identity() )
        {
          EwaldLongRangeForceComputeFunc<NullXForm> func = { {} , get_species() , *ewald_config , ewald_rho->rho.data() };
          compute_cell_particles( *grid , false , func , ewald_force_field_set , parallel_execution_context() );
        }
        else */
        {
          EwaldLongRangeForceComputeFunc<LinearXForm> func = { {xform} , get_species() , *ewald_config , ewald_rho->rho.data() };
          compute_cell_particles( *grid , false , func , ewald_force_field_set , parallel_execution_context() );
        }
      }
      else
      {
	      //lout << "Ewald force : CPU version" << std::endl;
#       pragma omp parallel
        {
          GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) nowait )
          {
            size_t i = grid_ijk_to_index( dims , loc + gl );
            size_t n = cells[i].size();

            const double* __restrict__ charges = cells[i].field_pointer_or_null(field::charge); ONIKA_ASSUME_ALIGNED(charges);
            const uint8_t* __restrict__ types = cells[i].field_pointer_or_null(field::type); ONIKA_ASSUME_ALIGNED(types);
            const double* __restrict__ rx = cells[i][field::rx]; ONIKA_ASSUME_ALIGNED(rx);
            const double* __restrict__ ry = cells[i][field::ry]; ONIKA_ASSUME_ALIGNED(ry);
            const double* __restrict__ rz = cells[i][field::rz]; ONIKA_ASSUME_ALIGNED(rz);
            //double* __restrict__ ep = cells[i][field::ep]; ONIKA_ASSUME_ALIGNED(ep);
            double* __restrict__ fx = cells[i][field::fx]; ONIKA_ASSUME_ALIGNED(fx);
            double* __restrict__ fy = cells[i][field::fy]; ONIKA_ASSUME_ALIGNED(fy);
            double* __restrict__ fz = cells[i][field::fz]; ONIKA_ASSUME_ALIGNED(fz);
   
            for(size_t j=0;j<n;j++)
            {
              double q = 2. * get_particle_charge(charges,types,j);
              Vec3d r = { rx[j] , ry[j] , rz[j] };
              r = xform * r;
              //ep[j] += re;

              double lfx = 0.0;
              double lfy = 0.0;
              double lfz = 0.0;
              for (size_t k=0; k<nk; ++k)
              {
                const auto gdata = Gdata[k];
                const double ps = r.x * gdata.Gx + r.y * gdata.Gy + r.z * gdata.Gz;
                const double al = q * gdata.Gc * ( ewald_rho->rho[k].r * sin(ps) - ewald_rho->rho[k].i * cos(ps) );
                lfx += al * gdata.Gx;
                lfy += al * gdata.Gy;
                lfz += al * gdata.Gz;
              }

              fx[j] += lfx;
              fy[j] += lfy;
              fz[j] += lfz;
            }
          }
          GRID_OMP_FOR_END
        }
      }
      
      ldbg<<"------------------------"<<std::endl<<std::flush;
      ldbg<<"End of long range energy"<<std::endl<<std::flush;
      ldbg<<"------------------------"<<std::endl<<std::flush;
    }

  };

  namespace tplhelper {
  template<class GridT> using EwaldPotentialOperatorName = ::exaStamp::EwaldPotentialOperatorName<GridT>;
  }

  // === register factories ===  
  EWALD_REGISTER_INIT()
  {  
    OperatorNodeFactory::instance()->register_factory( EwaldPotentialStr , make_grid_variant_operator< tplhelper::EwaldPotentialOperatorName > );
  }

}


