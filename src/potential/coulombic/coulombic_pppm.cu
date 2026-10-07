/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particles.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/flat_tuple.h>
#include <onika/cuda/cuda.h>
#include <onika/cuda/cuda_context.h>
#include <onika/parallel/parallel_for.h>
#include <mpi.h>

#include <exaStamp/potential/coulombic/pppm.h>
#include <exaStamp/potential/coulombic/pppm_fft.h>

namespace exaStamp
{
inline namespace coulombic_ewald // distinct symbols from the legacy ewald plugin (plugins are loaded RTLD_GLOBAL)
{
  using namespace exanb;

  template<bool PerAtomCharge, class ChargeOrTypeT>
  ONIKA_HOST_DEVICE_FUNC static inline double pppm_particle_charge( const ParticleSpecie* __restrict__ species, ChargeOrTypeT ct )
  {
    if constexpr ( PerAtomCharge ) return ct;
    else return species[ct].m_charge;
  }

  // LAMMPS PPPM::make_rho : charge density on the (global, replicated) mesh
  template<bool PerAtomCharge>
  struct PPPMSpreadFunc
  {
    ReadOnlyPPPMParameters p;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    double * __restrict__ m_density = nullptr;

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      const double q = pppm_particle_charge<PerAtomCharge>( m_species , ct );
      if( q == 0.0 ) return;
      PPPMStencil st;
      pppm_stencil( p , Vec3d{rx,ry,rz} , st );
      const double z0 = p.delvolinv * q;
      for( int n = 0 ; n < p.order ; n++ )
      {
        int mz = st.iz + n; if( mz >= p.nz ) mz -= p.nz;
        const double y0 = z0 * st.wz[n];
        for( int m = 0 ; m < p.order ; m++ )
        {
          int my = st.iy + m; if( my >= p.ny ) my -= p.ny;
          const double x0 = y0 * st.wy[m];
          double * __restrict__ row = m_density + ( size_t(mz) * p.ny + my ) * p.nx;
          for( int l = 0 ; l < p.order ; l++ )
          {
            int mx = st.ix + l; if( mx >= p.nx ) mx -= p.nx;
            ONIKA_CU_ATOMIC_ADD( row[mx] , x0 * st.wx[l] );
          }
        }
      }
    }
  };

  // real space meshes interpolated back to particles
  struct PPPMMeshes
  {
    const double * __restrict__ vdx = nullptr; // gradient of the potential
    const double * __restrict__ vdy = nullptr;
    const double * __restrict__ vdz = nullptr;
    const double * __restrict__ u = nullptr;   // potential (per particle energy)
    const double * __restrict__ v = nullptr;   // 6 virial meshes, xx yy zz xy xz yz, nfft apart
  };

  // LAMMPS PPPM::fieldforce_ik (+ fieldforce_peratom and per atom self/background correction)
  template<bool PerAtomCharge, bool ComputeEnergy, bool ComputeVirial>
  struct PPPMForceFunc
  {
    ReadOnlyPPPMParameters p;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    PPPMMeshes m;
    size_t nfft = 0;

    ONIKA_HOST_DEVICE_FUNC inline void compute( double q, double rx, double ry, double rz, Vec3d& f, double& ep, Mat3d& vir ) const
    {
      if( q == 0.0 ) return;
      PPPMStencil st;
      pppm_stencil( p , Vec3d{rx,ry,rz} , st );
      double ekx = 0.0, eky = 0.0, ekz = 0.0;
      double u = 0.0;
      double v0 = 0.0, v1 = 0.0, v2 = 0.0, v3 = 0.0, v4 = 0.0, v5 = 0.0;
      for( int n = 0 ; n < p.order ; n++ )
      {
        int mz = st.iz + n; if( mz >= p.nz ) mz -= p.nz;
        const double z0 = st.wz[n];
        for( int mm = 0 ; mm < p.order ; mm++ )
        {
          int my = st.iy + mm; if( my >= p.ny ) my -= p.ny;
          const double y0 = z0 * st.wy[mm];
          const size_t row = ( size_t(mz) * p.ny + my ) * p.nx;
          for( int l = 0 ; l < p.order ; l++ )
          {
            int mx = st.ix + l; if( mx >= p.nx ) mx -= p.nx;
            const size_t idx = row + mx;
            const double x0 = y0 * st.wx[l];
            ekx -= x0 * m.vdx[idx];
            eky -= x0 * m.vdy[idx];
            ekz -= x0 * m.vdz[idx];
            if constexpr ( ComputeEnergy ) u += x0 * m.u[idx];
            if constexpr ( ComputeVirial )
            {
              v0 += x0 * m.v[idx];
              v1 += x0 * m.v[nfft+idx];
              v2 += x0 * m.v[2*nfft+idx];
              v3 += x0 * m.v[3*nfft+idx];
              v4 += x0 * m.v[4*nfft+idx];
              v5 += x0 * m.v[5*nfft+idx];
            }
          }
        }
      }
      const double qfactor = COULOMB_CONSTANT * q;
      f.x += qfactor * ekx;
      f.y += qfactor * eky;
      f.z += qfactor * ekz;
      if constexpr ( ComputeEnergy )
      {
        ep += 0.5 * qfactor * u + pppm_self_energy( p , q );
        if constexpr ( ComputeVirial )
        {
          const double h = 0.5 * qfactor;
          vir.m11 += h * v0; vir.m22 += h * v1; vir.m33 += h * v2;
          vir.m12 += h * v3; vir.m21 += h * v3;
          vir.m13 += h * v4; vir.m31 += h * v4;
          vir.m23 += h * v5; vir.m32 += h * v5;
        }
      }
    }

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double & fx, double & fy, double & fz, double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      static_assert( ! ComputeEnergy && ! ComputeVirial );
      Vec3d f = {0.,0.,0.}; double ep = 0.0; Mat3d vir;
      compute( pppm_particle_charge<PerAtomCharge>( m_species , ct ) , rx, ry, rz, f, ep, vir );
      fx += f.x; fy += f.y; fz += f.z;
    }

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( double & fx, double & fy, double & fz, double & ep, Mat3d & virial, double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      static_assert( ComputeEnergy && ComputeVirial );
      Vec3d f = {0.,0.,0.}; double e = 0.0; Mat3d vir = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      compute( pppm_particle_charge<PerAtomCharge>( m_species , ct ) , rx, ry, rz, f, e, vir );
      fx += f.x; fy += f.y; fz += f.z; ep += e; virial += vir;
    }
  };

  // ------------- per mesh point operations (parallel_for, one thread per mesh point on GPU) -------------

  struct PPPMZeroFunc
  {
    double * __restrict__ a = nullptr;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const { a[i] = 0.0; }
  };

  // w = density (real) , then w *= scale * G(k) once transformed
  struct PPPMLoadDensityFunc
  {
    const double * __restrict__ density = nullptr;
    Complexd * __restrict__ w = nullptr;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const { w[i] = Complexd{ density[i] , 0.0 }; }
  };

  struct PPPMApplyGreenFunc
  {
    Complexd * __restrict__ w = nullptr;
    const double * __restrict__ greensfn = nullptr;
    double scaleinv = 1.0;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const
    {
      const double s = scaleinv * greensfn[i];
      w[i].r *= s;
      w[i].i *= s;
    }
  };

  // w2 = V(k) times a per mesh point factor, before a backward FFT
  enum class PPPMFactor { COPY , VIRIAL , IK };
  struct PPPMFactorFunc
  {
    const Complexd * __restrict__ w1 = nullptr;
    Complexd * __restrict__ w2 = nullptr;
    PPPMFactor mode = PPPMFactor::COPY;
    const double * __restrict__ coef = nullptr; // vg (6 per point) for VIRIAL, fk component for IK
    int j = 0;                                   // virial component for VIRIAL
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const
    {
      const Complexd c = w1[i];
      if( mode == PPPMFactor::COPY ) w2[i] = c;
      else if( mode == PPPMFactor::VIRIAL ) { const double s = coef[6*i+j]; w2[i] = Complexd{ c.r * s , c.i * s }; }
      else { const double k = coef[i]; w2[i] = Complexd{ - k * c.i , k * c.r }; }
    }
  };

  struct PPPMRealPartFunc
  {
    const Complexd * __restrict__ w = nullptr;
    double * __restrict__ out = nullptr;
    ONIKA_HOST_DEVICE_FUNC inline void operator () ( ssize_t i ) const { out[i] = w[i].r; }
  };

}
}

namespace exanb
{
  template<bool PerAtomCharge> struct ComputeCellParticlesTraits< exaStamp::PPPMSpreadFunc<PerAtomCharge> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };

  template<bool PerAtomCharge, bool ComputeEnergy, bool ComputeVirial>
  struct ComputeCellParticlesTraits< exaStamp::PPPMForceFunc<PerAtomCharge,ComputeEnergy,ComputeVirial> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace onika
{
  namespace parallel
  {
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMZeroFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMLoadDensityFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMApplyGreenFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMFactorFunc > { static inline constexpr bool CudaCompatible = true; };
    template<> struct ParallelForFunctorTraits< exaStamp::PPPMRealPartFunc > { static inline constexpr bool CudaCompatible = true; };
  }
}

namespace exaStamp
{
inline namespace coulombic_ewald
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class PPPMLongRangePC : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( PPPMParameters  , pppm_config          , INPUT , OPTIONAL );
    ADD_SLOT( GridT           , grid                 , INPUT_OUTPUT );
    ADD_SLOT( Domain          , domain               , INPUT , REQUIRED );
    ADD_SLOT( double          , rcut_max             , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( ParticleSpecies , species              , INPUT , REQUIRED );
    ADD_SLOT( bool            , per_atom_charge      , INPUT , true , DocString{"read charges from per particle charge field instead of species charges"} );
    ADD_SLOT( MPI_Comm        , mpi                  , INPUT );
    ADD_SLOT( bool            , trigger_thermo_state , INPUT , OPTIONAL );

    // mesh work buffers, kept between time steps
    onika::memory::CudaMMVector<double> m_density;
    onika::memory::CudaMMVector<double> m_vd;     // 3 gradient meshes
    onika::memory::CudaMMVector<double> m_u;
    onika::memory::CudaMMVector<double> m_v;      // 6 virial meshes
    onika::memory::CudaMMVector<Complexd> m_work1;
    onika::memory::CudaMMVector<Complexd> m_work2;
    PPPMFFT m_fft;

    // per mesh point loop : GPU (one thread per point) when the FFT runs there, OpenMP otherwise
    template<class FuncT>
    inline void mesh_for( size_t n, const FuncT& func )
    {
      onika::parallel::ParallelForOptions opts;
      opts.enable_gpu = m_fft.on_gpu();
      onika::parallel::parallel_for( n , func , parallel_execution_context() , opts );
    }

    // m_work2 = m_work1 * factor (copy, virial coefficient or i.k), backward FFT, real part into out
    inline void backward_to( double * __restrict__ out, PPPMFactor mode, const double* coef = nullptr, int j = 0 )
    {
      const size_t nfft = m_work1.size();
      mesh_for( nfft , PPPMFactorFunc{ m_work1.data() , m_work2.data() , mode , coef , j } );
      m_fft.backward( m_work2.data() );
      mesh_for( nfft , PPPMRealPartFunc{ m_work2.data() , out } );
    }

  public:
    inline void execute () override final
    {
      const bool log_energy = trigger_thermo_state.has_value() ? *trigger_thermo_state : false;

      if( ! pppm_config.has_value() )
      {
        ldbg << "pppm_config not set, skip coulombic_pppm" << std::endl;
        return;
      }
      const auto & p = *pppm_config;
      *rcut_max = std::max( *rcut_max , p.radius );

      if( grid->number_of_cells() == 0 ) return;

      const Mat3d cell = ewald_cell_matrix( domain->xform() , domain->bounds_size() );
      if( p.volume == 0.0 || ! ewald_same_cell( cell , p.cell ) )
      {
        fatal_error() << "coulombic_pppm : domain cell "<<cell<<" differs from the one PPPM was set up for "
                      << p.cell << ". Call coulombic_pppm_init before force computation when the cell changes." << std::endl;
      }

      const size_t nfft = p.nfft();
      if( m_density.size() != nfft )
      {
        m_density.resize( nfft );
        m_vd.resize( 3*nfft );
        m_work1.resize( nfft );
        m_work2.resize( nfft );
        m_u.clear(); m_v.clear();
      }
      // GPU path when a device is available : particle kernels, mesh loops and cuFFT all run there (unified memory)
      const bool gpu_available = ( global_cuda_ctx() != nullptr ) && global_cuda_ctx()->has_devices() && PPPMFFT::gpu_support();
      void* stream = nullptr;
#     ifdef EXASTAMP_PPPM_CUFFT
      if( gpu_available ) stream = global_cuda_ctx()->getThreadStream(0);
#     endif
      m_fft.resize( p.nx , p.ny , p.nz , gpu_available , stream );
      if( log_energy && m_u.size() != nfft )
      {
        m_u.resize( nfft );
        m_v.resize( 6*nfft );
      }

      int nprocs = 1;
      MPI_Comm_size( *mpi , &nprocs );

      const ReadOnlyPPPMParameters ro( p , domain->bounds().bmin , domain->bounds_size() );

      auto rx = grid->field_accessor( field::rx );
      auto ry = grid->field_accessor( field::ry );
      auto rz = grid->field_accessor( field::rz );
      auto fx = grid->field_accessor( field::fx );
      auto fy = grid->field_accessor( field::fy );
      auto fz = grid->field_accessor( field::fz );
      auto ep = grid->field_accessor( field::ep );
      auto virial = grid->field_accessor( field::virial );

      auto compute_with_charges = [&]( auto per_atom_charge_tag , auto charge_or_type )
      {
        static constexpr bool PerAtomCharge = decltype(per_atom_charge_tag)::value;

        // 1. charge density of local particles, summed over all ranks (replicated global mesh)
        double * __restrict__ density = m_density.data();
        mesh_for( nfft , PPPMZeroFunc{ density } );
        PPPMSpreadFunc<PerAtomCharge> spread_func = { ro , species->data() , density };
        compute_cell_particles( *grid , false , spread_func , onika::make_flat_tuple(rx,ry,rz,charge_or_type) , parallel_execution_context() );
        if( nprocs > 1 ) MPI_Allreduce( MPI_IN_PLACE , density , nfft , MPI_DOUBLE , MPI_SUM , *mpi );

        // 2. rho(k), then V(k) = G(k) rho(k) / N (LAMMPS poisson_ik)
        mesh_for( nfft , PPPMLoadDensityFunc{ density , m_work1.data() } );
        m_fft.forward( m_work1.data() );
        mesh_for( nfft , PPPMApplyGreenFunc{ m_work1.data() , p.greensfn.data() , 1.0 / double(nfft) } );

        // 3. per particle energy and virial meshes (LAMMPS poisson_peratom)
        if( log_energy )
        {
          backward_to( m_u.data() , PPPMFactor::COPY );
          for( int j = 0 ; j < 6 ; j++ ) backward_to( m_v.data() + j*nfft , PPPMFactor::VIRIAL , p.vg.data() , j );
        }

        // 4. gradient of the potential, i.k V(k) back to real space
        backward_to( m_vd.data()          , PPPMFactor::IK , p.fkx.data() );
        backward_to( m_vd.data() +   nfft , PPPMFactor::IK , p.fky.data() );
        backward_to( m_vd.data() + 2*nfft , PPPMFactor::IK , p.fkz.data() );

        // 5. interpolate field (and energy, virial) back to local particles
        PPPMMeshes meshes = { m_vd.data() , m_vd.data() + nfft , m_vd.data() + 2*nfft , m_u.data() , m_v.data() };
        if( log_energy )
        {
          PPPMForceFunc<PerAtomCharge,true,true> force_func = { ro , species->data() , meshes , nfft };
          compute_cell_particles( *grid , false , force_func , onika::make_flat_tuple(fx,fy,fz,ep,virial,rx,ry,rz,charge_or_type) , parallel_execution_context() );
        }
        else
        {
          PPPMForceFunc<PerAtomCharge,false,false> force_func = { ro , species->data() , meshes , nfft };
          compute_cell_particles( *grid , false , force_func , onika::make_flat_tuple(fx,fy,fz,rx,ry,rz,charge_or_type) , parallel_execution_context() );
        }
      };

      if( *per_atom_charge ) compute_with_charges( std::true_type{}  , grid->field_accessor( field::charge ) );
      else                   compute_with_charges( std::false_type{} , grid->field_accessor( field::type ) );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Reciprocal space part of PPPM long range coulomb (ik differentiation), set up by coulombic_pppm_init.
Same algorithm as LAMMPS kspace_style pppm, orthogonal and triclinic cells. Use with coulombic_ewald_short_range for
the real space part. Computes forces, and when trigger_thermo_state is true, per particle energy (reciprocal + self +
neutralizing background) and per particle reciprocal virial.
The mesh is global and replicated on every MPI rank (density summed with MPI_Allreduce). When a GPU is available,
particle kernels, mesh loops and FFTs (cuFFT) run on the GPU ; otherwise on CPU with pocketfft.
)EOF";
    }
  };

  template<class GridT> using PPPMLongRangePCTmpl = PPPMLongRangePC<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(coulombic_pppm)
  {
    OperatorNodeFactory::instance()->register_factory( "coulombic_pppm" , make_grid_variant_operator<PPPMLongRangePCTmpl> );
  }

}
}
