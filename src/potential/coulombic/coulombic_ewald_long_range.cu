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

#include <exaStamp/potential/coulombic/ewald.h>
#include <exaStamp/unit_system.h>

#include <exanb/core/config.h>

#include <exanb/core/parallel_grid_algorithm.h>
#include <onika/cuda/cuda.h>
#include <exanb/core/xform.h>
#include <onika/flat_tuple.h>
#include <mpi.h>

namespace exaStamp
{
inline namespace coulombic_ewald // distinct symbols from the legacy ewald plugin (plugins are loaded RTLD_GLOBAL)
{
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;

  template<bool PerAtomCharge, class ChargeOrTypeT>
  ONIKA_HOST_DEVICE_FUNC static inline double ewald_particle_charge( const ParticleSpecie* __restrict__ species, ChargeOrTypeT ct )
  {
    if constexpr ( PerAtomCharge ) return ct;
    else return species[ct].m_charge;
  }

  // structure factor S(k) = sum_i q_i exp(i k.r_i)
  template<class XFormT, bool PerAtomCharge>
  struct EwaldLongRangeRhoComputeFunc
  {
    const XFormT xform;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    ReadOnlyEwaldParameters p;
    Complexd* __restrict__ m_ewald_rho = nullptr;
    
    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      const double q = ewald_particle_charge<PerAtomCharge>( m_species , ct );
      const Vec3d r = xform.transformCoord( Vec3d{rx,ry,rz} );
      const unsigned int nk = p.nknz;
      for(unsigned int k=0;k<nk;k++)
      {
        const EwaldCoeffs& gdata = p.Gdata[k];
        const double ps = r.x * gdata.Gx + r.y * gdata.Gy + r.z * gdata.Gz;
        double s,c;
        sincos(ps,&s,&c);
        // all threads of all blocks accumulate into the same S(k) : must be a true atomic on CPU too
        // (ONIKA_CU_BLOCK_ATOMIC_ADD is a plain += on CPU, racing between OpenMP threads)
        ONIKA_CU_ATOMIC_ADD( m_ewald_rho[k].r , q * c );
        ONIKA_CU_ATOMIC_ADD( m_ewald_rho[k].i , q * s );
      }
    }
  };

  // reciprocal forces, and optionaly per particle energy (reciprocal + self + background) and reciprocal virial.
  // per particle energy  : e_i = q_i sum_k Gc Re( conj(S(k)) exp(i k.r_i) ) , sums up to sum_k Gc |S(k)|^2
  // per particle virial  : W_i = e_i,k ( I - Gv G (x) G ) summed over k (same as LAMMPS Ewald per atom virial)
  template<class XFormT, bool PerAtomCharge, bool ComputeEnergy, bool ComputeVirial>
  struct EwaldLongRangeForceComputeFunc
  {
    const XFormT xform;
    const ParticleSpecie* __restrict__ m_species = nullptr;
    ReadOnlyEwaldParameters p;
    const Complexd* __restrict__ m_ewald_rho = nullptr;

    ONIKA_HOST_DEVICE_FUNC
    inline void compute ( double q, double rx, double ry, double rz, Vec3d& f, double& ep, Mat3d& vir ) const
    {
      const Vec3d r = xform.transformCoord( Vec3d{rx,ry,rz} );
      const unsigned int nk = p.nknz;
      const double q2 = 2. * q;
      for(unsigned int k=0;k<nk;k++)
      {
        const EwaldCoeffs& gdata = p.Gdata[k];
        const double ps = r.x * gdata.Gx + r.y * gdata.Gy + r.z * gdata.Gz;
        double s,c;
        sincos(ps,&s,&c);
        const double rr = m_ewald_rho[k].r;
        const double ri = m_ewald_rho[k].i;
        const double al = q2 * gdata.Gc * ( rr * s - ri * c );
        f.x += al * gdata.Gx;
        f.y += al * gdata.Gy;
        f.z += al * gdata.Gz;
        if constexpr ( ComputeEnergy )
        {
          const double ek = q * gdata.Gc * ( rr * c + ri * s );
          ep += ek;
          if constexpr ( ComputeVirial )
          {
            const double w = ek * gdata.Gv;
            vir.m11 += ek - w * gdata.Gx * gdata.Gx;
            vir.m22 += ek - w * gdata.Gy * gdata.Gy;
            vir.m33 += ek - w * gdata.Gz * gdata.Gz;
            const double vxy = - w * gdata.Gx * gdata.Gy;
            const double vxz = - w * gdata.Gx * gdata.Gz;
            const double vyz = - w * gdata.Gy * gdata.Gz;
            vir.m12 += vxy; vir.m21 += vxy;
            vir.m13 += vxz; vir.m31 += vxz;
            vir.m23 += vyz; vir.m32 += vyz;
          }
        }
      }
      if constexpr ( ComputeEnergy )
      {
        ep += ewald_self_energy( p , q );
      }
    }

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( double & fx, double & fy, double & fz, double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      static_assert( ! ComputeEnergy && ! ComputeVirial );
      Vec3d f = {0.,0.,0.}; double ep = 0.0; Mat3d vir;
      compute( ewald_particle_charge<PerAtomCharge>( m_species , ct ) , rx, ry, rz, f, ep, vir );
      fx += f.x; fy += f.y; fz += f.z;
    }

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( double & fx, double & fy, double & fz, double & ep, double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      static_assert( ComputeEnergy && ! ComputeVirial );
      Vec3d f = {0.,0.,0.}; double e = 0.0; Mat3d vir;
      compute( ewald_particle_charge<PerAtomCharge>( m_species , ct ) , rx, ry, rz, f, e, vir );
      fx += f.x; fy += f.y; fz += f.z; ep += e;
    }

    template<class ChargeOrTypeT>
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( double & fx, double & fy, double & fz, double & ep, Mat3d & virial, double rx, double ry, double rz, ChargeOrTypeT ct ) const
    {
      static_assert( ComputeEnergy && ComputeVirial );
      Vec3d f = {0.,0.,0.}; double e = 0.0; Mat3d vir = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
      compute( ewald_particle_charge<PerAtomCharge>( m_species , ct ) , rx, ry, rz, f, e, vir );
      fx += f.x; fy += f.y; fz += f.z; ep += e; virial += vir;
    }
  };
  
}
}

namespace exanb
{
  template<class XFormT, bool PerAtomCharge> struct ComputeCellParticlesTraits< exaStamp::EwaldLongRangeRhoComputeFunc<XFormT,PerAtomCharge> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };

  template<class XFormT, bool PerAtomCharge, bool ComputeEnergy, bool ComputeVirial>
  struct ComputeCellParticlesTraits< exaStamp::EwaldLongRangeForceComputeFunc<XFormT,PerAtomCharge,ComputeEnergy,ComputeVirial> >
  {
    static inline constexpr bool RequiresBlockSynchronousCall = false;
    static inline constexpr bool CudaCompatible = true;
  };
}

namespace exaStamp
{
inline namespace coulombic_ewald // distinct symbols from the legacy ewald plugin (plugins are loaded RTLD_GLOBAL)
{
  using namespace exanb;

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
  class EwaldLongRangePC : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( EwaldParameters , ewald_config , INPUT , OPTIONAL );
    ADD_SLOT( GridT      , grid         , INPUT_OUTPUT );
    ADD_SLOT( Domain     , domain       , INPUT , REQUIRED );
    ADD_SLOT( double     , rcut_max     , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( EwaldRho   , ewald_rho    , INPUT_OUTPUT , EwaldRho{} );
    ADD_SLOT( ParticleSpecies  , species           , INPUT , REQUIRED );
    ADD_SLOT( bool       , per_atom_charge , INPUT , true , DocString{"read charges from per particle charge field instead of species charges"} );
    ADD_SLOT( MPI_Comm   , mpi          , INPUT );
    ADD_SLOT( bool       , trigger_thermo_state    , INPUT , OPTIONAL );

  public:
    // Operator execution
    inline void execute () override final
    {
      bool log_energy = false;
      if( trigger_thermo_state.has_value() )
      {
        log_energy = *trigger_thermo_state ;
      }
      else
      {
        ldbg << "trigger_thermo_state missing " << std::endl;
      }
      
      if( ! ewald_config.has_value() )
      {
        ldbg << "ewald_config not set, skip ewal_long_range computation" << std::endl;
        return ;
      }

      *rcut_max = std::max( *rcut_max , ewald_config->radius );
      
      if( grid->number_of_cells() == 0 ) return;

      // k vectors depend on the cell, refuse to compute with stale ones
      const Mat3d cell = ewald_cell_matrix( domain->xform() , domain->bounds_size() );
      if( ! ewald_same_cell( cell , ewald_config->cell ) )
      {
        fatal_error() << "coulombic_ewald_long_range : domain cell "<<cell<<" differs from the one used to build k vectors "
                      << ewald_config->cell << ". Call coulombic_ewald_init before force computation when the cell changes." << std::endl;
      }

      const size_t nk = ewald_config->nknz;
      if (ewald_rho->rho.size() < nk) {
        ewald_rho->rho.resize(nk);
      }
      ewald_rho->nk = nk;
      
      const Mat3d xform = domain->xform();
      const ReadOnlyEwaldParameters ro_params = *ewald_config;
      const bool gpu_available = ( global_cuda_ctx() != nullptr ) && global_cuda_ctx()->has_devices();
      if( gpu_available )
      {
        ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MEMSET( ewald_rho->rho.data(), 0, sizeof(Complexd)*nk, global_cuda_ctx()->getThreadStream(0) ) );
      }
      else
      {
        for(size_t k=0;k<nk;k++) ewald_rho->rho[k] = Complexd{ 0.0 , 0.0 };
      }

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

        EwaldLongRangeRhoComputeFunc<LinearXForm,PerAtomCharge> rho_func = { {xform} , species->data() , ro_params , ewald_rho->rho.data() };
        compute_cell_particles( *grid , false , rho_func , onika::make_flat_tuple(rx,ry,rz,charge_or_type) , parallel_execution_context() );
        static_assert( sizeof(Complexd) == 2*sizeof(double) );
        MPI_Allreduce(MPI_IN_PLACE, (double*) ewald_rho->rho.data(),nk*2,MPI_DOUBLE,MPI_SUM,*mpi);

        if( log_energy )
        {
          EwaldLongRangeForceComputeFunc<LinearXForm,PerAtomCharge,true,true> force_func = { {xform} , species->data() , ro_params , ewald_rho->rho.data() };
          compute_cell_particles( *grid , false , force_func , onika::make_flat_tuple(fx,fy,fz,ep,virial,rx,ry,rz,charge_or_type) , parallel_execution_context() );
        }
        else
        {
          EwaldLongRangeForceComputeFunc<LinearXForm,PerAtomCharge,false,false> force_func = { {xform} , species->data() , ro_params , ewald_rho->rho.data() };
          compute_cell_particles( *grid , false , force_func , onika::make_flat_tuple(fx,fy,fz,rx,ry,rz,charge_or_type) , parallel_execution_context() );
        }
      };

      if( *per_atom_charge ) compute_with_charges( std::true_type{}  , grid->field_accessor( field::charge ) );
      else                   compute_with_charges( std::false_type{} , grid->field_accessor( field::type ) );
    }

    inline std::string documentation() const override final
    {
      return R"EOF(
Reciprocal space part of the Ewald summation (direct sum over k vectors built by coulombic_ewald_init).
Computes forces, and when trigger_thermo_state is true, per particle energy (reciprocal + self + neutralizing background)
and per particle reciprocal virial.
)EOF";
    }

  };
  
  template<class GridT> using EwaldLongRangePCTmpl = EwaldLongRangePC<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_ewald_long_range)
  {  
    OperatorNodeFactory::instance()->register_factory( "coulombic_ewald_long_range" , make_grid_variant_operator<EwaldLongRangePCTmpl> );
  }

}
}
