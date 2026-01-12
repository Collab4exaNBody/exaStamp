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
#include <mpi.h>

namespace exaStamp
{
  using ewald_constants::fpe0;
  using ewald_constants::epsilonZero;  
  using namespace exanb;

  using onika::memory::DEFAULT_ALIGNMENT;
  
  template<class XFormT>
  struct EwaldLongRangeRhoComputeFunc
  {
    const XFormT xform;
    const ParticleSpecie * __restrict__ m_species = nullptr;
    ReadOnlyEwaldParameters p;
    Complexd* __restrict__ m_ewald_rho = nullptr;
    
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( double rx, double ry, double rz, unsigned int type ) const
    {
      const double q = m_species[type].m_charge;
      const Vec3d r = xform.transformCoord( Vec3d{rx,ry,rz} );
      const unsigned int nk = p.nknz;
      for(unsigned int k=0;k<nk;k++)
      {
        const auto gdata = p.Gdata[k];
        const double ps = r.x * gdata.Gx + r.y * gdata.Gy + r.z * gdata.Gz;
        ONIKA_CU_BLOCK_ATOMIC_ADD( m_ewald_rho[k].r , q * cos(ps) );
        ONIKA_CU_BLOCK_ATOMIC_ADD( m_ewald_rho[k].i , q * sin(ps) );
      }
    }
  };

  template<class XFormT>
  struct EwaldLongRangeForceComputeFunc
  {
    const XFormT xform;
    const ParticleSpecie* __restrict__ m_species = nullptr;
    ReadOnlyEwaldParameters p;
    const Complexd* __restrict__ m_ewald_rho = nullptr;
    
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( double & fx, double & fy, double & fz, double rx, double ry, double rz, unsigned int type  ) const
    {
      const double q = 2. * m_species[type].m_charge;
      const Vec3d r = xform.transformCoord( Vec3d{rx,ry,rz} );
      const unsigned int nk = p.nknz;
      double lfx = 0.0;
      double lfy = 0.0;
      double lfz = 0.0;
      for(unsigned int k=0;k<nk;k++)
      {
        const auto gdata = p.Gdata[k];
        const double ps = r.x * gdata.Gx + r.y * gdata.Gy + r.z * gdata.Gz;
        const double al = q * gdata.Gc * ( m_ewald_rho[k].r * sin(ps) - m_ewald_rho[k].i * cos(ps) );
        lfx += al * gdata.Gx;
        lfy += al * gdata.Gy;
        lfz += al * gdata.Gz;
      }
      fx += lfx;
      fy += lfy;
      fz += lfz;
    }
  };

  struct EwaldSelfForceComputeFunc
  {
    const ParticleSpecie* __restrict__ m_species = nullptr;
    ReadOnlyEwaldParameters p;
    
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( double & ep, unsigned int type  ) const
    {
      const double q = m_species[type].m_charge;
      ep -= 1. / fpe0 * p.g_ewald / std::sqrt(M_PI) * q * q;
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
  
  template<> struct ComputeCellParticlesTraits< exaStamp::EwaldSelfForceComputeFunc >
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
    ADD_SLOT( MPI_Comm   , mpi          , INPUT );
    ADD_SLOT( bool       , trigger_thermo_state    , INPUT , OPTIONAL );
    ADD_SLOT( double     , potential_energy_shift  , OUTPUT );

    using ewald_rho_field_set_t = FieldSet<field::_rx,field::_ry,field::_rz,field::_type>;
    using ewald_force_field_set_t = FieldSet<field::_fx,field::_fy,field::_fz,field::_rx,field::_ry,field::_rz,field::_type>;
    using ewald_self_force_field_set_t = FieldSet<field::_ep,field::_type>;    
    
    static constexpr ewald_rho_field_set_t ewald_rho_field_set = {};
    static constexpr ewald_force_field_set_t ewald_force_field_set = {};
    static constexpr ewald_self_force_field_set_t ewald_self_force_field_set = {};

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

      // Recupération du nombre de point k
      const size_t nk = ewald_config->nknz;
      ewald_rho->nk = nk;
      ewald_rho->rho.resize(nk);
      
      //      const EwaldCoeffs * __restrict__ Gdata = ewald_config->Gdata.data();
      Mat3d xform = domain->xform();
      ONIKA_CU_CHECK_ERRORS( ONIKA_CU_MEMSET( ewald_rho->rho.data(), 0, sizeof(Complexd)*nk, global_cuda_ctx()->getThreadStream(0) ) );

      EwaldLongRangeRhoComputeFunc<LinearXForm> rho_func = { {xform} , species->data() , *ewald_config , ewald_rho->rho.data() };
      compute_cell_particles( *grid , false , rho_func , ewald_rho_field_set , parallel_execution_context() );
      static_assert( sizeof(Complexd) == 2*sizeof(double) );
      MPI_Allreduce(MPI_IN_PLACE, (double*) ewald_rho->rho.data(),nk*2,MPI_DOUBLE,MPI_SUM,*mpi);

      EwaldLongRangeForceComputeFunc<LinearXForm> force_func = { {xform} , species->data() , *ewald_config , ewald_rho->rho.data() };
      compute_cell_particles( *grid , false , force_func , ewald_force_field_set , parallel_execution_context() );

      // Self energy + reciprocal energy are computed only and if only trigger_thermo_state = true
      if ( log_energy ) {
        *potential_energy_shift = 0.;
        EwaldSelfForceComputeFunc self_func = { species->data() , *ewald_config };
        compute_cell_particles( *grid , false , self_func , ewald_self_force_field_set , parallel_execution_context() );
        const size_t nk = ewald_rho->nk;
        double re = 0.;
  #     pragma omp parallel for schedule(static) reduction(+:re)
        for (size_t k=0; k<nk; ++k)
          {
            re +=  ewald_config->Gdata[k].Gc * complex_norm( ewald_rho->rho[k] ); // (totalRho_r[k]*totalRho_r[k] + totalRho_i[k]*totalRho_i[k]);
          }
        *potential_energy_shift += re;
      }
      
      ldbg<<"------------------------"<<std::endl<<std::flush;
      ldbg<<"End of long range energy"<<std::endl<<std::flush;
      ldbg<<"------------------------"<<std::endl<<std::flush;
    }

  };
  
  template<class GridT> using EwaldLongRangePCTmpl = EwaldLongRangePC<GridT>;

  // === register factories ===  
  ONIKA_AUTORUN_INIT(coulombic_ewald_long_rage)
  {  
    OperatorNodeFactory::instance()->register_factory( "coulombic_ewald_long_range" , make_grid_variant_operator<EwaldLongRangePCTmpl> );
  }

}


