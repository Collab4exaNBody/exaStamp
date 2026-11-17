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

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exaStamp/mlip_ye3t/ye3t_context.h>
#include <exaStamp/mlip_ye3t/ye3t_driver.h>
#include <exaStamp/mlip_utils/mlip_flat_bonds.h>
#include <exaStamp/mlip_utils/mlip_onika_launch.h>

#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/cuda/cuda_context.h>

#include <omp.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace exaStamp
{
  using namespace exanb;

  namespace mlip_ye3t
  {
    // GPU code paths run on the host (ye3t_emul.cpp)
    void ye3t_onika_emul( const Ye3tTables& T, Ye3tBuffers& buf, const mlip_utils::BondsView& bv, mlip_utils::OutputsView& out,
                          int nthreads, double mem_gb, int& n_rho_limit );
    void ye3t_tagged_emul( const Ye3tTaggedTables& T, Ye3tTaggedBuffers& buf, const mlip_utils::BondsView& bv, mlip_utils::OutputsView& out,
                           int nthreads, double mem_gb );
  }

  template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz, field::_type>>
  class Ye3tForce : public OperatorNode
  {
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT        , exanb::GridChunkNeighbors{} , DocString{"neighbor list (full, not symmetric)"} );
    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain          , INPUT        , REQUIRED );
    ADD_SLOT( Ye3tContext               , ye3t_ctx        , INPUT_OUTPUT , REQUIRED , DocString{"model built by ye3t_init"} );

    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

Linear atomic cluster expansion (.yace family of ye3t-lammps pair_style ye3t), native onika kernels: GPU if
available, else OpenMP. The model comes from ye3t_init. Central atoms and their ghost neighbours are gathered into
flat arrays, the kernels give forces on every atom, ghosts included, so forces and energies must be zeroed in the
ghost layer before and ghost contributions folded back afterwards:

compute_force_prolog:
  - zero_force_energy: { ghost: true }
compute_force_epilog:
  - update_force_energy_from_ghost
  - force_to_accel

compute_force:
  - ye3t_force

A tagged-Cauchy correction (composite model of ye3t_init) is computed after the linear ACE backbone, on the same bonds.

Per-atom energy: the atomic energy E_i of the model. Per-atom virial: the whole virial of each bond on its central atom.

)EOF";
    }

    inline void execute() override final
    {
      Ye3tContext& ctx = *ye3t_ctx;
      if( !ctx.tables && !ctx.tagged ) fatal_error() << "ye3t_force: no model, ye3t_init must run first" << std::endl;
      if( grid->number_of_cells() == 0 ) return;
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      const int nth = omp_get_max_threads();
      mlip_utils::BondsView bv;
      mlip_utils::OutputsView out;
      static const bool prof = std::getenv( "YE3T_PROF" ) != nullptr;
      static double tg = 0.0, ts = 0.0;
      auto T0 = std::chrono::steady_clock::now();
      mlip_utils::gather_bonds( *grid, *chunk_neighbors, *domain, ctx.species_type, ctx.rcut, has_virial_field, ctx.flat, bv, out, nth );
      if( prof ) tg += std::chrono::duration<double>( std::chrono::steady_clock::now() - T0 ).count();

      const bool gpu = onika::parallel::functor_gpu_support_v< mlip_ye3t::ye3t_onika::BasisFunctor >
                    && global_cuda_ctx() != nullptr && global_cuda_ctx()->has_devices();
      const mlip_utils::OnikaLauncher<Ye3tForce> pec{ this };
      int n_rho = 0;
      if( ctx.tables )
      {
        const mlip_ye3t::Ye3tTables& T = *ctx.tables;
        if( ctx.emulate > 0 ) mlip_ye3t::ye3t_onika_emul( T, ctx.buf, bv, out, ctx.emulate, ctx.mem_gb, n_rho );
        else mlip_ye3t::ye3t_onika_compute( T, ctx.buf, bv, out, gpu, ctx.mem_gb, pec, n_rho );
      }
      if( n_rho > 0 ) fatal_error() << "ye3t_force: " << n_rho << " atom(s) with a density beyond rho_core_cutoff - drho_core_cutoff (outside the validity of the linear model)" << std::endl;
      if( ctx.tagged )
      {
        if( ctx.emulate > 0 ) mlip_ye3t::ye3t_tagged_emul( *ctx.tagged, ctx.tbuf, bv, out, ctx.emulate, ctx.mem_gb );
        else mlip_ye3t::ye3t_tagged_compute( *ctx.tagged, ctx.tbuf, bv, out, gpu, ctx.mem_gb, pec );
      }
      T0 = std::chrono::steady_clock::now();
      mlip_utils::scatter_outputs( *grid, ctx.flat, nth );
      if( prof ) { ts += std::chrono::duration<double>( std::chrono::steady_clock::now() - T0 ).count(); std::fprintf( stderr, "YE3T_PROF gather %.3f s scatter %.3f s\n", tg, ts ); }
    }
  };

  template<class GridT> using Ye3tForceTmpl = Ye3tForce<GridT>;

  ONIKA_AUTORUN_INIT(ye3t_force)
  {
    OperatorNodeFactory::instance()->register_factory( "ye3t_force", make_grid_variant_operator<Ye3tForceTmpl> );
  }

}
