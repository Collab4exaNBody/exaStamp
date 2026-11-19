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

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exaStamp/unit_system.h>

#include <onika/log.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include "deepmd_force_op.h"
#include "deepmd_atten_op.h"
#include "deepmd_hybrid_op.h"
#include "deepmd_dpa3_driver.h"
#include "deepmd_dpa2.h"
#include "deepmd_dpa4.h"
#include "deepmd_model.h"

namespace exaStamp
{
  using namespace exanb;

  template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz, field::_type>>
  class DeepmdForce : public OperatorNode
  {
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT        , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
    ADD_SLOT( bool                      , ghost           , INPUT        , false , DocString{"compute central atoms in the ghost layer too (not needed with ghost force fold-back)"} );
    ADD_SLOT( GridT                     , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain                    , domain          , INPUT        , REQUIRED );
    ADD_SLOT( DeepmdContext             , deepmd_ctx      , INPUT        , REQUIRED , DocString{"DeePMD model built by deepmd_init"} );
    ADD_SLOT( GridParticleLocks         , particle_locks  , INPUT        , OPTIONAL , DocString{"particle spin locks (unused, updates are atomic)"} );
    ADD_SLOT( MPI_Comm                  , mpi             , INPUT        , MPI_COMM_WORLD );
    ADD_SLOT( GhostCommunicationScheme  , ghost_comm_scheme , INPUT      , OPTIONAL , DocString{"ghost scheme of the last ghost update (DPA3 / DPA2: node features exchanged between layers)"} );
    ADD_SLOT( double                    , dpa3_store_gb   , INPUT        , 4.0 , DocString{"DPA3: keep the edge / angle features of every layer (about 2 MB per atom for the 128/64/32 example model) instead of recomputing them, up to this size in GB per process (0: always recompute)"} );

    static constexpr bool UseWeights   = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights, UseNeighbors, DeepmdComputeBuffer, DefaultComputePairBufferAppendFunc,
                                             DP_MAX_NEIGHBORS, ComputePairBuffer2Weights, FieldSet<field::_type> >;

#   if defined(__CUDACC__) || defined(__HIPCC__)
    // one block per central atom: the compute buffer, with the per-atom scratch, is block shared memory
    static_assert( sizeof(ComputeBuffer) <= 46 * 1024 , "DeePMD per-atom scratch exceeds the static shared memory of a GPU block: "
                   "lower XSTAMP_DEEPMD_MAX_NEIGHBORS, XSTAMP_DEEPMD_MAX_DESCRIPTOR or XSTAMP_DEEPMD_MAX_FIT_* " );
#   endif

    static constexpr bool has_virial_field = GridHasField<GridT, field::_virial>::value;
    using ComputeFieldsWithoutVirial = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_type >;
    using ComputeFieldsWithVirial    = FieldSet< field::_ep, field::_fx, field::_fy, field::_fz, field::_type, field::_virial >;
    using ComputeFields = std::conditional_t< has_virial_field, ComputeFieldsWithVirial, ComputeFieldsWithoutVirial >;
    static constexpr ComputeFields compute_force_field_set{};

  public:

    inline std::string documentation() const override final
    {
      return R"EOF(

DeePMD machine-learning potential (se_e2_a, se_e2_r, se_e3, se_atten, se_t_tebd, hybrid, dpa3, dpa2, dpa4 descriptors), evaluated natively: no DeePMD-kit or
libtorch library at run time. The model comes from deepmd_init. On GPU, one block of threads
computes one atom (neighbors over threads for the embedding networks, neurons over threads for the
fitting network); model sizes are bounded by the XSTAMP_DEEPMD_MAX_* build options. Forces on neighbors are atomically
scattered, ghost particles included, so forces and energies must be zeroed in the ghost layer before
the force computation and ghost contributions folded back afterwards :

compute_force_prolog:
  - zero_force_energy: { ghost: true }
compute_force_epilog:
  - update_force_energy_from_ghost
  - force_to_accel

compute_force:
  - deepmd_force

Per-atom energy : the atomic energy of the model.

)EOF";
    }

    inline void execute() override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      if( grid->number_of_cells() == 0 ) return;

      if( deepmd_ctx->ntypes <= 0 || ( deepmd_ctx->fitting.empty() && deepmd_ctx->descriptor != DP_DPA4 ) )
      {
        fatal_error() << "deepmd_force: the DeePMD model is empty, deepmd_init must run first" << std::endl;
      }

      static constexpr double conv_energy = EXASTAMP_CONST_QUANTITY( 1. * eV );
      DeepmdForceOp force_op{ DeepmdModelView( *deepmd_ctx, conv_energy ) };
      const DeepmdContext& c = *deepmd_ctx;
      if( c.descriptor == DP_DPA4 )
      {
        dpa4_compute( *grid, *chunk_neighbors, *domain, *c.dpa4, DeepmdModelView( c, conv_energy ),
                      ghost_comm_scheme.has_value() ? ghost_comm_scheme.get_pointer() : nullptr, *mpi, *c.dpa4_buf );
        return;
      }
      if( c.descriptor == DP_DPA3 )
      {
        dpa3_compute( *grid, *chunk_neighbors, *domain, *c.dpa3, DeepmdModelView( c, conv_energy ),
                      ghost_comm_scheme.has_value() ? ghost_comm_scheme.get_pointer() : nullptr, *mpi, *c.dpa3_buf, *dpa3_store_gb );
        return;
      }

      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      ComputePairNullWeightIterator cp_weight{};
      LinearXForm cp_xform{ domain->xform() };
      ComputePairOptionalLocks<false> cp_locks{};
      auto force_buf = make_compute_pair_buffer<ComputeBuffer>();

      // neighbor species are packed into the compute buffer (buf.nbh_pt)
      auto optional = make_compute_pair_optional_args( nbh_it, cp_weight, cp_xform, cp_locks,
                                                       ComputePairTrivialCellFiltering{}, ComputePairTrivialParticleFiltering{},
                                                       grid->field_accessors_from_field_set( FieldSet<field::_type>{} ) );
      auto cp_fields = grid->field_accessors_from_field_set( compute_force_field_set );
      auto atten_view = [](const DeepmdContext& c)
      {
        DeepmdAttenView av;
        av.layers = c.attn_layers.data(); av.nlayers = c.attn_layers.size(); av.hidden = c.attn_hidden; av.tebd_dim = c.tebd_dim;
        av.scaling = c.attn_scaling; av.shift = c.attn_shift; av.normalize = c.attn_normalize; av.smooth = c.attn_smooth;
        av.dotr = c.attn_dotr; av.strip = c.strip; av.concat_output_tebd = c.concat_output_tebd;
        av.tebd = c.tebd.data(); av.strip_table = c.strip_table.data();
        return av;
      };
      if( c.descriptor == DP_DPA2 )
      {
        // pass A: repinit / three-body descriptors of the central atoms; repformers; pass C: their backward from dE/dD0
        const DPA2Model& M = *c.dpa2;
        DPA2Buffers& BB = *c.dpa2_buf;
        const size_t n_nodes = grid->number_of_particles();
        DPA2Timer TW;
        dpa2_prof().evals ++;
        BB.D0.assign( n_nodes * M.dim0, 0.0 );
        DPA2BlockOp op;
        op.rep = DeepmdAttenOp{ DeepmdForceOp{ DeepmdModelView( *M.rep, conv_energy ) }, atten_view( *M.rep ) };
        if( M.three ) { op.three = DeepmdForceOp{ DeepmdModelView( *M.three, conv_energy ) }; op.has_three = true; }
        op.dim_rep = M.dim_rep; op.dim_3b = M.dim_3b; op.dim0 = M.dim0;
        op.poff = grid->cell_particle_offset_data();
        op.D0 = BB.D0.data();
        const double rc0 = std::max( M.rep->rcut, M.three ? M.three->rcut : 0.0 );
        compute_cell_particle_pairs2( *grid, rc0, *ghost, optional, force_buf, op, cp_fields, DefaultPositionFields{}, parallel_execution_context() );
        TW.mark(0);
        dpa2_compute( *grid, *chunk_neighbors, *domain, M, DeepmdModelView( c, conv_energy ),
                      ghost_comm_scheme.has_value() ? ghost_comm_scheme.get_pointer() : nullptr, *mpi, BB );
        TW = DPA2Timer();
        op.backward = true;
        op.dD0 = BB.dD0.data();
        compute_cell_particle_pairs2( *grid, rc0, *ghost, optional, force_buf, op, cp_fields, DefaultPositionFields{}, parallel_execution_context() );
        TW.mark(5);
        return;
      }
      if( c.descriptor == DP_HYBRID )
      {
        std::vector<DeepmdForceOp> fops;
        std::vector<DeepmdAttenOp> aops;
        std::vector<int> kind, idx, dim;
        for(const auto& sub : c.subs)
        {
          DeepmdForceOp op{ DeepmdModelView( *sub, conv_energy ) };
          kind.push_back( sub->descriptor == DP_SE_ATTEN ? 1 : 0 );
          if( kind.back() == 1 ) { idx.push_back( aops.size() ); aops.push_back( DeepmdAttenOp{ op, atten_view(*sub) } ); }
          else { idx.push_back( fops.size() ); fops.push_back( op ); }
          dim.push_back( sub->dim );
        }
        DeepmdHybridOp hybrid_op{ DeepmdModelView( c, conv_energy ), fops.data(), aops.data(), kind.data(), idx.data(), dim.data(), int(c.subs.size()) };
        compute_cell_particle_pairs2( *grid, deepmd_ctx->rcut, *ghost, optional, force_buf, hybrid_op, cp_fields,
                                      DefaultPositionFields{}, parallel_execution_context() );
        return;
      }
      if( c.descriptor == DP_SE_ATTEN )
      {
        DeepmdAttenOp atten_op{ force_op, atten_view(c) };
        compute_cell_particle_pairs2( *grid, deepmd_ctx->rcut, *ghost, optional, force_buf, atten_op, cp_fields,
                                      DefaultPositionFields{}, parallel_execution_context() );
        return;
      }
      compute_cell_particle_pairs2( *grid, deepmd_ctx->rcut, *ghost, optional, force_buf, force_op, cp_fields,
                                    DefaultPositionFields{}, parallel_execution_context() );
    }
  };

  template<class GridT> using DeepmdForceTmpl = DeepmdForce<GridT>;

  ONIKA_AUTORUN_INIT(deepmd_force)
  {
    OperatorNodeFactory::instance()->register_factory( "deepmd_force", make_grid_variant_operator<DeepmdForceTmpl> );
  }

}
