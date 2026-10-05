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
#include <exanb/core/grid_fields.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>

#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <algorithm>
#include <string>
#include <vector>

#include "potential.h"
#include "k2b_descriptor_op.h"

// CPU/OpenMP-only descriptor-only pass for k2b, mirroring mlip-pod/compute_descriptor_pod.cu's
// shape: computes the coefficient-free per-atom 2-body kernel descriptor (and, optionally, its
// compact per-atom derivative aggregate) instead of the weighted energy/force
// k2b_compute_force_symetric computes -- no locks/scatter needed for the descriptor itself since
// each particle only ever writes its own output slot; the derivative aggregate uses
// atomic_add_contribution like POD/SNAP.
namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  class ComputeDescriptorK2b : public OperatorNode
  {
    ADD_SLOT( K2bPotentialParameters , parameters      , INPUT , REQUIRED );
    ADD_SLOT( double                 , rcut            , INPUT , REQUIRED );
    ADD_SLOT( double                 , rcut_max        , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( exanb::GridChunkNeighbors , chunk_neighbors , INPUT , exanb::GridChunkNeighbors{}, DocString{"neighbor list"} );
    ADD_SLOT( bool                   , ghost           , INPUT , false );
    ADD_SLOT( GridT                  , grid            , INPUT_OUTPUT );
    ADD_SLOT( Domain                 , domain          , INPUT , REQUIRED );

    ADD_SLOT( onika::memory::CudaMMVector<double>, k2b_descriptors, OUTPUT,
               DocString{"Flat per-particle k2b descriptor buffer: k2b_descriptors[ ncoeff*(cell_particle_offset[cell]+particle) + component ]"} );
    ADD_SLOT( long, ncoeff, OUTPUT, DocString{"Stride of the k2b_descriptors buffer = n_rbf"} );

    ADD_SLOT( bool , compute_derivative , INPUT , false ,
               DocString{"If true, also computes the compact per-atom descriptor-derivative aggregate: for atom a, n_rbf*3 values (k*3+xyz order), the sum over every atom i that has a as a neighbor (or i==a, the self term) of d(D_k of atom i)/d(r_a). Stored as dynamically-named generic-real grid fields (see deriv_agg_field_prefix), NOT a private buffer -- run update_opt_from_ghost on them on multi-rank runs."} );
    ADD_SLOT( std::string , deriv_agg_field_prefix , INPUT , std::string("k2bda_") ,
               DocString{"compute_derivative only: name prefix for the n_rbf*3 dynamically-named generic-real grid fields ('<prefix>0'..'<prefix>{n_rbf*3-1}') holding the derivative aggregate. KEEP THIS SHORT: dynamic field names are silently truncated to 15 characters + null (onika::soatl::FieldId's fixed char[16] m_name) -- this operator fatal_errors instead of silently colliding if prefix+max-index would overflow that limit."} );

    static constexpr bool UseWeights   = false;
    static constexpr bool UseNeighbors = true;
    using ComputeBuffer = ComputePairBuffer2<UseWeights, UseNeighbors>;
    static constexpr FieldSet<> compute_descriptor_field_set{};

  public:

    inline void execute() override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );

      *rcut_max = std::max( *rcut , *rcut_max );

      if (grid->number_of_cells() == 0) { *ncoeff = 0; return; }

      const int n_rbf = parameters->n_rbf;
      *ncoeff = n_rbf;

      const size_t total_particles = grid->number_of_particles();
      k2b_descriptors->clear();
      k2b_descriptors->resize( total_particles * n_rbf );
      std::fill( k2b_descriptors->begin(), k2b_descriptors->end(), 0.0 );

      onika::memory::CudaMMVector<double*> deriv_agg_ptrs;
      if (*compute_derivative)
      {
        const size_t nc3 = static_cast<size_t>(n_rbf) * 3;
        // onika::soatl::FieldId's dynamic-field name storage is a fixed char[16] (incl. null
        // terminator), silently strncpy-truncated -- a too-long prefix+index would alias
        // multiple components onto the same field with no error, so check instead of guessing.
        static constexpr size_t FIELD_NAME_MAX_LEN = 16;
        const size_t max_index_digits = std::to_string(nc3-1).size();
        if( deriv_agg_field_prefix->size() + max_index_digits + 1 > FIELD_NAME_MAX_LEN )
        {
          fatal_error() << "compute_descriptor_k2b: deriv_agg_field_prefix '"<<*deriv_agg_field_prefix<<"' is too long -- "
                        <<"prefix ("<<deriv_agg_field_prefix->size()<<" chars) + largest index ("<<max_index_digits<<" digits) "
                        <<"+ null terminator must fit in "<<FIELD_NAME_MAX_LEN<<" characters (onika::soatl::FieldId's fixed name buffer)" << std::endl;
        }
        deriv_agg_ptrs.resize( nc3 );
        for( size_t k=0; k<nc3; k++ )
        {
          double * const ptr = grid->flat_array_data( field::mk_generic_real( *deriv_agg_field_prefix + std::to_string(k) ) );
          std::fill_n( ptr, total_particles, 0.0 );
          deriv_agg_ptrs[k] = ptr;
        }
      }

      ComputePairNullWeightIterator cp_weight{};
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto descriptor_buf = make_compute_pair_buffer<ComputeBuffer>();
      LinearXForm cp_xform{ domain->xform() };
      ComputePairOptionalLocks<false> cp_locks{};

      K2bDescriptorOp descriptor_op{ K2bPotentialParametersRO(*parameters),
                                      grid->cell_particle_offset_data(), k2b_descriptors->data(),
                                      *compute_derivative ? deriv_agg_ptrs.data() : nullptr };
      compute_cell_particle_pairs(
          *grid, *rcut, *ghost,
          make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks),
          descriptor_buf, descriptor_op, compute_descriptor_field_set,
          parallel_execution_context());
    }

  };

  template<class GridT> using ComputeDescriptorK2bTmpl = ComputeDescriptorK2b<GridT>;

  ONIKA_AUTORUN_INIT(compute_descriptor_k2b)
  {
    OperatorNodeFactory::instance()->register_factory("compute_descriptor_k2b", make_grid_variant_operator<ComputeDescriptorK2bTmpl>);
  }

}
