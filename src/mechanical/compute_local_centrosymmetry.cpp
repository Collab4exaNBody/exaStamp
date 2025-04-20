#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>
#include <onika/yaml/yaml_enum.h>
#include <onika/string_utils.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/mechanical/compute_local_centrosymmetry.h>

namespace exaStamp
{
using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;

template< class GridT,
    class = AssertGridHasFields< GridT, field::_ep ,field::_fx ,field::_fy ,field::_fz >
    >
class ComputeCSPOperator : public OperatorNode
{    

  // ========= I/O slots =======================
  //  ADD_SLOT( double                , rcut_max            , INPUT_OUTPUT , 0.0 );
  ADD_SLOT( exanb::GridChunkNeighbors                   , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
  ADD_SLOT( bool                  , ghost               , INPUT , false );
  ADD_SLOT( GridT                 , grid                , INPUT );
  ADD_SLOT( Domain                , domain              , INPUT , REQUIRED );

  ADD_SLOT(double, rcut, INPUT, 0.);
  ADD_SLOT(uint64_t, nnn, INPUT, 12);
  ADD_SLOT(LocalCentroMethod, method, INPUT, LocalCentroMethod::GES);

  using ComputeBuffer = ComputePairBuffer2<false, false>;
  using ComputeFields = FieldSet<>;
  static constexpr ComputeFields compute_force_field_set{};

public:

  inline void execute () override final
  {
    assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
    bool has_chunk_neighbors = chunk_neighbors.has_value();
    if( !has_chunk_neighbors )
    {
      lerr << "No neighbors input data available" << std::endl;
      std::abort();
    }

    double cutoff = *rcut;
    lout << onika::format_string("\t- Computing per-atom centrosymmetry:\n")
         << onika::format_string("\t    rcut = %.5f\n", cutoff)
         << onika::format_string("\t    nnn  = %d\n", *nnn);

    // check if number of nearest neighbor is even
    if (( *nnn % 2 ) != 0) {
      lerr << "ComputeCentrosymmetryOperator: the number of nearest neighbors (nnn = " << *nnn << ") must be even" << std::endl;
      std::abort();
    }

    CentroSymmetryOp local_op{cutoff, *nnn, *method};

    auto field_csp = grid->field_accessor(field::csp);
    auto local_op_fields = make_field_tuple_from_field_set(FieldSet<>{}, field_csp);
    
    ComputePairNullWeightIterator cp_weight{};
    ComputePairOptionalLocks<false> cp_locks{};
    exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{*chunk_neighbors};
    auto local_op_buf = make_compute_pair_buffer<ComputeBuffer>();
    ComputePairTrivialCellFiltering cpu_cell_filter = {};

    if (domain->xform_is_identity()) {
      NullXForm cp_xform;
      auto optional = make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter);
      compute_cell_particle_pairs(*grid, cutoff, *ghost, optional, local_op_buf, local_op, local_op_fields,
                                  parallel_execution_context());
    } else {
      LinearXForm cp_xform{domain->xform()};
      auto optional = make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter);
      compute_cell_particle_pairs(*grid, cutoff, *ghost, optional, local_op_buf, local_op, local_op_fields,
                                  parallel_execution_context());
    }

    lout << onika::format_string("\t- Computing per-atom centrosymmetry END\n");    
  }

};

template <class GridT> using ComputeCSPOperatorTmpl = ComputeCSPOperator<GridT>;

// === register factories ===
ONIKA_AUTORUN_INIT(compute_local_csp) {
  OperatorNodeFactory::instance()->register_factory("compute_local_csp",
                                                    make_grid_variant_operator<ComputeCSPOperatorTmpl>);
}
} // namespace exaStamp
