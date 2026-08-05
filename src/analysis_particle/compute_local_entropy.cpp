#include <exanb/core/make_grid_variant_operator.h>
#include <onika/cpp_utils.h>
#include <onika/file_utils.h>
#include <onika/log.h>
#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/string_utils.h>

#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/particle_species/particle_specie.h>

#include <exaStamp/mechanical/compute_local_entropy.h>

namespace exaStamp
{

using namespace exanb;

template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz>>
class ComputeLocalEntropyOp : public OperatorNode {

  // ========= I/O slots =======================
  ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
  ADD_SLOT(bool, ghost, INPUT, false);
  ADD_SLOT(GridT, grid, INPUT);
  ADD_SLOT(Domain, domain, INPUT, REQUIRED);
  ADD_SLOT(ThermodynamicState, thermodynamic_state, INPUT, REQUIRED);

  ADD_SLOT(double, rcut, INPUT, 0.);
  ADD_SLOT(double, sigma, INPUT, REQUIRED);
  ADD_SLOT(size_t, nbins, INPUT, 0);
  ADD_SLOT(bool, local, INPUT, false);

  ADD_SLOT(std::string, name, INPUT, "entropy");
  ADD_SLOT(bool, verbose, INPUT, false);

  using ComputeBuffer = ComputePairBuffer2<false, false>;
  using ComputeFields = FieldSet<>;
  static constexpr ComputeFields compute_force_field_set{};

  static constexpr std::string_view tmpl_infos = "\t- Computing per-atom entropy:\n"
                                                 "\t\t rcut  = %.5f\n"
                                                 "\t\t sigma = %.5f\n"
                                                 "\t\t local = %s\n";

public:

  inline void execute () override final 
  {
    assert(chunk_neighbors->number_of_cells() == grid->number_of_cells());
    bool has_chunk_neighbors = chunk_neighbors.has_value();
    if (!has_chunk_neighbors) {
      lerr << "No neighbors input data available" << std::endl;
      std::abort();
    }

    if (*verbose)
      lout << onika::format_string(std::string(tmpl_infos), *rcut, *sigma, *local ? "true": "false");

    LocalEntropyOp local_op{
        .verbose = *verbose,
        .rcut = *rcut,
        .sigma = *sigma,
        .nbins = *nbins,
        .local = *local,
    };

    local_op.initialize(*thermodynamic_state);

    auto field_id = field::mk_generic_real(*name);
    auto field_acc = grid->field_accessor(field_id);
    auto local_op_fields = make_field_tuple_from_field_set(compute_force_field_set, field_acc);

    ComputePairNullWeightIterator cp_weight{};
    ComputePairOptionalLocks<false> cp_locks{};
    exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{*chunk_neighbors};
    auto local_op_buf = make_compute_pair_buffer<ComputeBuffer>();
    ComputePairTrivialCellFiltering cpu_cell_filter = {};

    if (domain->xform_is_identity()) {
      NullXForm cp_xform;
      auto optional = make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter);
      compute_cell_particle_pairs(*grid, *rcut, *ghost, optional, local_op_buf, local_op, local_op_fields,
                                  parallel_execution_context());
    } else {
      LinearXForm cp_xform{domain->xform()};
      auto optional = make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter);
      compute_cell_particle_pairs(*grid, *rcut, *ghost, optional, local_op_buf, local_op, local_op_fields,
                                  parallel_execution_context());
    }

    if (*verbose)
      lout << onika::format_string("\t- Computing per-atom entropy END\n");
  }
};

template <class GridT> using ComputeLocalEntropyOpTmpl = ComputeLocalEntropyOp<GridT>;

// === register factories ===
ONIKA_AUTORUN_INIT(compute_local_entropy) {
  OperatorNodeFactory::instance()->register_factory("compute_local_entropy",
                                                    make_grid_variant_operator<ComputeLocalEntropyOpTmpl>);
}
} // namespace exaStamp
