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
#include <exanb/core/grid_additional_fields.h>
#include <exanb/core/domain.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/mechanical/compute_local_centrosymmetry.h>


struct MyGenericType {
  std::array<double, 32> m_data;
};

XNB_DECLARE_DYNAMIC_FIELD(MyGenericType, generic_arr3, "....");

namespace exaStamp {
using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;

template <typename ComputePairBufferT>
struct alignas(DEFAULT_ALIGNMENT) DummyOp {

  template<typename CellParticlesT>
  inline void operator()(size_t n, ComputePairBufferT& buf, MyGenericType& field, CellParticlesT cells) const {

  }

};
} // namespace exaStamp

namespace exaStamp
{
using namespace exanb;

template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz>>
class ComputeSteinhardtOp : public OperatorNode {

  ADD_SLOT(exanb::GridChunkNeighbors, chunk_neighbors, INPUT, exanb::GridChunkNeighbors{}, DocString{"neighbor list"});
  ADD_SLOT(bool, ghost, INPUT, false);
  ADD_SLOT(GridT, grid, INPUT);
  ADD_SLOT(Domain, domain, INPUT, REQUIRED);

  ADD_SLOT(double, rcut, INPUT, 0.0);
  ADD_SLOT(size_t, nnn, INPUT, 12);

  ADD_SLOT(std::string, name, INPUT, "orientorder");
  ADD_SLOT(bool, verbose, INPUT, false);

  using ComputeBuffer = ComputePairBuffer2<false, false>;
  using ComputeFields = FieldSet<>;
  static constexpr ComputeFields compute_field_set{};
  using CellsAccessorT = std::remove_cv_t<std::remove_reference_t<decltype(grid->cells_accessor())> >;

  static constexpr std::string_view tmpl_infos = "\t- Computing per-atom Steinhardt parameters:\n"
                                                 "\t\t rcut = %.5f\n"
                                                 "\t\t nnn  = %d\n";

public:
  inline void execute() override final {

    assert(chunk_neighbors->number_of_cells() == grid->number_of_cells());

    if (!chunk_neighbors.has_value()) {
      lerr << "No neighbors input data available" << std::endl;
      onika::fatal_error();
    }

    if (*verbose)
      lout << onika::format_string(std::string(tmpl_infos), *rcut, *nnn);

    if ((*nnn % 2) != 0) {
      lerr << "compute_local_centrosymmetry: the number of nearest neighbors (nnn = " << *nnn << ") must be even.\n";
      onika::fatal_error();
    }

    DummyOp<ComputeBuffer> local_op{};

    auto csp_field_id = field::mk_generic_arr3(*name);
    auto csp_field_acc = grid->field_accessor(csp_field_id);
    auto local_op_fields = make_field_tuple_from_field_set(FieldSet<>{}, csp_field_acc);

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
      lout << onika::format_string("\t- Computing per-atom centrosymmetry END\n");
  }
};

template <class GridT> using ComputeSteinhardtOpTmpl = ComputeSteinhardtOp<GridT>;

// === register factories ===
ONIKA_AUTORUN_INIT(compute_local_steinhardt) {
  OperatorNodeFactory::instance()->register_factory("compute_local_steinhardt",
                                                    make_grid_variant_operator<ComputeSteinhardtOpTmpl>);
}
} // namespace exaStamp
