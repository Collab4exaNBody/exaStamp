#include <mpi.h>

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
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/mechanical/compute_local_entropy.h>
#include <exaStamp/particle_species/particle_specie.h>

namespace exaStamp {
using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;


struct alignas(DEFAULT_ALIGNMENT) DummyOp {

  size_t nbin;
  size_t nt;
  double dr, drinv;
  std::vector<double>& buffer;

  template <class CellParticlesT>
  inline void operator()(size_t n, ComputePairBuffer2<false, false>& buf, CellParticlesT) const {

    size_t tid = omp_get_thread_num();
    double* __restrict local_thread_buf = buffer.data() + tid * nbin;

    for (size_t i = 0; i < n; ++i) {
      double r = sqrt(buf.d2[i]);
      size_t ibin = static_cast<size_t>(r * drinv);
      if (ibin >= nbin)
        continue;
      local_thread_buf[ibin] += 1.0;
    }

    // for (size_t i = 0; i < nbin; ++i) {
    //   local_thread_buf[i] += 1.0;
    // }

    // double* __restrict local_buf = buffer.data() + tid * ncells * nbin;
    // double* __restrict local_buf = buffer.data() + ncells * nt * nbin + tid;

    // lout << tid << std::endl;
    // lock[buf.cell].lock();
    // for (size_t i = 0; i < 10; ++i) {
    //   buffer[buf.cell][i] += 1.0;
    // }
    // lock[buf.cell].unlock();
  }
};
} // namespace exaStamp

namespace exaStamp
{

using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;

template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz>>
class ComputeRadialDistributionFunctionOperator : public OperatorNode {

  // ========= I/O slots =======================
  ADD_SLOT( MPI_Comm              , mpi                , INPUT);
  ADD_SLOT( exanb::GridChunkNeighbors    , chunk_neighbors   , INPUT , exanb::GridChunkNeighbors{} , DocString{"neighbor list"} );
  ADD_SLOT( bool                  , ghost              , INPUT , false );
  ADD_SLOT( GridT                 , grid               , INPUT );
  ADD_SLOT( Domain                , domain             , INPUT , REQUIRED );
  ADD_SLOT( ThermodynamicState    , thermodynamic_state, INPUT, REQUIRED );

  ADD_SLOT(double, rcut , INPUT, 0.0);
  ADD_SLOT(size_t, nbins,  INPUT, 10);

  using ComputeBuffer = ComputePairBuffer2<false, false>;
  using ComputeFields = FieldSet<>;
  static constexpr ComputeFields compute_force_field_set{};

public:
  inline void execute() override final {

    size_t nbin = *nbins;
    int rank;
    MPI_Comm_rank(*mpi, &rank);

    assert(chunk_neighbors->number_of_cells() == grid->number_of_cells());
    bool has_chunk_neighbors = chunk_neighbors.has_value();
    if (!has_chunk_neighbors) {
      lerr << "No neighbors input data available" << std::endl;
      std::abort();
    }

    lout << onika::format_string("\t- Computing radial distribution function\n");

    size_t num_thread = omp_get_max_threads();

    double dr = *rcut / nbin;
    double drinv = 1.0 / dr;

    std::vector<double> hist(nbin);
    std::vector<double> per_rank_hist(nbin);
    std::vector<double> per_thread_hist(nbin * num_thread);
    DummyOp local_op{nbin, num_thread, dr, drinv, per_thread_hist};

    auto local_op_fields = make_field_tuple_from_field_set(FieldSet<>{});

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

    double volume = thermodynamic_state->volume();
    size_t particle_count = thermodynamic_state->particle_count();
    double density = static_cast<double>(particle_count) / volume;
    // std::vector<double> result(nbin);

    for (size_t n = 0; n < num_thread; ++n) {
      double* __restrict local_buf = per_thread_hist.data() + n * nbin;
      for (size_t i = 0; i < nbin; ++i) {
        per_rank_hist[i] += local_buf[i];
      }
    }

    MPI_Reduce(&per_rank_hist[0], &hist[0], nbin, MPI_DOUBLE, MPI_SUM, 0, *mpi);


    for (size_t i = 0; i < nbin; ++i) {
      double rlo = i * dr;
      double rhi = rlo + dr;
      double shell_volume = (4. / 3.) * M_PI * (rhi*rhi*rhi - rlo*rlo*rlo);
      double norm = 1.0 / (density * particle_count * shell_volume);
      hist[i] *= norm;
    }

    for (size_t i = 0; i < nbin; ++i) {
      lout << i << " " << hist[i] << std::endl;
    }

    // for (size_t c = 0; c < n_cells; c++) {
    //   if (!grid->is_ghost_cell(c)) {
    //     lout << "\t " << c << " " << cells[c].size() << std::endl;
    //     for (size_t i = 0; i < nbin; i++)  {
    //       lout << "\t\t " << test[c][i] << std::endl;
    //     }
    //   }
    // }

    lout << onika::format_string("\t- Computing radial distribution function END\n");
  }
};

template <class GridT>
using ComputeRadialDistributionFunctionOperatorTmpl = ComputeRadialDistributionFunctionOperator<GridT>;

// === register factories ===
ONIKA_AUTORUN_INIT(compute_local_entropy) {
  OperatorNodeFactory::instance()->register_factory(
      "compute_rdf", make_grid_variant_operator<ComputeRadialDistributionFunctionOperatorTmpl>);
}
} // namespace exaStamp
