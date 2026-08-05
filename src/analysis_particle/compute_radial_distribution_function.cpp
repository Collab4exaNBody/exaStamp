#include <mpi.h>
#include <filesystem>

#include <onika/cpp_utils.h>
#include <onika/file_utils.h>
#include <onika/print_utils.h>
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
#include <exanb/core/particle_type_properties.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exaStamp/mechanical/cmp_utils.h>

// TODO: 
// - File size might become an issue
// - Check for symmetric re-normalization

namespace exaStamp {
using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;

template <bool Symmetric, typename ComputePairBufferT, typename CellParticlesT>
struct alignas(DEFAULT_ALIGNMENT) PairDistFunctor {

  struct alignas(onika::memory::DEFAULT_ALIGNMENT) NeighEntry {
    size_t bin, pair_index;
  };

  const double m_rcut_sq{};
  const size_t m_nbin{};
  const size_t m_ntype{};
  const double m_inv_dr{};
  std::vector<double>& m_buffer;
  std::vector<int>& m_type_counter;

  inline void operator()(size_t n, ComputePairBufferT& buf, size_t type_i, CellParticlesT cells) const {

    size_t thread_id = omp_get_thread_num();
    size_t num_thread = omp_get_num_threads();
    std::array<NeighEntry, exanb::MAX_PARTICLE_NEIGHBORS> nbh_entries{};

    m_type_counter[type_i * num_thread + thread_id] += 1;

    double* __restrict data = m_buffer.data();

    size_t nvalid = 0;
    size_t cell_j, index_j;
    for (size_t j = 0; j < n; ++j) {

      if (buf.d2[j] > m_rcut_sq)
        continue;

      buf.nbh.get(j, cell_j, index_j);
      size_t type_j = cells[cell_j][field::type][index_j];

      double rij = sqrt(buf.d2[j]);
      size_t ibin = static_cast<size_t>(rij * m_inv_dr);

      if constexpr (Symmetric) {
        size_t ti = cmp::min(type_i, type_j);
        size_t tj = cmp::max(type_i, type_j);
        nbh_entries[nvalid++] = {ibin, cmp::upper_triangle_index(ti, tj)};
      } else {
        nbh_entries[nvalid++] = {ibin, type_i * m_ntype + type_j};
      }

    }

    std::sort(nbh_entries.begin(), nbh_entries.begin() + nvalid, [](const NeighEntry& a, const NeighEntry& b) {
      return (a.pair_index < b.pair_index) || (a.pair_index == b.pair_index && a.bin < b.bin);
    });

    for (size_t j = 0; j < nvalid; ++j) {
      const NeighEntry& nbh = nbh_entries[j];
      size_t offset = (nbh.pair_index * num_thread + thread_id) * m_nbin + nbh.bin;
      data[offset] += 1.0;
    }
  }
};
} // namespace exaStamp

namespace exaStamp {

using namespace exanb;

template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz, field::_type>>
class ComputeRadialDistributionFunctionOperator : public OperatorNode {

  // ========= I/O slots =======================
  ADD_SLOT(MPI_Comm, mpi, INPUT);
  ADD_SLOT(exanb::GridChunkNeighbors, chunk_neighbors, INPUT, exanb::GridChunkNeighbors{}, DocString{"neighbor list"});
  ADD_SLOT(bool, ghost, INPUT, false);
  ADD_SLOT(GridT, grid, INPUT);
  ADD_SLOT(double, rcut_max, INPUT_OUTPUT, 0.0);
  ADD_SLOT(Domain, domain, INPUT, REQUIRED);

  ADD_SLOT(ThermodynamicState, thermodynamic_state, INPUT, REQUIRED);
  ADD_SLOT(ParticleSpecies, species, INPUT, REQUIRED);
  ADD_SLOT(long, timestep, INPUT, REQUIRED);
  ADD_SLOT(double, physical_time, INPUT);

  ADD_SLOT(double, rcut, INPUT, 5.0);
  ADD_SLOT(size_t, nbin, INPUT, 100);
  ADD_SLOT(bool, symmetric, INPUT, false); // symmetrize pair interaction
  ADD_SLOT(std::string, dirname, INPUT, "dir.rdf");
  ADD_SLOT(bool, verbose, INPUT, false);

  static constexpr bool UseWeights = true;
  static constexpr bool UseNeighbours = true;
  using NeighFieldSet = FieldSet<field::_type>;
  using ComputeFieldsSet = FieldSet<field::_type>;

  using ComputeBuffer = ComputePairBuffer2<
    UseWeights,
    UseNeighbours,
    NoExtraStorage,
    DefaultComputePairBufferAppendFunc,
    exanb::MAX_PARTICLE_NEIGHBORS,
    ComputePairBuffer2Weights,
    NeighFieldSet
  >;

  using CellsAccessorT = std::remove_cv_t<std::remove_reference_t<decltype(grid->cells_accessor())> >;
  static constexpr std::true_type use_cells_accessor{};

  template<bool Symmetric> using RadialDistFnOp = PairDistFunctor<Symmetric, ComputeBuffer, CellsAccessorT>;

public:
  inline void execute() override final {

    int rank;
    MPI_Comm_rank(*mpi, &rank);
    size_t num_thread = omp_get_max_threads();

    assert(chunk_neighbors->number_of_cells() == grid->number_of_cells());
    bool has_chunk_neighbors = chunk_neighbors.has_value();
    if (!has_chunk_neighbors) {
      lerr << "No neighbors input data available" << std::endl;
      std::abort();
    }

    if (*verbose)
      lout << onika::format_string("\t- Computing radial distribution function\n");

    size_t ntype = static_cast<size_t>((*species).size());
    size_t npair = (*symmetric) ? static_cast<size_t>(ntype * (ntype + 1) / 2) : static_cast<size_t>(ntype * ntype);

    // init file output
    std::filesystem::path dirpath = std::filesystem::absolute(*dirname);
    if ((rank == 0) && (!std::filesystem::exists(dirpath))) {
      std::filesystem::create_directories(dirpath);
    }

    std::vector<std::string> filearray(npair + 1);
    filearray[npair] = onika::format_string("%s/rdf-total.dat", dirpath.c_str());

    if (*symmetric) {
      for (size_t i = 0; i < ntype; i++) {
        for (size_t j = i; j < ntype; j++) {
          filearray[cmp::upper_triangle_index(i, j)] = onika::format_string(
              "%s/rdf-%s-%s.dat", dirpath.c_str(), (*species).at(i).m_name, (*species).at(j).m_name);
        }
      }
    } else {
      for (size_t i = 0; i < ntype; ++i) {
        for (size_t j = 0; j < ntype; ++j) {
          filearray[i * ntype + j] = onika::format_string(
              "%s/rdf-%s-%s.dat", dirpath.c_str(), (*species).at(i).m_name, (*species).at(j).m_name);
        }
      }
    }

    *rcut_max = std::max(*rcut, *rcut_max);
    double dr = *rcut / static_cast<double>(*nbin);
    double inv_dr = 1.0 / dr;

    std::vector<double> hist((npair + 1) * (*nbin));
    std::vector<double> per_rank_buffer((npair + 1) * (*nbin));
    std::vector<double> per_thread_buffer((npair + 1) * num_thread * (*nbin));

    std::vector<int> type_counter(ntype);
    std::vector<int> rank_type_counter(ntype);
    std::vector<int> thread_type_counter(ntype * num_thread);

    auto local_op_buf = make_compute_pair_buffer<ComputeBuffer>();
    auto local_op_fields = make_field_tuple_from_field_set(ComputeFieldsSet{});

    ComputePairNullWeightIterator cp_weight{};
    ComputePairOptionalLocks<false> cp_locks{};
    exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{*chunk_neighbors};
    auto cell_acc = grid->field_accessors_from_field_set(NeighFieldSet{});
    ComputePairTrivialCellFiltering cell_filter{};
    ComputePairTrivialParticleFiltering partitcle_filter{};

    LinearXForm cp_xform{domain->xform()};
    auto optional =
        make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks, cell_filter, partitcle_filter, cell_acc);

    if (*symmetric) {
      RadialDistFnOp<true> local_op{(*rcut) * (*rcut), *nbin, ntype, inv_dr, per_thread_buffer, thread_type_counter};
      compute_cell_particle_pairs2(*grid, *rcut, false, optional, local_op_buf, local_op, local_op_fields,
                                   DefaultPositionFields{}, parallel_execution_context(), use_cells_accessor);
    } else {
      RadialDistFnOp<false> local_op{(*rcut) * (*rcut), *nbin, ntype, inv_dr, per_thread_buffer, thread_type_counter};
      compute_cell_particle_pairs2(*grid, *rcut, false, optional, local_op_buf, local_op, local_op_fields,
                                   DefaultPositionFields{}, parallel_execution_context(), use_cells_accessor);
    }

    // reduce per-thread data
    double* __restrict local_total_buf = per_rank_buffer.data() + (npair * (*nbin));
    for (size_t i = 0; i < npair; ++i) {
      for (size_t j = 0; j < num_thread; ++j) {
        double* __restrict local_thread_buf = per_thread_buffer.data() + (i * num_thread + j) * (*nbin);
        double* __restrict local_rank_buf = per_rank_buffer.data() + (i * (*nbin));
        #pragma omp simd
        for (size_t k = 0; k < (*nbin); ++k) {
          local_rank_buf[k] += local_thread_buf[k];
          local_total_buf[k] += local_thread_buf[k];
        }
      }
    }

    for (size_t i = 0; i < ntype; ++i) {
      for (size_t j = 0; j < num_thread; ++j) {
        rank_type_counter[i] += thread_type_counter[i * num_thread + j];
      }
    }

    // reduce per-rank data
    MPI_Reduce(&per_rank_buffer[0], &hist[0], (npair + 1) * (*nbin), MPI_DOUBLE, MPI_SUM, 0, *mpi);
    MPI_Reduce(&rank_type_counter[0], &type_counter[0], ntype, MPI_INT, MPI_SUM, 0, *mpi);

    if (rank != 0)
      return;

    // some normalization
    
    double volume = thermodynamic_state->volume();
    size_t particle_count = thermodynamic_state->particle_count();
    double density = static_cast<double>(particle_count) / volume;

    for (size_t n = 0; n < npair; ++n) {

      size_t i, j;
      if (*symmetric) {
        cmp::unflatten_lower(n, i, j);
      } else {
        i = n / ntype;
        j = n % ntype;
      }

      double* __restrict buf = hist.data() + (n * (*nbin));
      size_t Na = type_counter[i];
      size_t Nb = type_counter[j];

      double norm;
      if ((Na == 0) || (Nb == 0)) {
        norm = 0.0;
      } else {
        norm = 1.0 / ((Nb / volume) * Na);
      }

      for (size_t k = 0; k < *nbin; ++k) {
        double rlo = k * dr;
        double rhi = rlo + dr;
        double shell_volume = (4. / 3.) * M_PI * (rhi * rhi * rhi - rlo * rlo * rlo);
        buf[k] *= norm / shell_volume;
      }
    }

    double* __restrict buf = hist.data() + (npair * (*nbin));
    for (size_t k = 0; k < *nbin; ++k) {
      double rlo = k * dr;
      double rhi = rlo + dr;
      double shell_volume = (4. / 3.) * M_PI * (rhi * rhi * rhi - rlo * rlo * rlo);
      double norm = 1.0 / (density * particle_count * shell_volume);
      buf[k] *= norm;
    }

    // write results to files
    for (size_t i = 0; i < (npair + 1); ++i) {

      std::string line;
      line.reserve((9 + 18 + 5 + 18) + (*nbin) * 18 + 1);

      char tmp[32];

      snprintf(tmp, sizeof(tmp), "%9ld", *timestep);
      line.append(tmp);
      snprintf(tmp, sizeof(tmp), "%18.6e", *physical_time);
      line.append(tmp);
      snprintf(tmp, sizeof(tmp), "%5ld", *nbin);
      line.append(tmp);
      snprintf(tmp, sizeof(tmp), "%18.6e", *rcut);
      line.append(tmp);

      double* __restrict buf = hist.data() + (i * (*nbin));
      for (size_t k = 0; k < (*nbin); ++k) {
        snprintf(tmp, sizeof(tmp), "%18.6e", buf[k]);
        line.append(tmp);
      }
      line.append("\n");
      onika::FileAppendWriteBuffer::instance().append_to_file(filearray[i], line, false);
    }

    if (*verbose)
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
