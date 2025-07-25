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
#include <exanb/core/particle_type_properties.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>

#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/particle_species/particle_specie.h>

#include <cxxabi.h>

namespace exaStamp {
using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;

template<typename ComputePairBufferT, typename CellParticlesT>
struct alignas(DEFAULT_ALIGNMENT) PairDistFunctor {

  const double m_rcut_sq{};
  const size_t m_nbin{};
  const double m_dr{};
  const double m_inv_dr{};

  std::vector<double>& m_buffer;

  inline void operator()(size_t n, ComputePairBufferT& buf, int type_a, CellParticlesT cells) const {

    size_t cell, index;
    buf.nbh.get(0, cell, index);

    int v1 = cells[cell][field::type][index];
    int v2 = buf.nbh_pt[0][field::type];

    std::string name = typeid(cells[cell][field::type][index]).name();
    const char* name2 = abi::__cxa_demangle(name.c_str(), 0, 0, 0);
    lout << onika::format_string("> type = (v1: %d, v2: %d) %s\n", v1, v2, name2);

    // for (size_t i = 0; i < n; ++i) {
    //   lout << type_a << " " << "> " << buf.nbh_pt[i][field::type] << " <" << std::endl;
    // }

    // const int type_b = buf.nbh_pt[0][field::type];
    // lout << "> " << type_b << " > " << buf.nbh_pt[0][field::fx] << " <" << std::endl;

    // lout << buf.nbh[0][field::_type] << std::endl;

    // size_t thread_id = omp_get_thread_num();
    // double* __restrict local_thread_buf = m_buffer.data() + thread_id * m_nbin;

    // for (size_t i = 0; i < n; ++i) {

    //   if (buf.d2[i] <= m_rcut_sq) {
    //     double r = sqrt(buf.d2[i]);
    //     size_t j = static_cast<size_t>(r, m_inv_dr);
    //     if (j >= m_nbin) {

    //     }
    //   }


    //   double r = sqrt(buf.d2[i]);
    //   size_t ibin = static_cast<size_t>(r * drinv);
    //   if (ibin >= nbin)
    //     continue;
    //   local_thread_buf[ibin] += 1.0;
    // }

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

namespace exaStamp {

using namespace exanb;

using onika::memory::DEFAULT_ALIGNMENT;

template <class GridT, class = AssertGridHasFields<GridT, field::_ep, field::_fx, field::_fy, field::_fz, field::_type>>
class ComputeRadialDistributionFunctionOperator : public OperatorNode {


  // ========= I/O slots =======================
  ADD_SLOT(MPI_Comm, mpi, INPUT);
  ADD_SLOT(exanb::GridChunkNeighbors, chunk_neighbors, INPUT, exanb::GridChunkNeighbors{}, DocString{"neighbor list"});
  ADD_SLOT(bool, ghost, INPUT, false);
  ADD_SLOT(GridT, grid, INPUT);
  ADD_SLOT(double, rcut_max, INPUT_OUTPUT, 0.0, DocString{"Updated max rcut"});
  ADD_SLOT(Domain, domain, INPUT, REQUIRED);
  ADD_SLOT(ThermodynamicState, thermodynamic_state, INPUT, REQUIRED);

  ADD_SLOT(ParticleTypeProperties, particle_type_properties, INPUT, ParticleTypeProperties{});
  ADD_SLOT(double                , rcut                    , INPUT, 5.0                     );
  ADD_SLOT(size_t                , nbin                    , INPUT, 100                     );


  // using ComputeBuffer = ComputePairBuffer2<false, true>;
  
  static constexpr bool UseWeights = false;
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

public:
  inline void execute() override final {

    int rank;
    MPI_Comm_rank(*mpi, &rank);
    size_t num_thread = omp_get_max_threads();

    *rcut_max = std::max(*rcut, *rcut_max);

    assert(chunk_neighbors->number_of_cells() == grid->number_of_cells());
    bool has_chunk_neighbors = chunk_neighbors.has_value();
    if (!has_chunk_neighbors) {
      lerr << "No neighbors input data available" << std::endl;
      std::abort();
    }

    lout << onika::format_string("\t- Computing radial distribution function\n");

    double dr = *rcut / static_cast<double>(*nbin);
    double inv_dr = 1.0 / dr;

    std::vector<double> hist(*nbin);
    std::vector<double> per_rank_hist(*nbin);
    std::vector<double> per_thread_hist(*nbin * num_thread);
    
    PairDistFunctor<ComputeBuffer, CellsAccessorT> local_op{(*rcut)*(*rcut), *nbin, dr, inv_dr, per_thread_hist};

    auto local_op_fields = make_field_tuple_from_field_set(ComputeFieldsSet{});

    ComputePairNullWeightIterator cp_weight{};
    ComputePairOptionalLocks<false> cp_locks{};
    exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{*chunk_neighbors};
    auto local_op_buf = make_compute_pair_buffer<ComputeBuffer>();

    // using has_field_type_t = typename GridT:: template HasField<field::_type>;
    // constexpr bool has_field_type = has_field_type_t::value;
    // lout << "okok: " << has_field_type << std::endl;
    // if (particle_type_properties.has_value()) {
    //   for (const auto & it : particle_type_properties->m_scalars) {
    //     lout << " >> " << it.first << std::endl;
    //   }
    // }

    LinearXForm cp_xform{domain->xform()};
    auto optional =
        make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks, ComputePairTrivialCellFiltering{},
                                        ComputePairTrivialParticleFiltering{});

    compute_cell_particle_pairs2(*grid, *rcut, false, optional, local_op_buf, local_op, local_op_fields,
                                 DefaultPositionFields{}, parallel_execution_context(), use_cells_accessor);

    // using OkFields = onika::FlatTuple< onika::soatl::FieldId<field::_rx> , onika::soatl::FieldId<field::_ry> , onika::soatl::FieldId<field::_rz>, onika::soatl::FieldId<field::_type> >;
    // static constexpr OkFields posfields;
    
    // if (domain->xform_is_identity()) {
    //   NullXForm cp_xform;
    //   // auto optional = make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks, cpu_cell_filter);
    //   auto optional = make_compute_pair_optional_args(
    //       nbh_it, cp_weight, cp_xform, cp_locks, ComputePairTrivialCellFiltering{},
    //       ComputePairTrivialParticleFiltering{}, grid->field_accessors_from_field_set(FieldSet<field::_type>{}));
    //   // compute_cell_particle_pairs(*grid, *rcut, *ghost, optional, local_op_buf, local_op, local_op_fields,
    //                               // parallel_execution_context());
    //   compute_cell_particle_pairs2(*grid, *rcut, *ghost, optional, local_op_buf, local_op, local_op_fields, DefaultPositionFields{}, parallel_execution_context(), use_cells_accessor);
    // } else {
    //   LinearXForm cp_xform{domain->xform()};
    //   // auto optional = make_compute_pair_optional_args(nbh_it, cp_weight, cp_xform, cp_locks);
    //   // compute_cell_particle_pairs(*grid, *rcut, *ghost, optional, local_op_buf, local_op, local_op_fields,
    //   //                             parallel_execution_context());
    //   // compute_cell_particle_pairs2(*grid, *rcut, *ghost, optional, local_op_buf, local_op, local_op_fields, posfields, parallel_execution_context(), use_cells_accessor);
    //   auto optional = make_compute_pair_optional_args(
    //       nbh_it, cp_weight, cp_xform, cp_locks, ComputePairTrivialCellFiltering{},
    //       ComputePairTrivialParticleFiltering{}, grid->field_accessors_from_field_set(FieldSet<field::_type>{}));
    //   // compute_cell_particle_pairs(*grid, *rcut, *ghost, optional, local_op_buf, local_op, local_op_fields,
    //                               // parallel_execution_context());
    //   compute_cell_particle_pairs2(*grid, *rcut, *ghost, optional, local_op_buf, local_op, local_op_fields, DefaultPositionFields{}, parallel_execution_context(), use_cells_accessor);
    // }

    // double volume = thermodynamic_state->volume();
    // size_t particle_count = thermodynamic_state->particle_count();
    // double density = static_cast<double>(particle_count) / volume;
    // std::vector<double> result(nbin);

    // for (size_t n = 0; n < num_thread; ++n) {
    //   double* __restrict local_buf = per_thread_hist.data() + n * nbin;
    //   for (size_t i = 0; i < nbin; ++i) {
    //     per_rank_hist[i] += local_buf[i];
    //   }
    // }

    // MPI_Reduce(&per_rank_hist[0], &hist[0], nbin, MPI_DOUBLE, MPI_SUM, 0, *mpi);


    // for (size_t i = 0; i < nbin; ++i) {
    //   double rlo = i * dr;
    //   double rhi = rlo + dr;
    //   double shell_volume = (4. / 3.) * M_PI * (rhi*rhi*rhi - rlo*rlo*rlo);
    //   double norm = 1.0 / (density * particle_count * shell_volume);
    //   hist[i] *= norm;
    // }

    // for (size_t i = 0; i < nbin; ++i) {
    //   lout << i << " " << hist[i] << std::endl;
    // }

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

// template <class GridT>
// using ComputeRadialDistributionFunctionOperatorTmpl = ComputeRadialDistributionFunctionOperator<GridT>;

// === register factories ===
ONIKA_AUTORUN_INIT(compute_local_entropy) {
  OperatorNodeFactory::instance()->register_factory(
      "compute_rdf", make_grid_variant_operator<ComputeRadialDistributionFunctionOperator>);
}
} // namespace exaStamp
