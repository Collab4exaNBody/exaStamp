#include <string>

#include <mpi.h>

#include <onika/file_utils.h>
#include <onika/log.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/math/basic_types_yaml.h>
#include <onika/parallel/random.h>
#include <onika/print_utils.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <onika/string_utils.h>

// clang-format off
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/check_particles_inside_cell.h>
// clang-format on

#include <exaStamp/io/external_atom_format_v2.h>
#include <exaStamp/particle_species/particle_specie.h>

namespace exaStamp {

using namespace exanb;

struct PBC {
  bool x;
  bool y;
  bool z;
};

struct exaStampBox {

  Mat3d H = make_identity_matrix();
  Vec3d origin{};

  Mat3d inv_H{};
  Mat3d D1{};
  Mat3d D2{};
  Mat3d invD1{};
  Mat3d invD2{};

  Mat3d Xform = make_identity_matrix();
  Mat3d invXform = make_identity_matrix();
  Vec3d boxlen{};

  bool uniform_scale;

  void init_from_context(aio::io::IOContext& ctx, Domain& domain, ReadBoundsSelectionMode bounds_mode) {

    H = ctx.data.cell;
    origin = ctx.data.origin;

    Vec3d a{H.m11, H.m21, H.m31};
    Vec3d b{H.m12, H.m22, H.m32};
    Vec3d c{H.m13, H.m23, H.m33};

    boxlen.x = norm(a);
    boxlen.y = norm(b);
    boxlen.z = norm(c);

    lout << "a = " << a << onika::format_string(" , norm = %10.5f ", boxlen.x) << std::endl;
    lout << "b = " << b << onika::format_string(" , norm = %10.5f ", boxlen.y) << std::endl;
    lout << "c = " << c << onika::format_string(" , norm = %10.5f ", boxlen.z) << std::endl;

    D1 = diag_matrix(boxlen);
    invD1 = inverse(D1);

    // check for grid dims. If grid dims is not provided by the user,
    // it will be deducde from box lengths and cellsize.
    double cellsize = domain.cell_size();
    IJK dim = domain.grid_dimension();
    size_t nx = static_cast<ssize_t>(boxlen.x / cellsize);
    size_t ny = static_cast<ssize_t>(boxlen.y / cellsize);
    size_t nz = static_cast<ssize_t>(boxlen.z / cellsize);
    nx = (dim.i > 0 && cmp::ne(dim.i, nx)) ? dim.i : nx;
    ny = (dim.j > 0 && cmp::ne(dim.j, ny)) ? dim.j : ny;
    nz = (dim.k > 0 && cmp::ne(dim.k, nz)) ? dim.k : nz;

    D2 = diag_matrix(Vec3d{nx * cellsize, ny * cellsize, nz * cellsize});
    invD2 = inverse(D2);

    Mat3d F1 = H * invD1;
    Mat3d F2 = D1 * invD2;

    Xform = F1 * F2;
    invXform = inverse(Xform);

    AABB bounds{{0., 0., 0.}, {D2.m11, D2.m22, D2.m33}};
    compute_domain_bounds(domain, bounds_mode, 0.0, bounds, bounds, true);
    uniform_scale = onika::math::is_uniform_scale(Xform);
    domain.set_xform(Xform);
  }

  inline void wrap(Vec3d& r) {
    r.x = cmp::wrap(r.x, D2.m11);
    r.y = cmp::wrap(r.y, D2.m22);
    r.z = cmp::wrap(r.z, D2.m33);
  }
};

template <class GridT, class = AssertGridHasFields<GridT, field::_rx, field::_ry, field::_rz, field::_vx, field::_vy,
                                                   field::_vz, field::_id, field::_type>>
class ReadExternalAtomFormatV2Node : public OperatorNode {
  ADD_SLOT(MPI_Comm, mpi, INPUT, MPI_COMM_WORLD);
  ADD_SLOT(Domain, domain, INPUT_OUTPUT);
  ADD_SLOT(GridT, grid, INPUT_OUTPUT);

  ADD_SLOT(ParticleSpecies, species, INPUT);
  ADD_SLOT(ReadBoundsSelectionMode, bounds_mode, INPUT, ReadBoundsSelectionMode::COMPUTED_BOUNDS);

  ADD_SLOT(std::string, file, INPUT, REQUIRED);
  ADD_SLOT(std::string, format, INPUT, "");
  ADD_SLOT(std::string, compression, INPUT, "");
  ADD_SLOT(PBC, pbc, INPUT, PBC{true, true, true});

  ADD_SLOT(bool, remap_atom, INPUT, false);

  // ADD_SLOT(std::string, units       , INPUT, "metal" );
  // ADD_SLOT(std::string, style       , INPUT, "atomic");

public:
  inline void execute() override final {

    lout << "======== Read external file format v2 ========\n";

    int rank = 0, np = 1;
    MPI_Comm_rank(*mpi, &rank);
    MPI_Comm_size(*mpi, &np);

    Domain& domain = *(this->domain);
    GridT& grid = *(this->grid);
    ParticleSpecies species = *(this->species);
    PBC& pbc = *(this->pbc);

    domain.set_expandable(false);
    domain.set_periodic_boundary(pbc.x, pbc.y, pbc.z);
    domain.set_xform(make_identity_matrix());

    // if the number of particule > 0, then the system has already be defined
    if (!grid.number_of_particles() == 0) {
      lerr << "System already defined. Skipping.";
      return;
    }

    exaStampBox box;
    aio::io::IOContext io_ctx{};
    io_ctx.set<aio::io::IOContext::REMAP_ATOM>(*remap_atom);

    if (this->species.has_value()) {
      for (size_t i = 0; i < species.size(); ++i) {
        io_ctx.species.add_species(species.at(i).m_name, i);
      }
    }

    // only read on rank 0
    if (rank == 0) {
      if (!read_external_atom_file(io_ctx, *file, *format, *compression)) {
        PANIC("Fail to read.... PANIC!")
      }

      linfo("Parsing done in %.0f ms | %.0f µs | %.0f ns\n", io_ctx.timer.x, io_ctx.timer.y, io_ctx.timer.z);

      box.init_from_context(io_ctx, domain, *bounds_mode);

      lout << std::endl;
      lout << "Uniform scale    = " << std::boolalpha << box.uniform_scale << std::endl;
      lout << "Particles        = " << io_ctx.data.particles.size() << std::endl;
      lout << "Domain XForm     = " << domain.xform() << std::endl;
      lout << "Domain bounds    = " << domain.bounds() << std::endl;
      lout << "Domain size      = " << bounds_size(domain.bounds()) << std::endl;
      lout << "Real size        = "
           << bounds_size(domain.bounds()) * Vec3d{domain.xform().m11, domain.xform().m22, domain.xform().m33}
           << std::endl;
      lout << "Cell size        = " << domain.cell_size() << std::endl;
      lout << "Grid dimensions  = " << domain.grid_dimension() << " (" << grid_cell_count(domain.grid_dimension())
           << " cells)" << std::endl;
    }

    // Send bounds and box_size to all cores
    MPI_Bcast(&domain, sizeof(Domain), MPI_CHARACTER, 0, *mpi);
    grid.set_origin(domain.bounds().bmin);
    grid.set_offset(IJK{0, 0, 0});
    grid.set_cell_size(domain.cell_size());
    grid.set_dimension(domain.grid_dimension());

    // Assign particles on the grid
    if (rank == 0) {
      auto cells = grid.cells();

      for (size_t i = 0; i < io_ctx.data.particles.size(); ++i) {
        aio::policy::ExaStampTraits::ParticleTupleIO& p = io_ctx.data.particles[i];
        Vec3d r{p[field::rx], p[field::ry], p[field::rz]};

        r = box.invXform * r;
        box.wrap(r);

        p[field::rx] = r.x;
        p[field::ry] = r.y;
        p[field::rz] = r.z;

        IJK loc = domain_periodic_location(domain, r);
        size_t cell_index = grid_ijk_to_index(domain.grid_dimension(), loc);
        cells[cell_index].push_back(p);
      }
    }

    // collective built of particle offset
    grid.rebuild_particle_offsets();

    // Sanity check on the grid
    if (!check_domain(domain)) {
      lout << "domain = " << domain << std::endl;
      fatal_error() << "Invalid domain configuration" << std::endl;
    }

    lout << "===========================================\n\n";
  }
};

template <class GridT> using ReadExternalAtomFormatNodeV2Tmpl = ReadExternalAtomFormatV2Node<GridT>;

// === register factories ===
ONIKA_AUTORUN_INIT(read_external_atom_format) {
  OperatorNodeFactory::instance()->register_factory("read_external_atom_format_v2",
                                                    make_grid_variant_operator<ReadExternalAtomFormatNodeV2Tmpl>);
}

} // namespace exaStamp
