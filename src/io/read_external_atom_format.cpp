#include <ctime>
#include <filesystem>
#include <string>

#include <mpi.h>

#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/file_utils.h>
#include <onika/print_utils.h>
#include <onika/string_utils.h>
#include <onika/parallel/random.h>
#include <onika/log.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/check_particles_inside_cell.h>

#include <exaStamp/particle_species/particle_specie.h>

#include <exaStamp/io/external_atom_format.h>

namespace exaStamp {

using namespace exanb;
using namespace ioextra;

template<
  class GridT,
  class = AssertGridHasFields< GridT, field::_rx, field::_ry, field::_rz, field::_vx, field::_vy, field::_vz, field::_id, field::_type >
>
class ReadExternalAtomFormatNode : public OperatorNode
{
  ADD_SLOT( MPI_Comm        , mpi          , INPUT , MPI_COMM_WORLD  );
  ADD_SLOT( Domain          , domain       , INPUT_OUTPUT );
  ADD_SLOT( GridT           , grid         , INPUT_OUTPUT );

  // optional. if no species given, type ids are allocated automatically
  ADD_SLOT( ParticleSpecies , species      , INPUT ); 

  ADD_SLOT( ReadBoundsSelectionMode, bounds_mode     , INPUT , ReadBoundsSelectionMode::FILE_BOUNDS );

  ADD_SLOT(std::string, file        , INPUT, REQUIRED );
  ADD_SLOT(std::string, format      , INPUT, "");
  ADD_SLOT(std::string, compression , INPUT, "");
  // ADD_SLOT(std::string, units       , INPUT, "metal" );
  // ADD_SLOT(std::string, style       , INPUT, "atomic");

  // ADD_SLOT(size_t     , step        , INPUT, 0);

  // ADD_SLOT(double     , noise        , INPUT , 0.0 );
  // ADD_SLOT(double     , noise_cutoff , INPUT , OPTIONAL );

public:

  inline void execute () override final {

    lout << "======== Read external file format ========" << std::endl;

    // MPI initialization
    int rank = 0, np = 1;
    MPI_Comm_rank(*mpi, &rank);
    MPI_Comm_size(*mpi, &np);

    Domain& domain = *(this->domain);
    GridT& grid = *(this->grid);
    ParticleSpecies species = *(this->species);
    ParticleData particle_data;

    // if the number of particule > 0, then the system has already be defined
    if (!grid.number_of_particles() == 0) {
      lerr << "System already defined. Skipping.";
      return;
    }

    // Ensure the file exists
    std::string filepath = std::string(std::filesystem::absolute( *file ));
    if (! std::filesystem::exists( filepath )) {
      lerr << onika::format_string("Input file doest not exists: '%s'\n", filepath.c_str());
      std::abort();
    }

    IOContext ctx{.particle_data = particle_data};

    if (this->species.has_value()) {
      for (size_t i = 0; i < species.size(); i++) {
        ctx.species.push_back(species.at(i).m_name);
      }
    }

    if ( rank == 0 ) {

      // create reader object
      FileInfos infos = get_file_infos(filepath, *format, *compression);
      std::unique_ptr<Parser> parser = infos.format.creator(filepath, FileMode::READ, infos.compression);

      if (parser == nullptr) {
        lerr << "Something went wrong... I AM PANICKING" << std::endl;
        std::abort();
      }

      lout << std::endl;

      // need to dereference the ptr to use call operator
      if (!(*parser)(ctx)) {
        lerr << "Something went wrong when parsing the file" << std::endl;
        // lerr << parser->current_line() << std::endl;
        std::abort();
      }

      lout << onika::format_string("\nReading done in %ld ms | %ld µs | %ld ns\n", ctx.elapsed_time_ms,
                                   ctx.elapsed_time_ys, ctx.elapsed_time_ns)
           << std::endl;

      // ----------------------------------------- 
      // Now do exastamp stuff
      
      ctx.init_domain(domain, *bounds_mode);

      lout << std::endl;
      lout << "Uniform scale    = " << std::boolalpha << ctx.uniform_scale << std::endl;
      lout << "Particles        = " << particle_data.size() << std::endl;
      lout << "Domain XForm     = " << domain.xform() << std::endl;
      lout << "Domain bounds    = " << domain.bounds() << std::endl;
      lout << "Domain size      = " << bounds_size(domain.bounds()) << std::endl;
      lout << "Real size        = " << bounds_size(domain.bounds()) * Vec3d{ domain.xform().m11, domain.xform().m22, domain.xform().m33 } << std::endl;
      lout << "Cell size        = " << domain.cell_size() << std::endl;
      lout << "Grid dimensions  = " << domain.grid_dimension() << " (" << grid_cell_count(domain.grid_dimension()) << " cells)" << std::endl;

    }

    // send bounds and box_size to all cores
    MPI_Bcast(&domain, sizeof(Domain), MPI_CHARACTER, 0, *mpi);
    grid.set_origin(domain.bounds().bmin);
    grid.set_offset(IJK{0, 0, 0});
    grid.set_cell_size(domain.cell_size());
    grid.set_dimension(domain.grid_dimension());

    if (rank == 0) {

      spin_mutex_array cell_locks;
      cell_locks.resize(grid.number_of_cells());
      auto cells = grid.cells();

      #pragma omp parallel for
      for (size_t i = 0; i < particle_data.size(); ++i) {

        ParticleTupleIO& p = particle_data[i];

        Vec3d r = {p[field::rx], p[field::ry], p[field::rz]};
        r = ctx.invM * r;

        if (r.x < 0.)
          r.x += 1.;
        if (r.y < 0.)
          r.y += 1.;
        if (r.z < 0.)
          r.z += 1.;

        if (r.x >= 1.)
          r.x -= 1.;
        if (r.y >= 1.)
          r.y -= 1.;
        if (r.z >= 1.)
          r.z -= 1.;

        r.x *= ctx.boxlen.x;
        r.y *= ctx.boxlen.y;
        r.z *= ctx.boxlen.z;

        r = ctx.invXform * r;

        if (ctx.uniform_scale)
          r = r * ctx.domainXform.m11;

        p[field::rx] = r.x;
        p[field::ry] = r.y;
        p[field::rz] = r.z;

        IJK loc = domain_periodic_location(domain, r);
        size_t cell_index = grid_ijk_to_index(domain.grid_dimension(), loc);

        cell_locks[cell_index].lock();
        cells[cell_index].push_back(p);
        cell_locks[cell_index].unlock();
      }
    }

    grid.rebuild_particle_offsets();
#ifndef NDEBUG
    bool particles_inside_cell = check_particles_inside_cell(grid);
    assert(particles_inside_cell);
#endif
    lout << "===========================================" << std::endl << std::endl;
  }
};

template <class GridT> using ReadExternalAtomFormatNodeTmpl = ReadExternalAtomFormatNode<GridT>;

// === register factories ===
ONIKA_AUTORUN_INIT(read_external_atom_format) {
  OperatorNodeFactory::instance()->register_factory("read_external_atom_format",
                                                    make_grid_variant_operator<ReadExternalAtomFormatNodeTmpl>);
}

} // exaStamp

