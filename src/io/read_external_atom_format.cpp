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

    // declare holders that will be fill by parsing files
    Domain& domain = *(this->domain);
    GridT& grid = *(this->grid);
    ParticleSpecies species = *(this->species);
    ParticleData particle_data;
    Mat3d H = make_identity_matrix();
    size_t n_particles = 0;

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

    Context ctx{
        .domain = domain,
        .H = H,
        .particle_data = particle_data,
        .n_particles = n_particles,
    };

    if ( this->species.has_value() ) {
      for (size_t i = 0; i < species.size(); i++) {
        ctx.species.insert({ species.at(i).m_name, i });
      }
      ctx.next_type_id = species.size();
    }

    // from here only one proc needs to work
    Mat3d Ht, inv_Ht, xform, inv_xform;
    double box_size_x, box_size_y, box_size_z;
    bool uniform_scale;

    if ( rank == 0 ) {

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
        std::abort();
      }

      lout << onika::format_string("\nParsing done in %ld ms | %ld µs | %ld ns\n", ctx.elapsed_time_ms, ctx.elapsed_time_ys, ctx.elapsed_time_ns) << std::endl;

      // ----------------------------------------- 
      // Now do exastamp stuff

      Vec3d a{H.m11, H.m12, H.m13};
      Vec3d b{H.m21, H.m22, H.m23};
      Vec3d c{H.m31, H.m32, H.m33};

      box_size_x = norm(a);
      box_size_y = norm(b);
      box_size_z = norm(c);

      lout << "a = " << a << onika::format_string(" , norm = %10.5f ", box_size_x) << std::endl;
      lout << "b = " << b << onika::format_string(" , norm = %10.5f ", box_size_y) << std::endl;
      lout << "c = " << c << onika::format_string(" , norm = %10.5f ", box_size_z) << std::endl;

      if (!domain.xform_is_identity()) {
        lout << std::endl << "Init initial xform to identity" << std::endl;
        domain.set_xform(make_identity_matrix());
      }

      AABB domain_bounds = {{0., 0., 0.}, {box_size_x, box_size_y, box_size_z}};
      compute_domain_bounds(domain, *bounds_mode, 0.0, domain_bounds, domain_bounds, true);

      Ht = transpose(H);
      inv_Ht = inverse(Ht);
      inv_xform = (domain.xform_is_identity()) ? make_identity_matrix() : domain.inv_xform();

      // a partir de la je comprend plus
      Mat3d D = diag_matrix(Vec3d{box_size_x, box_size_y, box_size_z});
      Mat3d G = Ht * inverse(D);
      xform = G * domain.xform();

      lout << Ht << std::endl;
      lout << G << std::endl;
      lout << D << std::endl;
      lout << xform << std::endl;

      uniform_scale = is_uniform_scale(xform);

      if (uniform_scale) {
        domain.set_xform(make_identity_matrix());
        domain.set_cell_size(domain.cell_size() * xform.m11);
        domain.set_bounds({domain.origin() * xform.m11, domain.extent() * xform.m11});
      } else {
        domain.set_xform(xform);
      }

      lout << std::endl;
      lout << "Uniform scale    = " << std::boolalpha << uniform_scale << std::endl;
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

    // add the particles to the grid
    using particle_tuple_t = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_vx, field::_vy,
                                                      field::_vz, field::_id, field::_type>;

    if (rank == 0) {
      for (size_t i = 0; i < particle_data.size(); ++i) {

        ParticleTupleIO& p = particle_data[i];

        Vec3d r = {p[field::rx], p[field::ry], p[field::rz]};
        r = inv_Ht * r;

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

        r.x *= box_size_x;
        r.y *= box_size_y;
        r.z *= box_size_z;

        r = inv_xform * r;

        if (uniform_scale)
          r = r * xform.m11;

        p[field::rx] = r.x;
        p[field::ry] = r.y;
        p[field::rz] = r.z;

        IJK loc = domain_periodic_location(domain, r);
        particle_tuple_t t = p;
        grid.cell(loc).push_back(t);
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

