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

#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_stream.h>
#include <exanb/core/domain.h>
#include <onika/string_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/grid_fields.h>

#include <mpi.h>
#include <string>
#include <cstring>
#include <iomanip>
#include <cmath>

namespace exaStamp
{

  using namespace exanb;

  template< class GridT >
  class WriteLMP : public OperatorNode
  {
    // field type
    using has_type_field_t = typename GridT::CellParticles::template HasField< field::_type >;
    static constexpr bool has_type_field = has_type_field_t::value;

    // field id
    using has_id_field_t = typename GridT::CellParticles::template HasField< field::_id >;
    static constexpr bool has_id_field = has_id_field_t::value;

    ADD_SLOT( MPI_Comm        , mpi          , INPUT , REQUIRED );
    ADD_SLOT( GridT           , grid         , INPUT , REQUIRED );
    ADD_SLOT( Domain          , domain       , INPUT , REQUIRED );
    ADD_SLOT( std::string     , filename     , INPUT , std::string("output.lmp") , DocString{"Output file name"} );
    ADD_SLOT( long            , timestep     , INPUT , REQUIRED );
    ADD_SLOT( ParticleSpecies , species      , INPUT , REQUIRED );
    ADD_SLOT( bool            , triclinic         , INPUT , false , DocString{"Write triclinic box (xy xz yz tilt factors). False = orthogonal (or auto-detect)."} );
    ADD_SLOT( bool            , write_velocities  , INPUT , true  , DocString{"Write the Velocities section. Set to false to omit velocities from the output file."} );

  public:

    inline void execute () override final
    {
      // ------------------------------------------------------------------
      // Fixed-width format strings.
      // atom-id and atom-type are written first, then x y z.
      // The trailing space is intentional: we later replace it with '\n'.
      // ------------------------------------------------------------------
      static const char* data_line_format = "%12llu %5u  % .16e % .16e % .16e ";
      static const std::string data_line_ref =
          onika::format_string(data_line_format, 1ULL, 1U, -1.0e103, -1.0e103, -1.0e103);
      static const size_t data_line_size = data_line_ref.length();

      static const char* vel_line_format = "%12llu % .16e % .16e % .16e ";
      static const std::string vel_line_ref =
          onika::format_string(vel_line_format, 1ULL, -1.0e103, -1.0e103, -1.0e103);
      static const size_t vel_line_size = vel_line_ref.length();

      // ------------------------------------------------------------------
      // MPI setup  –– use *mpi everywhere, never MPI_COMM_WORLD directly
      // ------------------------------------------------------------------
      int rank = 0, nprocs = 1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &nprocs);

      size_t n_cells = grid->number_of_cells();
      auto   cells   = grid->cells();

      // ------------------------------------------------------------------
      // Count local / total particles
      // ------------------------------------------------------------------
      unsigned long long local_nb_particles = 0;
      for (size_t c = 0; c < n_cells; ++c)
        if (!grid->is_ghost_cell(c))
          local_nb_particles += cells[c].size();

      unsigned long long total_nb_particles = 0;
      MPI_Allreduce(&local_nb_particles, &total_nb_particles, 1,
                    MPI_UNSIGNED_LONG_LONG, MPI_SUM, *mpi);

      // Exclusive prefix sum → first global atom index owned by this rank
      unsigned long long global_atom_start = 0;
      MPI_Scan(&local_nb_particles, &global_atom_start, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, *mpi);
      global_atom_start -= local_nb_particles;

      // ------------------------------------------------------------------
      // Box geometry
      //
      // LAMMPS requires an orthogonal or "restricted triclinic" description:
      //
      //   xlo_bound  xhi_bound   (+ optional xy)
      //   ylo_bound  yhi_bound   (+ optional xz)
      //   zlo_bound  zhi_bound   (+ optional yz)
      //
      // For a general triclinic cell the transformation matrix H maps
      // fractional → Cartesian.  LAMMPS stores it in lower-triangular
      // "a, b, c" form (a along x, b in xy-plane, c general).
      //
      // The domain xform() gives the 3×3 matrix A such that
      //   r_cart = A * r_frac
      // with columns a1=(Axx,Ayx,Azx), a2=(Axy,Ayy,Azy), a3=(Axz,Ayz,Azz).
      //
      // LAMMPS tilt factors:
      //   xy = dot(a2, a1_hat)
      //   xz = dot(a3, a1_hat)
      //   yz = dot(a3, b2_hat)
      // where a1_hat = a1/|a1|, b2 = a2 – xy*a1_hat (orthogonalised).
      //
      // For a purely orthogonal box all off-diagonal elements are 0.
      // ------------------------------------------------------------------

      const auto& xform = domain->xform(); // 3x3 matrix: r_cart = xform * r_frac
      const Vec3d origin = xform * domain->origin();

      // Lattice vectors: transform each axis of the domain extent into Cartesian space.
      // domain->extent() holds the box size along each fractional axis, so mapping
      // (Lx,0,0), (0,Ly,0), (0,0,Lz) through xform gives the true edge vectors.
      const Vec3d ext = domain->extent();
      const Vec3d a1 = xform * Vec3d{ext.x, 0.0,   0.0  };  // first  lattice vector
      const Vec3d a2 = xform * Vec3d{0.0,   ext.y, 0.0  };  // second lattice vector
      const Vec3d a3 = xform * Vec3d{0.0,   0.0,   ext.z};  // third  lattice vector

      // LAMMPS "restricted triclinic" parameters
      const double lx = std::sqrt(a1.x*a1.x + a1.y*a1.y + a1.z*a1.z);
      const double xy = (a2.x*a1.x + a2.y*a1.y + a2.z*a1.z) / lx;
      const double ly = std::sqrt(std::max(0.0, a2.x*a2.x + a2.y*a2.y + a2.z*a2.z - xy*xy));

      const double xz = (a3.x*a1.x + a3.y*a1.y + a3.z*a1.z) / lx;
      const double yz = (ly > 0.0)
          ? ((a3.x*a2.x + a3.y*a2.y + a3.z*a2.z) - xy*xz) / ly
          : 0.0;
      const double lz = std::sqrt(std::max(0.0,
          a3.x*a3.x + a3.y*a3.y + a3.z*a3.z - xz*xz - yz*yz));

      // Decide whether to write triclinic header
      const double tilt_tol = 1.0e-12 * std::max({lx, ly, lz});
      const bool is_triclinic = *triclinic ||
          (std::abs(xy) > tilt_tol || std::abs(xz) > tilt_tol || std::abs(yz) > tilt_tol);

      // LAMMPS bounding box (shifted from origin)
      const double xlo = origin.x;
      const double ylo = origin.y;
      const double zlo = origin.z;

      // For triclinic, LAMMPS expects "bound" values:
      //   xlo_bound = xlo + min(0, xy, xz, xy+xz)
      //   xhi_bound = xhi + max(0, xy, xz, xy+xz)   etc.
      // For orthogonal boxes these corrections vanish.
      const double xhi_raw = xlo + lx;
      const double yhi_raw = ylo + ly;
      const double zhi_raw = zlo + lz;

      double xlo_b = xlo, xhi_b = xhi_raw;
      double ylo_b = ylo, yhi_b = yhi_raw;
      double zlo_b = zlo, zhi_b = zhi_raw;

      if (is_triclinic)
      {
        xlo_b += std::min({0.0, xy, xz, xy + xz});
        xhi_b += std::max({0.0, xy, xz, xy + xz});
        ylo_b += std::min(0.0, yz);
        yhi_b += std::max(0.0, yz);
        // zlo / zhi unchanged
      }

      // ------------------------------------------------------------------
      // Build header string (rank 0 will write it)
      // ------------------------------------------------------------------
      // Rule: first line is always a comment in LAMMPS data files.
      // Ovito and LAMMPS both require this.
      std::string domain_header = onika::format_string(
          "LAMMPS data file via exaStamp V3, timestep = %ld\n"
          "\n"
          "%12llu  atoms\n"
          "%12llu  atom types\n"
          "\n"
          "%.12f\t%.12f  xlo xhi\n"
          "%.12f\t%.12f  ylo yhi\n"
          "%.12f\t%.12f  zlo zhi\n",
          *timestep,
          static_cast<unsigned long long>(total_nb_particles),
          static_cast<unsigned long long>(species->size()),
          xlo_b, xhi_b,
          ylo_b, yhi_b,
          zlo_b, zhi_b);

      if (is_triclinic)
      {
        domain_header += onika::format_string(
            "%.12f\t%.12f\t%.12f  xy xz yz\n",
            xy, xz, yz);
      }

      domain_header += "\n";

      std::string species_header = "Masses\n\n";
      int sp_count = 0;
      for (const auto& s : *species)
        species_header += onika::format_string("%12d %.12f # %s\n",
                                               ++sp_count, s.m_mass, s.name());
      species_header += "\nAtoms # atomic\n\n";

      const std::string header = domain_header + species_header;
      const unsigned long long offset_header = static_cast<unsigned long long>(header.length());

      // ------------------------------------------------------------------
      // Open file  –– use *mpi (not MPI_COMM_WORLD)
      // ------------------------------------------------------------------
      MPI_File mpiFile;
      MPI_Status status;
      MPI_File_open(*mpi, filename->c_str(),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiFile);

      // Rank 0 writes the header
      if (rank == 0)
      {
        ldbg << "LMP: write header" << std::endl;
        MPI_File_write_at(mpiFile, 0,
                          header.data(), static_cast<int>(header.length()),
                          MPI_CHAR, &status);
      }

      // ------------------------------------------------------------------
      // Write atom positions
      // Each line has a fixed byte width (data_line_size) so every rank
      // can independently compute its file offset.
      // ------------------------------------------------------------------
      ldbg << "LMP: write positions @" << offset_header << std::endl;
      {
        unsigned long long atom_count = 0;
        for (size_t c = 0; c < n_cells; ++c)
        {
          if (!grid->is_ghost_cell(c))
          {
            const auto*    rx    = cells[c][field::rx];
            const auto*    ry    = cells[c][field::ry];
            const auto*    rz    = cells[c][field::rz];
            const uint8_t* types = cells[c].field_pointer_or_null(field::type);
            const size_t*  ids   = cells[c].field_pointer_or_null(field::id);

            const int np_cell = static_cast<int>(cells[c].size());
            for (int p = 0; p < np_cell; ++p)
            {
              Vec3d r = domain->xform() * Vec3d{rx[p], ry[p], rz[p]};

              // atom type: 1-based; default to 1 if field absent
              unsigned int atom_type = (types != nullptr) ? static_cast<unsigned int>(types[p]) + 1u : 1u;

              // atom id: use field if available, otherwise use global sequential index (1-based)
              unsigned long long atom_id = (ids != nullptr)
                  ? static_cast<unsigned long long>(ids[p]) + 1ULL
                  : global_atom_start + atom_count + 1ULL;

              auto data_line = onika::format_string(data_line_format,
                                                    atom_id, atom_type,
                                                    r.x, r.y, r.z);
              assert(data_line.length() <= data_line_size);
              assert(!data_line.empty());
              assert(data_line.back() == ' ');
              data_line.resize(data_line_size, ' ');
              data_line.back() = '\n';

              const unsigned long long offset =
                  offset_header + (global_atom_start + atom_count) * data_line_size;
              MPI_File_write_at(mpiFile, static_cast<MPI_Offset>(offset),
                                data_line.data(), static_cast<int>(data_line.length()),
                                MPI_CHAR, &status);
              ++atom_count;
            }
          }
        }
      }

      // ------------------------------------------------------------------
      // Velocities section  –– skipped entirely if write_velocities is false
      // ------------------------------------------------------------------
      if (!*write_velocities)
      {
        MPI_File_close(&mpiFile);
        return;
      }

      // Velocities section header  –– only rank 0 writes it
      const std::string velocity_header = "\nVelocities\n\n";
      const unsigned long long offset_velocities_header =
          offset_header + total_nb_particles * data_line_size;
      ldbg << "LMP: write velocities header @" << offset_velocities_header << std::endl;

      if (rank == 0)   // BUG FIX: was written by ALL ranks
      {
        MPI_File_write_at(mpiFile,
                          static_cast<MPI_Offset>(offset_velocities_header),
                          velocity_header.data(),
                          static_cast<int>(velocity_header.length()),
                          MPI_CHAR, &status);
      }

      const unsigned long long offset_velocities =
          offset_velocities_header + static_cast<unsigned long long>(velocity_header.length());

      // ------------------------------------------------------------------
      // Write atom velocities
      // The atom-id written here must match the one used in Atoms section.
      // ------------------------------------------------------------------
      ldbg << "LMP: write velocities @" << offset_velocities << std::endl;
      {
        unsigned long long atom_count = 0;
        for (size_t c = 0; c < n_cells; ++c)
        {
          if (!grid->is_ghost_cell(c))
          {
            const auto*   vx   = cells[c][field::vx];
            const auto*   vy   = cells[c][field::vy];
            const auto*   vz   = cells[c][field::vz];
            const size_t* ids  = cells[c].field_pointer_or_null(field::id);

            const int np_cell = static_cast<int>(cells[c].size());
            for (int p = 0; p < np_cell; ++p)
            {
              // Must use the same id scheme as in the Atoms section
              unsigned long long atom_id = (ids != nullptr)
                  ? static_cast<unsigned long long>(ids[p]) + 1ULL
                  : global_atom_start + atom_count + 1ULL;

              auto data_line = onika::format_string(vel_line_format,
                                                    atom_id,
                                                    vx[p], vy[p], vz[p]);
              assert(data_line.length() <= vel_line_size);
              assert(!data_line.empty());
              assert(data_line.back() == ' ');
              data_line.resize(vel_line_size, ' ');
              data_line.back() = '\n';

              const unsigned long long offset =
                  offset_velocities + (global_atom_start + atom_count) * vel_line_size;
              MPI_File_write_at(mpiFile, static_cast<MPI_Offset>(offset),
                                data_line.data(), static_cast<int>(data_line.length()),
                                MPI_CHAR, &status);
              ++atom_count;
            }
          }
        }
      }

      MPI_File_close(&mpiFile);
    }
  };

  template<class GridT> using WriteLMPTmpl = WriteLMP<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(write_lmp)
  {
    OperatorNodeFactory::instance()->register_factory(
        "write_lmp", make_grid_variant_operator< WriteLMPTmpl >);
  }

} // namespace exaStamp
