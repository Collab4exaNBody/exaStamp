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
#include <onika/physics/units.h>
#include <onika/string_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/grid_fields.h>

#include <exanb/defbox/deformation.h>
#include <exanb/defbox/deformation_stream.h>
#include <exanb/defbox/deformation_yaml.h>
#include <exanb/defbox/deformation_math.h>

#include <onika/soatl/packed_field_arrays.h>
#include <memory>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <string>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <stdexcept>

namespace exaStamp
{
  using namespace exanb;

  template< class GridT >
  class WriteXYZSkinfilesOperator : public OperatorNode
  {
    // field type
    using has_type_field_t = typename GridT::CellParticles::template HasField< field::_type >;
    static constexpr bool has_type_field = has_type_field_t::value;

    // field id
    using has_id_field_t = typename GridT::CellParticles::template HasField< field::_id >;
    static constexpr bool has_id_field = has_id_field_t::value;

    // field id mol
    using has_id_mol_field_t = typename GridT::CellParticles::template HasField< field::_idmol >;
    static constexpr bool has_id_mol_field = has_id_mol_field_t::value;

    ADD_SLOT( MPI_Comm        , mpi           , INPUT );
    ADD_SLOT( GridT           , grid          , INPUT );
    ADD_SLOT( Domain          , domain        , INPUT );
    ADD_SLOT( Deformation     , defbox        , INPUT );
    ADD_SLOT( long            , timestep      , INPUT );
    ADD_SLOT( bool            , is_ghosts     , INPUT , false );
    ADD_SLOT( ParticleSpecies , species       , INPUT , REQUIRED );
    ADD_SLOT( GridT           , grid_t0       , INPUT , OPTIONAL );
    ADD_SLOT( double          , skin_distance , INPUT , 5.0 ,
              DocString{"Physical skin thickness in the same length units as the simulation box. "
                        "Particles whose reference-frame (t0) position lies within this distance "
                        "of any face of the undeformed box are written to the skin file."} );

  public:

    inline void execute () override final
    {
      namespace fs = std::filesystem;
      using std::string;
      using std::ostringstream;

      // ------------------------------------------------------------------
      // Guard: grid_t0 is OPTIONAL — bail out cleanly if not connected
      // ------------------------------------------------------------------
      if (!grid_t0.has_value())
      {
        lerr << "write_xyz_skin_file: grid_t0 slot not connected, skipping output." << std::endl;
        return;
      }

      GridT& grid_ref    = *(this->grid);
      GridT& grid_t0_ref = *(this->grid_t0);

      ParticleSpecies& species_ref = *(this->species);
      const bool write_ghosts = *(this->is_ghosts);
      const Mat3d xform = domain->xform();

      // ------------------------------------------------------------------
      // Reference-frame lattice matrices
      //
      // lot    : row-major lattice in the *current*   frame (deformed)
      // lot_t0 : row-major lattice in the *reference* frame (t0, undeformed)
      //
      // Since grid_t0 stores Cartesian positions in the undeformed box, the
      // t0 lattice is a pure diagonal matrix: identity orientation, box size.
      // We transpose it to match the column-major Lattice= convention used in
      // the extended XYZ header (and to stay consistent with 'lot').
      // ------------------------------------------------------------------
      const Vec3d box_size  = domain->extent() - domain->origin();
      const Mat3d lattice   = diag_matrix(box_size);
      const Mat3d lot       = transpose(xform * lattice);   // current frame, row-major
      const Mat3d lot_t0    = transpose(lattice);            // t0 frame, row-major (consistent)
      const Mat3d inv_lot_t0 = inverse(lot_t0);

      ldbg << "lattice : " << lattice  << std::endl;
      ldbg << "lot     : " << lot      << std::endl;
      ldbg << "lot_t0  : " << lot_t0   << std::endl;

      // ------------------------------------------------------------------
      // Convert physical skin_distance → fractional coordinate thresholds.
      //
      // For an orthorhombic undeformed box the fractional thickness along
      // each axis i is:  delta_i = skin_distance / box_size_i
      // We take the most conservative (largest) fractional thickness so
      // that the physical distance is guaranteed on every axis.
      //
      // frac_skin is clamped to [0, 0.5] so it never exceeds half the box.
      // ------------------------------------------------------------------
      const double sd = *skin_distance;
      if (sd < 0.0)
      {
        lerr << "write_xyz_skin_file: skin_distance must be >= 0, got "
             << sd << ". Skipping output." << std::endl;
        return;
      }

      // Per-axis fractional skin width (safe against zero box dimension)
      const double frac_x = (box_size.x > 0.0) ? (sd / box_size.x) : 0.5;
      const double frac_y = (box_size.y > 0.0) ? (sd / box_size.y) : 0.5;
      const double frac_z = (box_size.z > 0.0) ? (sd / box_size.z) : 0.5;

      // Per-axis skin bounds (clamped so low < high)
      const double skin_low_x  = std::min(frac_x, 0.5);
      const double skin_high_x = 1.0 - skin_low_x;
      const double skin_low_y  = std::min(frac_y, 0.5);
      const double skin_high_y = 1.0 - skin_low_y;
      const double skin_low_z  = std::min(frac_z, 0.5);
      const double skin_high_z = 1.0 - skin_low_z;

      ldbg << "skin_distance=" << sd
           << "  frac=(" << frac_x << ", " << frac_y << ", " << frac_z << ")" << std::endl;

      // ------------------------------------------------------------------
      // Cell count — use the same grid for both passes to guarantee that
      // the particle counts per cell seen in Pass 1 match those in Pass 2.
      // We iterate over grid_ref for the ghost-cell test (authoritative),
      // but read particle counts from cells_t0 in Pass 1 and cells in Pass 2.
      // To be safe we take min(np_t0, np_cur) so we never overrun either.
      // ------------------------------------------------------------------
      const size_t n_cells = grid_ref.number_of_cells();
      auto cells    = grid_ref.cells();
      auto cells_t0 = grid_t0_ref.cells();

      // ------------------------------------------------------------------
      // MPI setup — always use *mpi, never MPI_COMM_WORLD
      // ------------------------------------------------------------------
      int rank = 0, nprocs = 1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &nprocs);

      // ------------------------------------------------------------------
      // Output directory — recreated once per run on rank 0.
      // Use a member flag (not a static local) so multiple operator
      // instances each manage their own directory lifecycle.
      // ------------------------------------------------------------------
      ldbg << "create xyz_datas_skin dir" << std::endl;
      if (rank == 0 && !m_dir_created)
      {
        m_dir_created = true;
        fs::remove_all("xyz_datas_skin");
        std::error_code ec;
        fs::create_directory("xyz_datas_skin", ec);
        if (ec)
        {
          lerr << "write_xyz_skin_file: failed to create output directory: "
               << ec.message() << std::endl;
        }
      }

      // All ranks must wait until rank 0 has created the directory
      MPI_Barrier(*mpi);

      // ------------------------------------------------------------------
      // Build filename
      // ------------------------------------------------------------------
      ostringstream filename_xyz;
      filename_xyz << "xyz_datas_skin/atomic_structure_"
                   << std::setfill('0') << std::setw(15) << std::to_string(*timestep)
                   << ".xyz";

      // ------------------------------------------------------------------
      // Build a fixed-width format reference line so every particle line
      // occupies exactly the same number of bytes in the file.
      // The element name field is padded to 10 chars; coords use % .10e.
      // ------------------------------------------------------------------
      static const char* xyz_line_format = "%-10s % .10e % .10e % .10e\n";
      static const std::string xyz_line_ref =
          onika::format_string(xyz_line_format, "XX", -1.0e103, -1.0e103, -1.0e103);
      static const unsigned long long xyz_line_size =
          static_cast<unsigned long long>(xyz_line_ref.length());

      // Helper lambda: returns true if the t0 fractional position is in the skin
      // (i.e. within skin_distance of any face of the undeformed box).
      auto in_skin = [&](double rx0, double ry0, double rz0) -> bool
      {
        const Vec3d frac = inv_lot_t0 * Vec3d{rx0, ry0, rz0};
        return frac.x < skin_low_x  || frac.x > skin_high_x ||
               frac.y < skin_low_y  || frac.y > skin_high_y ||
               frac.z < skin_low_z  || frac.z > skin_high_z;
      };

      // ------------------------------------------------------------------
      // Pass 1: count skin particles on this rank.
      // Use cells_t0 particle count as the reference; both grids should
      // have the same particle layout at t0, but we guard with the minimum.
      // ------------------------------------------------------------------
      unsigned long long local_skin_count = 0;
      for (size_t c = 0; c < n_cells; ++c)
      {
        if (!grid_ref.is_ghost_cell(c) || write_ghosts)
        {
          const auto* rx0 = cells_t0[c][field::rx];
          const auto* ry0 = cells_t0[c][field::ry];
          const auto* rz0 = cells_t0[c][field::rz];
          // Guard: use the smaller of the two cell particle counts
          const int np_t0  = static_cast<int>(cells_t0[c].size());
          const int np_cur = static_cast<int>(cells[c].size());
          const int np     = std::min(np_t0, np_cur);
          for (int p = 0; p < np; ++p)
          {
            if (in_skin(rx0[p], ry0[p], rz0[p]))
              ++local_skin_count;
          }
        }
      }

      // Global total and exclusive prefix sum (offset for this rank)
      unsigned long long total_skin_count = 0;
      MPI_Allreduce(&local_skin_count, &total_skin_count, 1,
                    MPI_UNSIGNED_LONG_LONG, MPI_SUM, *mpi);

      unsigned long long skin_offset = 0;   // first global index owned by this rank
      MPI_Scan(&local_skin_count, &skin_offset, 1,
               MPI_UNSIGNED_LONG_LONG, MPI_SUM, *mpi);
      skin_offset -= local_skin_count;

      // ------------------------------------------------------------------
      // Build XYZ header (rank 0 will write it).
      // Extended XYZ format: first line = atom count, second line = comment
      // with Lattice= key for Ovito / VESTA etc.
      // ------------------------------------------------------------------
      const string header = onika::format_string(
          "%llu\n"
          "Lattice=\"%.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e %.12e\"\n",
          static_cast<unsigned long long>(total_skin_count),
          lot.m11, lot.m12, lot.m13,
          lot.m21, lot.m22, lot.m23,
          lot.m31, lot.m32, lot.m33);

      const unsigned long long offset_header =
          static_cast<unsigned long long>(header.length());

      // ------------------------------------------------------------------
      // Open file — use *mpi, not MPI_COMM_WORLD.
      // Truncate any pre-existing file to the exact expected size so no
      // stale bytes from a previous (longer) file survive.
      // ------------------------------------------------------------------
      MPI_File mpiFile;
      MPI_Status status;
      MPI_File_open(*mpi, filename_xyz.str().c_str(),
                    MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiFile);

      // Truncate to the final expected file size to remove stale data
      const MPI_Offset expected_file_size =
          static_cast<MPI_Offset>(offset_header + total_skin_count * xyz_line_size);
      MPI_File_set_size(mpiFile, expected_file_size);

      // Rank 0 writes the header using MPI_File_write_at (positional, not collective)
      if (rank == 0)
      {
        ldbg << "XYZ skin: write header" << std::endl;
        MPI_File_write_at(mpiFile, 0,
                          header.data(), static_cast<int>(header.length()),
                          MPI_CHAR, &status);
      }

      // ------------------------------------------------------------------
      // Pass 2: write skin particle lines at their correct file offsets.
      // skin_write_count tracks how many skin particles THIS rank has
      // written so far — it is independent of the counting pass above.
      // ------------------------------------------------------------------
      ldbg << "XYZ skin: write positions @" << offset_header << std::endl;
      unsigned long long skin_write_count = 0;

      for (size_t c = 0; c < n_cells; ++c)
      {
        if (!grid_ref.is_ghost_cell(c) || write_ghosts)
        {
          const auto*    rx    = cells[c][field::rx];
          const auto*    ry    = cells[c][field::ry];
          const auto*    rz    = cells[c][field::rz];
          const auto*    rx0   = cells_t0[c][field::rx];
          const auto*    ry0   = cells_t0[c][field::ry];
          const auto*    rz0   = cells_t0[c][field::rz];
          const uint8_t* types = cells[c].field_pointer_or_null(field::type);
          ONIKA_ASSUME_ALIGNED(types);

          // Mirror the same guard used in Pass 1 so counts stay in sync
          const int np_t0  = static_cast<int>(cells_t0[c].size());
          const int np_cur = static_cast<int>(cells[c].size());
          const int np     = std::min(np_t0, np_cur);

          for (int p = 0; p < np; ++p)
          {
            if (!in_skin(rx0[p], ry0[p], rz0[p])) continue;

            // Current Cartesian position (deformed frame)
            const Vec3d pos = xform * Vec3d{rx[p], ry[p], rz[p]};

            const char* type_name = "XX";
            if (types != nullptr) type_name = species_ref[types[p]].m_name;

            auto line = onika::format_string(xyz_line_format,
                                             type_name, pos.x, pos.y, pos.z);

            // Strip the trailing '\n' produced by xyz_line_format before
            // padding. Without this, resize() would leave the original '\n'
            // in the middle of the slot and line.back()='\n' would add a
            // second one, producing a blank line in the output file.
            if (!line.empty() && line.back() == '\n')
              line.pop_back();

            // Runtime check: content (without the newline) must fit in
            // (xyz_line_size - 1) bytes, leaving room for exactly one '\n'.
            if (line.length() > xyz_line_size - 1)
            {
              lerr << "write_xyz_skin_file: formatted line exceeds fixed slot size ("
                   << line.length() << " > " << (xyz_line_size - 1)
                   << "). Truncating to avoid file corruption." << std::endl;
              line.resize(xyz_line_size - 1);
            }

            // Pad content with spaces to (xyz_line_size - 1), then close
            // with exactly one '\n' — every slot is xyz_line_size bytes total.
            line.resize(xyz_line_size - 1, ' ');
            line += '\n';

            const unsigned long long file_offset =
                offset_header + (skin_offset + skin_write_count) * xyz_line_size;

            MPI_File_write_at(mpiFile,
                              static_cast<MPI_Offset>(file_offset),
                              line.data(), static_cast<int>(line.length()),
                              MPI_CHAR, &status);
            ++skin_write_count;
          }
        }
      }

      MPI_File_close(&mpiFile);

      ldbg << "XYZ skin: wrote " << skin_write_count
           << " particles (rank " << rank << ")" << std::endl;
    }

  private:
    // Member flag instead of static local — each operator instance tracks its own state
    bool m_dir_created = false;
  };

  template<class GridT> using WriteXYZSkinfilesOperatorTmpl = WriteXYZSkinfilesOperator<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(write_xyz_skin_file)
  {
    OperatorNodeFactory::instance()->register_factory(
        "write_xyz_skin_file", make_grid_variant_operator< WriteXYZSkinfilesOperatorTmpl >);
  }

} // namespace exaStamp
