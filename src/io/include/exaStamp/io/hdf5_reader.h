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

// HDF5 atom reader (read_hdf5, aio ADR 0016, 0017).
//
// GridBindSink stages atoms in contiguous arrays (r[3n] or x[n] y[n] z[n],
// type[n], ids[n], properties; one (N,3) velocity dataset goes straight into
// GridFiller), so aio reads each dataset with one H5Dread straight into them
// (aio ADR 0015). end() then places each atom with GridPlacer and stages it in
// GridFiller, as GridStreamSink does for the text formats (atom_file_reader.h).
//
// Flows (read_into_grid, shared with the text readers):
//   - MPIO (parallel HDF5, mpio: true): every rank reads its slab of rows and
//     fills its own grid (full-domain geometry, as the other ranks get after a
//     rank-0 read). The deck's load balancing / migration redistributes the
//     atoms afterwards, as after any read.
//   - serial: rank 0 reads. Used when HDF5 has no parallel I/O (with a warning
//     if several ranks run) or with mpio: false.
//   - stream_sink: true (debug / benchmark): read_atom_file with GridStreamSink.
// File ids are checked for range and uniqueness (GridBindSink::global_check).
#pragma once

#include <cstdint>
#include <cstring>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <mpi.h>
#include <yaml-cpp/yaml.h>

#include <onika/log.h>

#include <exaStamp/io/atom_file_reader.h>

namespace exaStamp::atom_io
{

  // The `schema` slot of read_hdf5: dataset paths of the aio HDF5 reader.
  struct Hdf5Schema
  {
    aio::formats::hdf5::Schema s;
  };

  inline std::ostream &operator<<(std::ostream &out, const aio::formats::hdf5::Vec3Path &v)
  {
    if (!v.path.empty())
      return out << v.path;
    return out << "[" << v.xyz[0] << ", " << v.xyz[1] << ", " << v.xyz[2] << "]";
  }

  inline std::ostream &operator<<(std::ostream &out, const Hdf5Schema &h)
  {
    const auto &s = h.s;
    out << "{ positions: " << s.positions;
    if (!s.velocities.empty())
      out << ", velocities: " << s.velocities;
    for (const auto &[k, v] : {std::pair<const char *, const std::string *>{"ids", &s.ids}, {"types", &s.types}, {"cell", &s.cell}, {"origin", &s.origin}, {"pbc", &s.pbc}})
      if (!v->empty())
        out << ", " << k << ": " << *v;
    if (!s.properties.empty())
    {
      out << ", fields: {";
      for (std::size_t k = 0; k < s.properties.size(); ++k)
        out << (k ? ", " : " ") << s.properties[k].name << ": " << s.properties[k].path;
      out << " }";
    }
    return out << " }";
  }

  // GridBindSink: contiguous staging for direct HDF5 reads (aio ADR 0015), then
  // the placement and GridFiller staging shared with GridStreamSink.
  template <class GridT> class GridBindSink
  {
  public:
    GridBindSink(Domain &domain, GridT &grid, const ReaderSettings &settings, bool interleaved_pos, bool interleaved_vel, bool ids_from_file) : m_domain(domain), m_grid(grid), m_species(settings.species), m_bounds_mode(settings.bounds_mode), m_read_velocities(settings.read_velocities), m_file_has_pbc(settings.file_has_pbc), m_props(settings.properties, reserved_field_names(typename GridT::Fields{})), m_interleaved_pos(interleaved_pos), m_interleaved_vel(interleaved_vel), m_ids_from_file(ids_from_file) {}

    aio::Status begin(const aio::atoms::FrameInfo &info)
    {
      m_natoms = info.natoms;
      m_natoms_global = info.natoms_global ? info.natoms_global : info.natoms;
      m_first_row = info.first_row;
      if (aio::Status st = m_place.begin(info, m_domain, m_bounds_mode, m_file_has_pbc, "schema: cell"); !st.ok())
        return st;
      set_grid_geometry(m_grid, m_domain);
      m_filler.reset(m_natoms);
      return {};
    }

    // Contiguous staging, allocated for the fields the schema has (aio binds those only).
    // Not zero-filled: aio's direct reads write every value. One (N,3) velocity
    // dataset is read straight into GridFiller's velocity array.
    aio::atoms::FieldBinding bind(aio::atoms::Field f)
    {
      using F = aio::atoms::Field;
      using B = aio::atoms::FieldBinding;
      switch (f)
      {
      case F::Id:
        if (!m_ids_from_file)
          return {};
        return B::of(alloc(m_ids, m_natoms));
      case F::Type:
        return B::of(alloc(m_type, m_natoms));
      case F::X:
      case F::Y:
      case F::Z:
        return vec_binding(m_pos, m_interleaved_pos, static_cast<int>(f) - static_cast<int>(F::X));
      case F::VX:
      case F::VY:
      case F::VZ:
        if (!m_read_velocities)
          return {};
        if (m_interleaved_vel)
          return B::of(m_filler.velocity_storage(false), 3, static_cast<int>(f) - static_cast<int>(F::VX));
        return vec_binding(m_vel, false, static_cast<int>(f) - static_cast<int>(F::VX));
      default:
        return {};
      }
    }

    aio::atoms::FieldBinding bind_property(std::string_view name, std::size_t ncomp, std::size_t comp, aio::atoms::ScalarType) { return m_props.bind(m_filler, name, ncomp, comp); }

    // All datasets read: type map, id range, placement, then GridFiller staging.
    aio::Status end()
    {
      if (aio::Status st = m_props.check(m_msg, "not in the schema or read as a grid field"); !st.ok())
        return st;
      // Types: t -> species t - 1 (no labels in HDF5). Without species, the largest type sets the count.
      std::size_t ntypes = m_species ? m_species->size() : 1;
      if (!m_species && m_type)
        for (std::size_t i = 0; i < m_natoms; ++i)
          ntypes = std::max<std::size_t>(ntypes, static_cast<std::size_t>(std::max<std::int32_t>(m_type[i], 1)));
      if (!m_types.freeze(ntypes, m_species, m_msg))
        return aio::make_error(aio::Errc::Sink, 0, 0, "atom types");
      const aio::Status st = m_interleaved_pos ? stage_all<true>() : stage_all<false>();
      for (auto &a : m_pos)
        a.reset();
      for (auto &a : m_vel)
        a.reset();
      m_type.reset();
      if (st.ok())
        m_count = m_natoms;
      return st;
    }

    // After a successful read on every reading rank (collective when parallel).
    // File ids are in [1, N] (end()) and the rows add up to N, so they are unique
    // iff no id repeats. Sorted files (id = row + 1 everywhere) are not checked.
    // Serial: a bitset of N bits. Parallel: id k goes to rank k / B, B = ceil(N / P),
    // which tests it in a bitset of B bits (4 B sent per atom).
    std::string global_check(MPI_Comm comm, bool parallel) const
    {
      if (!m_ids_from_file)
        return {}; // remap_atom: id = global row
      int sorted = m_ids_are_rows ? 1 : 0;
      if (parallel)
        MPI_Allreduce(MPI_IN_PLACE, &sorted, 1, MPI_INT, MPI_LAND, comm);
      if (sorted)
        return {};
      int nranks = 1;
      if (parallel)
        MPI_Comm_size(comm, &nranks);
      std::uint64_t dup = kNoId;
      if (nranks > 1)
      {
        dup = repeated_id_parallel(comm);
      }
      else
      {
        aio::formats::common::SeenSet seen(m_natoms_global);
        for (std::size_t i = 0; i < m_natoms && dup == kNoId; ++i)
          if (seen.test_and_set(static_cast<std::size_t>(m_ids[i] - 1)))
            dup = static_cast<std::uint64_t>(m_ids[i]);
      }
      if (dup == kNoId)
        return {};
      return "atom id " + std::to_string(dup) + " appears more than once in the file";
    }

    // Puts the staged atoms into the grid, with id = file id - 1, or the global row
    // (remap_atom, or file ids equal to the rows).
    void commit()
    {
      const bool by_row = !m_ids_from_file || m_ids_are_rows;
      m_filler.commit(m_grid, m_first_row, by_row ? nullptr : m_ids.get());
      m_ids.reset();
    }

    std::vector<GenericFieldInfo> fill_properties() { return m_filler.fill_properties(m_grid); }
    const std::vector<std::string> &skipped_properties() const { return m_props.skipped(); }
    const std::string &message() const { return m_msg; }
    std::size_t count() const { return m_count; }
    const BoxTransform &box() const { return m_place.box(); }

  private:
    static constexpr std::uint64_t kNoId = std::numeric_limits<std::uint64_t>::max();

    template <bool kInterleavedPos> aio::Status stage_all()
    {
      const bool split_vel = static_cast<bool>(m_vel[0]);
      for (std::size_t i = 0; i < m_natoms; ++i)
      {
        const std::int64_t t = m_type ? m_type[i] : 1; // no types in the schema: every atom is type 1
        if (t < 1 || static_cast<std::size_t>(t) > m_types.size())
        {
          m_msg = "atom type " + std::to_string(t) + " at row " + std::to_string(m_first_row + i) + " outside [1, " + std::to_string(m_types.size()) + "]";
          return aio::make_error(aio::Errc::Sink, 0, 0, "atom type out of range");
        }
        if (m_ids_from_file)
        {
          if (m_ids[i] < 1 || static_cast<std::uint64_t>(m_ids[i]) > m_natoms_global)
          {
            m_msg = "atom id " + std::to_string(m_ids[i]) + " at row " + std::to_string(m_first_row + i) + " outside [1, " + std::to_string(m_natoms_global) + "]";
            return aio::make_error(aio::Errc::IdOutOfRange, 0, 0, "atom id out of range");
          }
          m_ids_are_rows = m_ids_are_rows && static_cast<std::uint64_t>(m_ids[i]) == m_first_row + i + 1;
        }
        const Vec3d p = kInterleavedPos ? Vec3d{m_pos[0][3 * i], m_pos[0][3 * i + 1], m_pos[0][3 * i + 2]} : Vec3d{m_pos[0][i], m_pos[1][i], m_pos[2][i]};
        Vec3d r;
        std::uint32_t cell;
        if (!m_place.place(p, r, cell))
        {
          m_msg = "atom at row " + std::to_string(m_first_row + i) + " outside the domain on a non-periodic axis";
          return aio::make_error(aio::Errc::Sink, 0, 0, "atom outside the domain");
        }
        m_filler.stage(i, cell, r, m_types[static_cast<std::size_t>(t - 1)]);
        if (split_vel)
          m_filler.stage_velocity(i, m_vel[0][i], m_vel[1][i], m_vel[2][i]);
      }
      return {};
    }

    // shortcut: MPI int counts limit this to 2^31 atoms per rank; use MPI_Alltoallv_c (MPI 4) beyond.
    std::uint64_t repeated_id_parallel(MPI_Comm comm) const
    {
      int rank = 0, nranks = 1;
      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &nranks);
      const std::size_t P = static_cast<std::size_t>(nranks);
      const std::uint64_t B = (m_natoms_global + P - 1) / P; // ids per owner; < 2^32 (rows per rank are)
      // k / B without a 64-bit division per atom: reciprocal, corrected by one (exact for k < 2^53)
      const double inv_b = 1.0 / static_cast<double>(B);
      auto owner_of = [B, inv_b](std::uint64_t k)
      {
        std::uint64_t q = static_cast<std::uint64_t>(static_cast<double>(k) * inv_b);
        if (q * B > k)
          --q;
        else if ((q + 1) * B <= k)
          ++q;
        return q;
      };
      std::vector<int> scount(P, 0), rcount(P), sdispl(P + 1, 0), rdispl(P + 1, 0);
      for (std::size_t i = 0; i < m_natoms; ++i)
        ++scount[owner_of(static_cast<std::uint64_t>(m_ids[i] - 1))];
      MPI_Alltoall(scount.data(), 1, MPI_INT, rcount.data(), 1, MPI_INT, comm);
      for (std::size_t r = 0; r < P; ++r)
      {
        sdispl[r + 1] = sdispl[r] + scount[r];
        rdispl[r + 1] = rdispl[r] + rcount[r];
      }
      std::vector<std::uint32_t> sbuf(m_natoms), rbuf(static_cast<std::size_t>(rdispl[P]));
      {
        std::vector<int> next(sdispl.begin(), sdispl.end() - 1);
        for (std::size_t i = 0; i < m_natoms; ++i)
        {
          const std::uint64_t k = static_cast<std::uint64_t>(m_ids[i] - 1);
          const std::uint64_t owner = owner_of(k);
          sbuf[static_cast<std::size_t>(next[owner]++)] = static_cast<std::uint32_t>(k - owner * B);
        }
      }
      MPI_Alltoallv(sbuf.data(), scount.data(), sdispl.data(), MPI_UINT32_T, rbuf.data(), rcount.data(), rdispl.data(), MPI_UINT32_T, comm);
      std::uint64_t dup = kNoId;
      aio::formats::common::SeenSet seen(static_cast<std::size_t>(B));
      for (const std::uint32_t v : rbuf)
        if (seen.test_and_set(v))
        {
          dup = static_cast<std::uint64_t>(rank) * B + v + 1;
          break;
        }
      MPI_Allreduce(MPI_IN_PLACE, &dup, 1, MPI_UINT64_T, MPI_MIN, comm);
      return dup;
    }

    template <typename T> T *alloc(std::unique_ptr<T[]> &p, std::size_t n)
    {
      if (!p)
        p.reset(new T[n > 0 ? n : 1]);
      return p.get();
    }

    // Component c of a vector field: r[3n] interleaved (one (N,3) dataset) or x[n] y[n] z[n].
    aio::atoms::FieldBinding vec_binding(std::unique_ptr<double[]> (&a)[3], bool interleaved, int c)
    {
      if (interleaved)
        return aio::atoms::FieldBinding::of(alloc(a[0], 3 * m_natoms), 3, c);
      return aio::atoms::FieldBinding::of(alloc(a[c], m_natoms));
    }

    Domain &m_domain;
    GridT &m_grid;
    const ParticleSpecies *m_species;
    ReadBoundsSelectionMode m_bounds_mode;
    bool m_read_velocities = true;
    bool m_file_has_pbc = false;
    PropertyFilter m_props;
    bool m_interleaved_pos, m_interleaved_vel, m_ids_from_file;
    bool m_ids_are_rows = true; // file ids == global row + 1 (sorted files): no id lookup, no uniqueness check
    GridPlacer m_place;
    TypeMap m_types;
    GridFiller<GridT> m_filler;
    std::unique_ptr<double[]> m_pos[3], m_vel[3]; // m_vel: split x/y/z velocities only
    std::unique_ptr<std::int32_t[]> m_type;
    std::unique_ptr<std::int64_t[]> m_ids;
    std::size_t m_natoms = 0, m_natoms_global = 0, m_first_row = 0, m_count = 0;
    std::string m_msg;
  };

  // read_hdf5 flows (see the header comment). `mpio`: every rank reads its slab (needs
  // AIO_HDF5_HAVE_MPIO); otherwise rank 0 reads everything.
  template <class GridT> void read_hdf5_into_grid(MPI_Comm comm, Domain &domain, GridT &grid, const ReaderSettings &settings, const std::string &path, const aio::formats::hdf5::Schema &schema, bool remap_atom, bool mpio)
  {
#if !AIO_HDF5_HAVE_MPIO
    int nranks = 1;
    MPI_Comm_size(comm, &nranks);
    if (mpio && nranks > 1)
      lout << "Warning: read_hdf5: HDF5 built without parallel I/O (H5_HAVE_PARALLEL): reading on rank 0 only, then "
              "the atoms are distributed as after any rank-0 read"
           << std::endl;
    mpio = false;
#endif
    const aio::formats::hdf5::Options opt; // row order: ids are handled by the sink
    read_into_grid(
        "read_hdf5", comm, domain, grid, path, mpio, [&] { return GridBindSink<GridT>(domain, grid, settings, !schema.positions.path.empty(), !schema.velocities.path.empty(), !remap_atom); },
        [&](const char *p, GridBindSink<GridT> &sink)
        {
#if AIO_HDF5_HAVE_MPIO
          if (mpio)
            return aio::formats::hdf5::read_file_mpio(comm, p, schema, sink, opt);
#endif
          return aio::formats::hdf5::read_file(p, schema, sink, opt);
        });
  }

} // namespace exaStamp::atom_io

namespace YAML
{

  // schema: a map of dataset paths (see read_hdf5's documentation). A scalar names a preset (none yet).
  template <> struct convert<exaStamp::atom_io::Hdf5Schema>
  {
    static bool decode(const Node &node, exaStamp::atom_io::Hdf5Schema &h)
    {
      auto &s = h.s;
      s = {};
      if (node.IsScalar())
      {
        onika::lerr << "read_hdf5: unknown schema preset '" << node.as<std::string>() << "' (no preset defined yet)" << std::endl;
        return false;
      }
      if (!node.IsMap())
      {
        onika::lerr << "read_hdf5: schema must be a map of dataset paths" << std::endl;
        return false;
      }
      auto vec3 = [](const Node &v, aio::formats::hdf5::Vec3Path &p, const std::string &key)
      {
        if (v.IsScalar())
        {
          p.path = v.as<std::string>();
          return true;
        }
        if (v.IsSequence() && v.size() == 3)
        {
          for (std::size_t c = 0; c < 3; ++c)
            p.xyz[c] = v[c].as<std::string>();
          return true;
        }
        onika::lerr << "read_hdf5: schema." << key << " must be one (N,3) dataset path or a list of 3 (N) paths" << std::endl;
        return false;
      };
      for (const auto &kv : node)
      {
        const std::string key = kv.first.as<std::string>();
        const Node &v = kv.second;
        if (key == "positions" || key == "velocities")
        {
          if (!vec3(v, key == "positions" ? s.positions : s.velocities, key))
            return false;
        }
        else if (key == "ids" || key == "types" || key == "cell" || key == "origin" || key == "pbc")
        {
          std::string &dst = key == "ids" ? s.ids : key == "types" ? s.types : key == "cell" ? s.cell : key == "origin" ? s.origin : s.pbc;
          dst = v.as<std::string>();
        }
        else if (key == "fields")
        {
          if (!v.IsMap())
          {
            onika::lerr << "read_hdf5: schema.fields must be a map { name: dataset path }" << std::endl;
            return false;
          }
          for (const auto &f : v)
            s.properties.push_back({f.first.as<std::string>(), f.second.as<std::string>()});
        }
        else
        {
          onika::lerr << "read_hdf5: unknown schema key '" << key << "' (positions, velocities, ids, types, cell, origin, pbc, fields)" << std::endl;
          return false;
        }
      }
      return true;
    }
  };

} // namespace YAML
