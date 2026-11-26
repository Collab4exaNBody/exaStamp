/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/

// HDF5 atom reader (read_hdf5, aio ADR 0016).
//
// Kept apart from atom_file_reader.h on purpose: the text readers are the
// reference and their code paths are not touched. This header only *uses*
// its helpers (BoxTransform, TypeMap, GridFiller, ...).
//
// GridBindSink stages atoms in contiguous arrays (r[3n] or x[n] y[n] z[n],
// type[n], ids[n], properties), so aio reads each dataset with one H5Dread
// straight into them (aio ADR 0015). end() then does the per-atom work of
// GridStreamSink::atom() (box transform, wrap, cell, type map) and stages the
// atoms in GridFiller, whose commit() fills the grid as for the text readers.
//
// Flows (read_hdf5_into_grid):
//   - MPIO (parallel HDF5, mpio: true): every rank reads its slab of rows and
//     fills its own grid (full-domain geometry, as the other ranks get after a
//     rank-0 read). The deck's load balancing / migration redistributes the
//     atoms afterwards, as after any read.
//   - serial: rank 0 reads, then the rank-0 flow of read_atom_file. Used when
//     HDF5 has no parallel I/O (with a warning if several ranks run) or with
//     mpio: false.
//   - stream_sink: true (debug / benchmark): read_atom_file with GridStreamSink.
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

namespace exaStamp::atom_io {

// The `schema` slot of read_hdf5: dataset paths of the aio HDF5 reader.
struct Hdf5Schema {
  aio::formats::hdf5::Schema s;
};

inline std::ostream& operator<<(std::ostream& out, const aio::formats::hdf5::Vec3Path& v) {
  if (!v.path.empty())
    return out << v.path;
  return out << "[" << v.xyz[0] << ", " << v.xyz[1] << ", " << v.xyz[2] << "]";
}

inline std::ostream& operator<<(std::ostream& out, const Hdf5Schema& h) {
  const auto& s = h.s;
  out << "{ positions: " << s.positions;
  if (!s.velocities.empty())
    out << ", velocities: " << s.velocities;
  for (const auto& [k, v] : {std::pair<const char*, const std::string*>{"ids", &s.ids},
                             {"types", &s.types}, {"cell", &s.cell}, {"origin", &s.origin}, {"pbc", &s.pbc}})
    if (!v->empty())
      out << ", " << k << ": " << *v;
  if (!s.properties.empty()) {
    out << ", fields: {";
    for (std::size_t k = 0; k < s.properties.size(); ++k)
      out << (k ? ", " : " ") << s.properties[k].name << ": " << s.properties[k].path;
    out << " }";
  }
  return out << " }";
}

// GridBindSink: contiguous staging for direct HDF5 reads, then GridFiller.
template <class GridT>
class GridBindSink {
public:
  GridBindSink(Domain& domain, GridT& grid, const ReaderSettings& settings, bool interleaved_pos,
               bool interleaved_vel, bool ids_from_file)
      : m_domain(domain), m_grid(grid), m_species(settings.species), m_bounds_mode(settings.bounds_mode),
        m_read_velocities(settings.read_velocities), m_file_has_pbc(settings.file_has_pbc),
        m_props(settings.properties), m_reserved(reserved_field_names(typename GridT::Fields{})),
        m_interleaved_pos(interleaved_pos), m_interleaved_vel(interleaved_vel), m_ids_from_file(ids_from_file) {}

  // shortcut: begin(), bind_property() and cell_of() are copies of GridStreamSink's (private
  // there; the text readers stay untouched). Merge them once read_hdf5 is validated.
  aio::Status begin(const aio::atoms::FrameInfo& info) {
    m_natoms = info.natoms;
    m_natoms_global = info.natoms_global ? info.natoms_global : info.natoms;
    m_first_row = info.first_row;
    if (m_natoms > std::numeric_limits<std::uint32_t>::max())
      return aio::make_error(aio::Errc::Unsupported, 0, 0, "more than 2^32 atoms on one rank");
    const auto& h = info.cell.h;
    const double det = h(0, 0) * (h(1, 1) * h(2, 2) - h(1, 2) * h(2, 1)) -
                       h(0, 1) * (h(1, 0) * h(2, 2) - h(1, 2) * h(2, 0)) +
                       h(0, 2) * (h(1, 0) * h(2, 1) - h(1, 1) * h(2, 0));
    if (!(std::abs(det) > 0.0))
      return aio::make_error(aio::Errc::Unsupported, 0, 0, "the file defines no simulation cell (schema: cell)");
    if (m_file_has_pbc) {
      const bool dp[3] = {m_domain.periodic_boundary_x(), m_domain.periodic_boundary_y(), m_domain.periodic_boundary_z()};
      if (dp[0] != info.cell.pbc[0] || dp[1] != info.cell.pbc[1] || dp[2] != info.cell.pbc[2])
        lout << "Warning: file pbc = " << std::boolalpha << info.cell.pbc[0] << "," << info.cell.pbc[1] << ","
             << info.cell.pbc[2] << " differs from the domain's periodic flags; the domain's are used" << std::endl;
    }
    m_box.init(info.cell, m_domain, m_bounds_mode);
    m_grid.set_offset(IJK{0, 0, 0});
    m_grid.set_origin(m_domain.bounds().bmin);
    m_grid.set_cell_size(m_domain.cell_size());
    m_grid.set_dimension(m_domain.grid_dimension());
    m_dims = m_domain.grid_dimension();
    m_grid_origin = m_domain.origin();
    m_cell = m_domain.cell_size();
    m_dims_d[0] = static_cast<double>(m_dims.i);
    m_dims_d[1] = static_cast<double>(m_dims.j);
    m_dims_d[2] = static_cast<double>(m_dims.k);
    m_periodic[0] = m_domain.periodic_boundary_x();
    m_periodic[1] = m_domain.periodic_boundary_y();
    m_periodic[2] = m_domain.periodic_boundary_z();
    m_filler.reset(m_natoms);
    return {};
  }

  // Contiguous staging, allocated for the fields the schema has (aio binds those only).
  // Not zero-filled: aio's direct reads write every value.
  aio::atoms::FieldBinding bind(aio::atoms::Field f) {
    using F = aio::atoms::Field;
    using B = aio::atoms::FieldBinding;
    switch (f) {
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
      return vec_binding(m_vel, m_interleaved_vel, static_cast<int>(f) - static_cast<int>(F::VX));
    default:
      return {};
    }
  }

  aio::atoms::FieldBinding bind_property(std::string_view name, std::size_t ncomp, std::size_t comp,
                                         aio::atoms::ScalarType) {
    if (!m_props.selects(name))
      return {};
    const bool named = m_props.named(name);
    if (comp == 0) {
      m_seen.emplace_back(name);
      std::string why;
      for (const auto& [fname, kind] : GridFiller<GridT>::field_names(std::string(name), ncomp)) {
        if (fname.size() > kGenericNameMax)
          why = "generic field name '" + fname + "' longer than " + std::to_string(kGenericNameMax) + " characters";
        else if (std::find(m_reserved.begin(), m_reserved.end(), fname) != m_reserved.end())
          why = "'" + fname + "' is a field of the grid";
      }
      if (!why.empty()) {
        if (named)
          m_prop_error += (m_prop_error.empty() ? "" : "; ") + ("property '" + std::string(name) + "': " + why);
        else
          m_skipped.push_back(std::string(name) + " (" + why + ")");
        m_rejected.emplace_back(name);
        return {};
      }
    } else if (std::find(m_rejected.begin(), m_rejected.end(), name) != m_rejected.end()) {
      return {};
    }
    double* d = m_filler.property_storage(name, ncomp);
    if (!d)
      return {};
    return aio::atoms::FieldBinding::of(d, static_cast<std::ptrdiff_t>(ncomp), static_cast<std::ptrdiff_t>(comp));
  }

  // All datasets read: the per-atom work of GridStreamSink::atom(), then GridFiller staging.
  aio::Status end() {
    if (!m_prop_error.empty()) {
      m_msg = m_prop_error;
      return aio::make_error(aio::Errc::Sink, 0, 0, "invalid generic property selection");
    }
    std::string missing;
    for (const std::string& n : m_props.names)
      if (std::find(m_seen.begin(), m_seen.end(), n) == m_seen.end())
        missing += (missing.empty() ? "" : ", ") + n;
    if (!missing.empty()) {
      m_msg = "not in the schema or read as a grid field: " + missing;
      return aio::make_error(aio::Errc::Sink, 0, 0, "selected properties not found");
    }

    // Types: t -> species t - 1 (no labels in HDF5). Without species, the largest type sets the count.
    std::size_t ntypes = m_species ? m_species->size() : 1;
    if (!m_species && m_type)
      for (std::size_t i = 0; i < m_natoms; ++i)
        ntypes = std::max<std::size_t>(ntypes, static_cast<std::size_t>(std::max<std::int32_t>(m_type[i], 1)));
    if (!m_types.freeze(ntypes, m_species, m_msg))
      return aio::make_error(aio::Errc::Sink, 0, 0, "atom types");

    for (std::size_t i = 0; i < m_natoms; ++i) {
      const std::int64_t t = m_type ? m_type[i] : 1; // no types in the schema: every atom is type 1
      if (t < 1 || static_cast<std::size_t>(t) > m_types.size()) {
        m_msg = "atom type " + std::to_string(t) + " at row " + std::to_string(m_first_row + i) + " outside [1, " +
                std::to_string(m_types.size()) + "]";
        return aio::make_error(aio::Errc::Sink, 0, 0, "atom type out of range");
      }
      if (m_ids_from_file && (m_ids[i] < 1 || static_cast<std::uint64_t>(m_ids[i]) > m_natoms_global)) {
        m_msg = "atom id " + std::to_string(m_ids[i]) + " at row " + std::to_string(m_first_row + i) +
                " outside [1, " + std::to_string(m_natoms_global) + "]";
        return aio::make_error(aio::Errc::IdOutOfRange, 0, 0, "atom id out of range");
      }
      if (m_ids_from_file)
        m_ids_are_rows = m_ids_are_rows && static_cast<std::uint64_t>(m_ids[i]) == m_first_row + i + 1;
      const Vec3d p = component(m_pos, m_interleaved_pos, i);
      // grid frame, relative to the domain origin (as GridStreamSink::atom)
      Vec3d r = m_box.inv_xform * (p - m_box.origin);
      if (m_periodic[0])
        r.x = wrap(r.x, m_box.grid_len.x);
      if (m_periodic[1])
        r.y = wrap(r.y, m_box.grid_len.y);
      if (m_periodic[2])
        r.z = wrap(r.z, m_box.grid_len.z);
      r = r + m_box.grid_origin;
      const IJK loc = cell_of(r);
      if (loc.i < 0 || loc.j < 0 || loc.k < 0 || loc.i >= m_dims.i || loc.j >= m_dims.j || loc.k >= m_dims.k) {
        m_msg = "atom at row " + std::to_string(m_first_row + i) + " outside the domain on a non-periodic axis";
        return aio::make_error(aio::Errc::Sink, 0, 0, "atom outside the domain");
      }
      m_filler.stage(i, static_cast<std::uint32_t>(grid_ijk_to_index(m_dims, loc)), r,
                     m_types[static_cast<std::size_t>(t - 1)]);
      if (m_vel[0]) {
        const Vec3d v = component(m_vel, m_interleaved_vel, i);
        m_filler.stage_velocity(i, v.x, v.y, v.z);
      }
    }
    m_count = m_natoms;
    for (auto& a : m_pos)
      a.reset();
    for (auto& a : m_vel)
      a.reset();
    m_type.reset();
    return {};
  }

  // Puts the staged atoms into the grid. GridFiller sets id = local index i: rewritten to the
  // global id (file id - 1, or first_row + i with remap_atom).
  // shortcut: ids are checked in [1, N] but not for uniqueness (it would need communication in
  // MPIO); read_lammps_data / read_extxyz check it, read_hdf5 trusts the file.
  void commit() {
    m_filler.commit(m_grid);
    const bool by_row = !m_ids_from_file || m_ids_are_rows; // id = first_row + i: no lookup
    if (by_row && m_first_row == 0)
      return; // ids are already the indices
    auto cells = m_grid.cells();
    const std::size_t ncells = m_grid.number_of_cells();
    for (std::size_t c = 0; c < ncells; ++c) {
      const std::size_t n = cells[c].size();
      auto* id = cells[c][field::id];
      for (std::size_t p = 0; p < n; ++p) {
        const std::size_t i = static_cast<std::size_t>(id[p]);
        id[p] = by_row ? m_first_row + i : static_cast<std::uint64_t>(m_ids[i] - 1);
      }
    }
    m_ids.reset();
  }

  std::vector<GenericFieldInfo> fill_properties() { return m_filler.fill_properties(m_grid); }
  const std::vector<std::string>& skipped_properties() const { return m_skipped; }
  const std::string& message() const { return m_msg; }
  std::size_t count() const { return m_count; }
  const BoxTransform& box() const { return m_box; }

private:
  template <typename T>
  T* alloc(std::unique_ptr<T[]>& p, std::size_t n) {
    if (!p)
      p.reset(new T[n > 0 ? n : 1]);
    return p.get();
  }

  // Component c of a vector field: r[3n] interleaved (one (N,3) dataset) or x[n] y[n] z[n].
  aio::atoms::FieldBinding vec_binding(std::unique_ptr<double[]> (&a)[3], bool interleaved, int c) {
    if (interleaved)
      return aio::atoms::FieldBinding::of(alloc(a[0], 3 * m_natoms), 3, c);
    return aio::atoms::FieldBinding::of(alloc(a[c], m_natoms));
  }
  static Vec3d component(const std::unique_ptr<double[]> (&a)[3], bool interleaved, std::size_t i) {
    if (interleaved)
      return Vec3d{a[0][3 * i], a[0][3 * i + 1], a[0][3 * i + 2]};
    return Vec3d{a[0][i], a[1][i], a[2][i]};
  }

  IJK cell_of(Vec3d& r) const {
    const double qx = (r.x - m_grid_origin.x) / m_cell;
    const double qy = (r.y - m_grid_origin.y) / m_cell;
    const double qz = (r.z - m_grid_origin.z) / m_cell;
    if (qx >= 0.0 && qy >= 0.0 && qz >= 0.0 && qx < m_dims_d[0] && qy < m_dims_d[1] && qz < m_dims_d[2])
      return IJK{static_cast<ssize_t>(qx), static_cast<ssize_t>(qy), static_cast<ssize_t>(qz)};
    return domain_periodic_location(m_domain, r);
  }

  Domain& m_domain;
  GridT& m_grid;
  const ParticleSpecies* m_species;
  ReadBoundsSelectionMode m_bounds_mode;
  bool m_read_velocities = true;
  bool m_file_has_pbc = false;
  PropertySelection m_props;
  std::vector<std::string> m_reserved;
  std::vector<std::string> m_seen, m_rejected, m_skipped;
  std::string m_prop_error;
  bool m_interleaved_pos, m_interleaved_vel, m_ids_from_file;
  bool m_ids_are_rows = true; // file ids == row + 1 (sorted files): no id lookup in commit()
  BoxTransform m_box;
  TypeMap m_types;
  GridFiller<GridT> m_filler;
  std::unique_ptr<double[]> m_pos[3], m_vel[3];
  std::unique_ptr<std::int32_t[]> m_type;
  std::unique_ptr<std::int64_t[]> m_ids;
  IJK m_dims{};
  Vec3d m_grid_origin{};
  double m_cell = 1.0;
  double m_dims_d[3] = {0, 0, 0};
  bool m_periodic[3] = {true, true, true};
  std::size_t m_natoms = 0, m_natoms_global = 0, m_first_row = 0, m_count = 0;
  std::string m_msg;
};

// The domain fields every rank must share after an MPIO read.
inline bool same_domain(const Domain& a, const Domain& b) {
  return a.bounds().bmin == b.bounds().bmin && a.bounds().bmax == b.bounds().bmax &&
         a.cell_size() == b.cell_size() && a.grid_dimension() == b.grid_dimension() && a.xform() == b.xform();
}

// read_hdf5 flows (see the header comment). `mpio`: every rank reads its slab (needs
// AIO_HDF5_HAVE_MPIO); otherwise rank 0 reads everything.
template <class GridT>
void read_hdf5_into_grid(MPI_Comm comm, Domain& domain, GridT& grid, const ReaderSettings& settings,
                         const std::string& path, const aio::formats::hdf5::Schema& schema, bool remap_atom, bool mpio) {
  const char* op_name = "read_hdf5";
  lout << "======== " << op_name << " ========" << std::endl;
  int rank = 0, nranks = 1;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nranks);
  if (grid.number_of_particles() != 0) {
    lerr << op_name << ": system already defined, skipping" << std::endl;
    return;
  }
#if !AIO_HDF5_HAVE_MPIO
  if (mpio && nranks > 1)
    lout << "Warning: " << op_name
         << ": HDF5 built without parallel I/O (H5_HAVE_PARALLEL): reading on rank 0 only, then the atoms are "
            "distributed as after any rank-0 read"
         << std::endl;
  mpio = false;
#endif

  aio::formats::hdf5::Options opt; // row order: ids are handled by the sink
  const bool reads = mpio || rank == 0;
  ReadOutcome out;
  std::vector<GenericFieldInfo> generic;
  std::size_t count = 0;
  double ms_read = 0, ms_grid = 0;
  {
    GridBindSink<GridT> sink(domain, grid, settings, !schema.positions.path.empty(),
                             !schema.velocities.path.empty(), !remap_atom);
    aio::Status st;
    using clock = std::chrono::steady_clock;
    MPI_Barrier(comm); // times the read only (MPI's first collectives may set up connections)
    const auto t0 = clock::now();
    if (reads) {
#if AIO_HDF5_HAVE_MPIO
      st = mpio ? aio::formats::hdf5::read_file_mpio(comm, path.c_str(), schema, sink, opt)
                : aio::formats::hdf5::read_file(path.c_str(), schema, sink, opt);
#else
      st = aio::formats::hdf5::read_file(path.c_str(), schema, sink, opt);
#endif
    }
    const auto t1 = clock::now();
    if (reads && !st.ok()) {
      out.ok = 0;
      std::snprintf(out.msg, sizeof out.msg, "%s: %s%s%s", path.c_str(), st.error().msg,
                    sink.message().empty() ? "" : " -- ", sink.message().c_str());
    }
    // Every rank learns the outcome (MPIO: the first failing rank's message).
    if (mpio) {
      int bad = out.ok ? nranks : rank;
      MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MIN, comm);
      if (bad < nranks)
        MPI_Bcast(&out, sizeof(ReadOutcome), MPI_BYTE, bad, comm);
    } else {
      MPI_Bcast(&out, sizeof(ReadOutcome), MPI_BYTE, 0, comm);
    }
    if (!out.ok)
      fatal_error() << op_name << ": " << out.msg << std::endl;
    const auto t2 = clock::now();
    if (reads) {
      sink.commit();
      grid.rebuild_particle_offsets(); // flat arrays are indexed by these offsets
      generic = sink.fill_properties();
      count = sink.count();
      if (rank == 0)
        for (const std::string& sk : sink.skipped_properties())
          lout << "Skipped property = " << sk << std::endl;
    }
    ms_read = std::chrono::duration<double, std::milli>(t1 - t0).count();
    ms_grid = std::chrono::duration<double, std::milli>(clock::now() - t2).count(); // commit + offsets + fields
  }

  // Domain: rank 0's for everyone. MPIO ranks computed it from the same file cell: it must match.
  Domain mine = domain;
  MPI_Bcast(&domain, sizeof(Domain), MPI_BYTE, 0, comm);
  if (mpio) {
    int differs = same_domain(mine, domain) ? 0 : 1;
    MPI_Allreduce(MPI_IN_PLACE, &differs, 1, MPI_INT, MPI_MAX, comm);
    if (differs)
      fatal_error() << op_name << ": ranks computed different domains from the file cell" << std::endl;
  } else if (rank != 0) {
    grid.set_offset(IJK{0, 0, 0});
    grid.set_origin(domain.bounds().bmin);
    grid.set_cell_size(domain.cell_size());
    grid.set_dimension(domain.grid_dimension());
  }
  grid.rebuild_particle_offsets();

  // Same generic fields on every rank (MPIO: each rank created them in fill_properties).
  int ngeneric = static_cast<int>(generic.size());
  MPI_Bcast(&ngeneric, 1, MPI_INT, 0, comm);
  if (ngeneric > 0 && !mpio) {
    generic.resize(static_cast<std::size_t>(ngeneric));
    MPI_Bcast(generic.data(), static_cast<int>(ngeneric * sizeof(GenericFieldInfo)), MPI_BYTE, 0, comm);
    if (rank != 0)
      for (const GenericFieldInfo& f : generic)
        generic_field_data(grid, f);
  }

  unsigned long long total = count;
  MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
  MPI_Allreduce(MPI_IN_PLACE, &ms_read, 1, MPI_DOUBLE, MPI_MAX, comm);
  MPI_Allreduce(MPI_IN_PLACE, &ms_grid, 1, MPI_DOUBLE, MPI_MAX, comm);
  lout << onika::format_string("Read %llu atoms in %.1f ms (read %.1f ms, grid insert %.1f ms; %s on %d rank%s)",
                               total, ms_read + ms_grid, ms_read, ms_grid, mpio ? "MPIO" : "rank 0",
                               mpio ? nranks : 1, (mpio && nranks > 1) ? "s" : "")
       << std::endl;
  lout << "Domain XForm     = " << domain.xform() << std::endl;
  lout << "Domain bounds    = " << domain.bounds() << std::endl;
  lout << "Cell size        = " << domain.cell_size() << std::endl;
  lout << "Grid dimensions  = " << domain.grid_dimension() << " (" << grid_cell_count(domain.grid_dimension())
       << " cells)" << std::endl;
  lout << "Domain periodic  = " << std::boolalpha << domain.periodic_boundary_x() << "," << domain.periodic_boundary_y()
       << "," << domain.periodic_boundary_z() << std::endl;
  lout << "Domain expandable = " << std::boolalpha << domain.expandable() << std::endl;
  for (const GenericFieldInfo& f : generic)
    lout << "Generic field    = " << f.name << " (" << generic_kind_name(f.kind) << ")" << std::endl;
  if (!check_domain(domain)) {
    lout << "domain = " << domain << std::endl;
    fatal_error() << op_name << ": invalid domain configuration" << std::endl;
  }
  lout << "===========================================" << std::endl << std::endl;
}

} // namespace exaStamp::atom_io

namespace YAML {

// schema: a map of dataset paths (see read_hdf5's documentation). A scalar names a preset (none yet).
template <>
struct convert<exaStamp::atom_io::Hdf5Schema> {
  static bool decode(const Node& node, exaStamp::atom_io::Hdf5Schema& h) {
    auto& s = h.s;
    s = {};
    if (node.IsScalar()) {
      onika::lerr << "read_hdf5: unknown schema preset '" << node.as<std::string>() << "' (no preset defined yet)"
                  << std::endl;
      return false;
    }
    if (!node.IsMap()) {
      onika::lerr << "read_hdf5: schema must be a map of dataset paths" << std::endl;
      return false;
    }
    auto vec3 = [](const Node& v, aio::formats::hdf5::Vec3Path& p, const std::string& key) {
      if (v.IsScalar()) {
        p.path = v.as<std::string>();
        return true;
      }
      if (v.IsSequence() && v.size() == 3) {
        for (std::size_t c = 0; c < 3; ++c)
          p.xyz[c] = v[c].as<std::string>();
        return true;
      }
      onika::lerr << "read_hdf5: schema." << key << " must be one (N,3) dataset path or a list of 3 (N) paths"
                  << std::endl;
      return false;
    };
    for (const auto& kv : node) {
      const std::string key = kv.first.as<std::string>();
      const Node& v = kv.second;
      if (key == "positions" || key == "velocities") {
        if (!vec3(v, key == "positions" ? s.positions : s.velocities, key))
          return false;
      } else if (key == "ids" || key == "types" || key == "cell" || key == "origin" || key == "pbc") {
        std::string& dst = key == "ids" ? s.ids : key == "types" ? s.types : key == "cell" ? s.cell
                                                                  : key == "origin" ? s.origin : s.pbc;
        dst = v.as<std::string>();
      } else if (key == "fields") {
        if (!v.IsMap()) {
          onika::lerr << "read_hdf5: schema.fields must be a map { name: dataset path }" << std::endl;
          return false;
        }
        for (const auto& f : v)
          s.properties.push_back({f.first.as<std::string>(), f.second.as<std::string>()});
      } else {
        onika::lerr << "read_hdf5: unknown schema key '" << key
                    << "' (positions, velocities, ids, types, cell, origin, pbc, fields)" << std::endl;
        return false;
      }
    }
    return true;
  }
};

} // namespace YAML
