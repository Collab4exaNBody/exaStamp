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

// Shared code of the atom file readers built on aio (exaStamp/io/aio.hpp).
//
// One operator per file format (read_lammps_data, ...). Each one declares
// EXASTAMP_ATOM_READER_SLOTS() and calls read_atom_file() with a function
// that reads its format into a sink. Compression is detected by aio from
// the file content (aio::text::with_source).
//
// Flow (read_into_grid, shared with read_hdf5; text: rank 0 reads, as the
// previous read_external_atom_format_v2):
//   1. header parsed -> box, domain and grid geometry set (begin(), GridPlacer)
//   2. each atom line: grid-frame position and cell computed (GridPlacer::place),
//      staged by index (GridFiller); a Velocities section fills the same staging slots
//   3. GridFiller::commit: counting sort by cell, one resize per cell,
//      sequential writes (much faster than one push_back per atom, ADR 0007)
//   4. generic per-atom properties (formats with named columns: extxyz,
//      dump), if selected: staged by index while reading, then copied into
//      grid flat arrays (mk_generic_real/vec3/mat3) in grid order (ADR 0011)
//   5. {status, domain, generic field list} broadcast: every rank stops on a
//      read error; other ranks set their grid geometry and create the same
//      (empty) generic fields; offsets rebuilt; domain checked
//
// Generic fields are valid right after reading only: exanb does not move flat
// arrays with the particles, so they no longer match once particles move or
// migrate (first parallel_update_particles). Use them in setup_system, e.g.
// read -> project on a grid -> dump.
#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include <mpi.h>

#include <onika/log.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/print_utils.h>

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>

#include <exaStamp/io/aio.hpp>
#include <exaStamp/particle_species/particle_specie.h>

// Slots common to every atom file reader.
#define EXASTAMP_ATOM_READER_SLOTS()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           \
  ADD_SLOT(MPI_Comm, mpi, INPUT, MPI_COMM_WORLD);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              \
  ADD_SLOT(Domain, domain, INPUT_OUTPUT);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      \
  ADD_SLOT(GridT, grid, INPUT_OUTPUT);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         \
  ADD_SLOT(ParticleSpecies, species, INPUT, OPTIONAL, DocString{"If set, file type labels are matched to species names; otherwise type t -> species t-1"});                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    \
  ADD_SLOT(ReadBoundsSelectionMode, bounds_mode, INPUT, ReadBoundsSelectionMode::COMPUTED_BOUNDS);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             \
  ADD_SLOT(std::string, file, INPUT, REQUIRED, DocString{"Input file (.gz/.bz2/.xz/.zst detected from content)"});                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             \
  ADD_SLOT(bool, remap_atom, INPUT, false, DocString{"true: particle ids follow the file order; false: id = file id - 1"});                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    \
  ADD_SLOT(bool, read_velocities, INPUT, true, DocString{"false: velocities in the file are not parsed; every velocity is set to zero"})

// Values of the common slots, for read_atom_file (use inside execute()).
#define EXASTAMP_ATOM_READER_SETTINGS()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        \
  ::exaStamp::atom_io::ReaderSettings { species.has_value() ? species.get_pointer() : nullptr, *bounds_mode, *read_velocities }

namespace exaStamp::atom_io
{

  using namespace exanb;

  // Per-atom properties to store as generic fields: none (default), all ("*"),
  // or a list of names. Formats with named columns only (extxyz, dump).
  struct PropertySelection
  {
    bool all = false;
    std::vector<std::string> names;

    PropertySelection() = default;
    explicit PropertySelection(const std::vector<std::string> &list)
    {
      for (const std::string &s : list)
      {
        if (s == "*")
          all = true;
        else
          names.push_back(s);
      }
    }
    bool any() const { return all || !names.empty(); }
    bool named(std::string_view n) const { return std::find(names.begin(), names.end(), n) != names.end(); }
    bool selects(std::string_view n) const { return all || named(n); }
  };

  // ExaStamp takes id, type, x y z and (optionally) vx vy vz from a file into
  // the grid: charges and masses come from the species definitions. The sink
  // asks aio for these fields only, so every other column is skipped without
  // conversion, whatever the format, unless it is selected as a generic
  // property (`properties`).
  struct ReaderSettings
  {
    const ParticleSpecies *species = nullptr;
    ReadBoundsSelectionMode bounds_mode = ReadBoundsSelectionMode::COMPUTED_BOUNDS;
    bool read_velocities = true;
    PropertySelection properties{};
    bool file_has_pbc = false; // the format states periodicity (extxyz pbc=): warn if it differs from the deck
  };

  // A generic field created by the reader (POD: broadcast to every rank).
  enum class GenericKind : std::uint8_t
  {
    Real,
    Vec3,
    Mat3
  };
  inline constexpr std::size_t kGenericNameMax = 15; // exanb dynamic field names: 16 bytes with '\0'
  struct GenericFieldInfo
  {
    char name[kGenericNameMax + 1] = {};
    GenericKind kind = GenericKind::Real;
  };

  inline const char *generic_kind_name(GenericKind k) { return k == GenericKind::Real ? "scalar" : k == GenericKind::Vec3 ? "Vec3" : "Mat3"; }

  // Creates (or resizes) the flat array of `f` on this rank and returns its data.
  template <class GridT> inline void *generic_field_data(GridT &grid, const GenericFieldInfo &f)
  {
    switch (f.kind)
    {
    case GenericKind::Real:
      return grid.flat_array_data(field::mk_generic_real(f.name));
    case GenericKind::Vec3:
      return grid.flat_array_data(field::mk_generic_vec3(f.name));
    case GenericKind::Mat3:
      return grid.flat_array_data(field::mk_generic_mat3(f.name));
    }
    return nullptr;
  }

  // Names a generic field may not take: the grid's own fields (short names) and
  // the names exanb tools (write_xyz) derive from them, so no output column is
  // duplicated.
  template <class... ids> inline std::vector<std::string> reserved_field_names(FieldSet<ids...>) { return {std::string(onika::soatl::FieldId<ids>::short_name())..., "position", "velocity", "force", "processor_id", "pos", "vel"}; }

  // Periodic wrap into [0, len). Fast path without floor(): ExaStamp builds
  // without SSE4.1, so floor() is a libm call (ADR 0007).
  inline double wrap(double x, double len)
  {
    if (x >= 0.0 && x < len)
      return x;
    double r = x - len * std::floor(x / len);
    if (r < 0.0)
      r += len;
    if (r >= len)
      r -= len;
    return r;
  }

  // Box -> domain transform (ADR 0010; first a port of v2's init_from_context).
  // The file cell H is the physical system and is never changed. The deck's
  // domain only chooses how it is represented: the grid is the axis-aligned box
  // [bmin, bmin + grid_len) of the domain, and Xform maps it onto H:
  //   r_real - file origin = Xform * (r_grid - bmin).
  // The deck xform Xu sets the grid-space size of the box (inv(Xu) * H), hence
  // the number of cells; the final Xform is always fitted to H. The box origin
  // is kept relative (ADR 0007): the file origin goes to the domain origin.
  struct BoxTransform
  {
    Mat3d xform = make_identity_matrix();
    Mat3d inv_xform = make_identity_matrix();
    Vec3d origin{};      // file origin (xlo ylo zlo)
    Vec3d grid_origin{}; // domain bounds().bmin
    Vec3d grid_len{};    // grid-space size of the file cell
    bool uniform_scale = true;

    void init(const aio::atoms::Cell &cell, Domain &domain, ReadBoundsSelectionMode bounds_mode)
    {
      Mat3d H;
      H.m11 = cell.h(0, 0), H.m12 = cell.h(0, 1), H.m13 = cell.h(0, 2);
      H.m21 = cell.h(1, 0), H.m22 = cell.h(1, 1), H.m23 = cell.h(1, 2);
      H.m31 = cell.h(2, 0), H.m32 = cell.h(2, 1), H.m33 = cell.h(2, 2);
      origin = Vec3d{cell.origin.x, cell.origin.y, cell.origin.z};

      const Vec3d a{H.m11, H.m21, H.m31};
      const Vec3d b{H.m12, H.m22, H.m32};
      const Vec3d c{H.m13, H.m23, H.m33};
      const Vec3d boxlen{norm(a), norm(b), norm(c)};
      lout << "a = " << a << onika::format_string(" , norm = %10.5f ", boxlen.x) << std::endl;
      lout << "b = " << b << onika::format_string(" , norm = %10.5f ", boxlen.y) << std::endl;
      lout << "c = " << c << onika::format_string(" , norm = %10.5f ", boxlen.z) << std::endl;

      // Requested grid-space box: the file cell seen through the deck xform.
      Vec3d boxlen_g = boxlen;
      if (!domain.xform_is_identity())
      {
        const Mat3d G = domain.inv_xform() * H;
        boxlen_g = Vec3d{norm(Vec3d{G.m11, G.m21, G.m31}), norm(Vec3d{G.m12, G.m22, G.m32}), norm(Vec3d{G.m13, G.m23, G.m33})};
      }

      // Grid dims: from the domain if given, else grid-space length / cell size, at least 1.
      const double cs = domain.cell_size();
      const IJK dim = domain.grid_dimension();
      auto ncells = [cs](ssize_t given, double len) { return given > 0 ? given : std::max<ssize_t>(1, static_cast<ssize_t>(len / cs)); };
      const ssize_t nx = ncells(dim.i, boxlen_g.x), ny = ncells(dim.j, boxlen_g.y), nz = ncells(dim.k, boxlen_g.z);
      const AABB req{{0., 0., 0.}, {nx * cs, ny * cs, nz * cs}};

      // Bounds, dims and cell size as exanb settles them (bounds_mode, periodic
      // fit). The xform it may set is replaced below: only the grid box counts.
      compute_domain_bounds(domain, bounds_mode, 0.0, req, req, true);
      const AABB B = domain.bounds();
      grid_origin = B.bmin;
      grid_len = B.bmax - B.bmin;

      // On a non-periodic axis exanb may add one cell layer beyond the bounds
      // (for atoms on the upper face): extend the bounds to the grid so the
      // domain passes check_domain. The file cell still maps onto grid_len.
      const IJK d = domain.grid_dimension();
      const double cs_final = domain.cell_size();
      Vec3d bmax = B.bmax;
      if (!domain.periodic_boundary_x() && d.i * cs_final > grid_len.x)
        bmax.x = B.bmin.x + d.i * cs_final;
      if (!domain.periodic_boundary_y() && d.j * cs_final > grid_len.y)
        bmax.y = B.bmin.y + d.j * cs_final;
      if (!domain.periodic_boundary_z() && d.k * cs_final > grid_len.z)
        bmax.z = B.bmin.z + d.k * cs_final;
      domain.set_bounds(AABB{B.bmin, bmax});

      const Mat3d D1 = diag_matrix(boxlen);
      const Mat3d D2 = diag_matrix(grid_len);
      xform = (H * inverse(D1)) * (D1 * inverse(D2));
      inv_xform = inverse(xform);
      uniform_scale = onika::math::is_uniform_scale(xform);
      domain.set_xform(xform);
    }
  };

  // File type (1..ntypes) -> ExaStamp species index. With species and type
  // labels: by name. Otherwise type - 1 (as v2).
  class TypeMap
  {
  public:
    void set_label(std::int64_t type, std::string_view label)
    {
      if (m_frozen)
      {
        m_late_label = true;
        return;
      }
      if (type < 1)
        return;
      if (m_labels.size() < static_cast<std::size_t>(type))
        m_labels.resize(static_cast<std::size_t>(type));
      m_labels[static_cast<std::size_t>(type - 1)] = std::string(label);
    }

    // Called once before the first atom. Returns false with a message on error.
    bool freeze(std::size_t ntypes, const ParticleSpecies *species, std::string &msg)
    {
      m_frozen = true;
      m_map.assign(ntypes, 0);
      for (std::size_t t = 0; t < ntypes; ++t)
      {
        std::size_t s = t;
        const std::string *label = t < m_labels.size() && !m_labels[t].empty() ? &m_labels[t] : nullptr;
        if (species && label)
        {
          s = species->size();
          for (std::size_t k = 0; k < species->size(); ++k)
            if (species->at(k).m_name == *label)
              s = k;
          if (s == species->size())
          {
            msg = "type " + std::to_string(t + 1) + " label '" + *label + "' matches no species";
            return false;
          }
        }
        else if (species && s >= species->size())
        {
          msg = "file has " + std::to_string(ntypes) + " atom types but only " + std::to_string(species->size()) + " species are defined";
          return false;
        }
        if (s > 255)
        {
          msg = "type index " + std::to_string(s) + " does not fit the particle type field";
          return false;
        }
        m_map[t] = static_cast<std::uint8_t>(s);
      }
      return true;
    }

    bool late_label() const { return m_late_label; }
    std::size_t size() const { return m_map.size(); }
    std::uint8_t operator[](std::size_t t) const { return m_map[t]; }

  private:
    std::vector<std::string> m_labels;
    std::vector<std::uint8_t> m_map;
    bool m_frozen = false;
    bool m_late_label = false;
  };

  // The generic fields a property with `ncomp` components becomes: ncomp 1 ->
  // scalar, 3 -> Vec3, 9 -> Mat3 (row-major, as write_xyz writes it), otherwise
  // one scalar per component, `name_c`.
  inline std::vector<std::pair<std::string, GenericKind>> generic_field_names(const std::string &name, std::size_t ncomp)
  {
    if (ncomp == 1)
      return {{name, GenericKind::Real}};
    if (ncomp == 3)
      return {{name, GenericKind::Vec3}};
    if (ncomp == 9)
      return {{name, GenericKind::Mat3}};
    std::vector<std::pair<std::string, GenericKind>> v;
    for (std::size_t c = 0; c < ncomp; ++c)
      v.emplace_back(name + "_" + std::to_string(c), GenericKind::Real);
    return v;
  }

  // Grid geometry of `domain` (every rank, before particles are added).
  template <class GridT> inline void set_grid_geometry(GridT &grid, const Domain &domain)
  {
    grid.set_offset(IJK{0, 0, 0});
    grid.set_origin(domain.bounds().bmin);
    grid.set_cell_size(domain.cell_size());
    grid.set_dimension(domain.grid_dimension());
  }

  // GridFiller: generic, format-independent grid insertion (rank 0).
  // Readers stage atoms by index in any order (stage / stage_velocity), then
  // commit() places them: count per cell, one resize per cell with the grid's
  // allocator, sequential writes. Inside a cell, atoms are in index order.
  // Memory: 32 B per atom (staging) + 24 B per atom if velocities are given +
  // 4 B per atom (order), all freed by commit(); plus 8 B per component of each
  // selected generic property, kept with `order` until fill_properties(). The staging array is not
  // zero-filled: readers stage every index in [0, n) exactly once.
  template <class GridT> class GridFiller
  {
  public:
    using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_vx, field::_vy, field::_vz, field::_id, field::_type>;

    void reset(std::size_t natoms)
    {
      m_n = natoms;
      m_staged.reset(new Staged[natoms]); // uninitialized on purpose
      m_vel.reset();
      m_props.clear();
      m_order.clear();
    }
    std::size_t size() const { return m_n; }

    AIO_INLINE void stage(std::size_t index, std::uint32_t cell, const Vec3d &r, std::uint8_t type)
    {
      Staged &s = m_staged[index];
      s.r[0] = r.x, s.r[1] = r.y, s.r[2] = r.z;
      s.cell = cell;
      s.type = type;
    }
    AIO_INLINE void stage_velocity(std::size_t index, double vx, double vy, double vz)
    {
      double *v = &velocity_storage(true)[3 * index]; // first velocity: atoms without one keep 0
      v[0] = vx, v[1] = vy, v[2] = vz;
    }
    // Velocities as vx vy vz per atom. `zero`: false when the caller writes all 3n values
    // (a direct read), so the array is not filled twice.
    AIO_INLINE double *velocity_storage(bool zero)
    {
      if (AIO_UNLIKELY(!m_vel))
        m_vel.reset(zero ? new double[3 * m_n]() : new double[3 * m_n]);
      return m_vel.get();
    }

    // Staging of a generic property with `ncomp` components per atom (one array
    // per property name, zero-filled: atoms of a section without it keep 0).
    // The reader stores component c of atom i at [ncomp * i + c].
    double *property_storage(std::string_view name, std::size_t ncomp)
    {
      for (StagedProperty &p : m_props)
        if (p.name == name)
          return p.ncomp == ncomp ? p.data.get() : nullptr;
      m_props.push_back({std::string(name), ncomp, std::unique_ptr<double[]>(new double[ncomp * m_n]())});
      return m_props.back().data.get();
    }
    bool has_properties() const { return !m_props.empty(); }

    static void fields_of(const std::string &name, std::size_t ncomp, std::vector<GenericFieldInfo> &out)
    {
      for (const auto &[n, k] : generic_field_names(name, ncomp))
      {
        GenericFieldInfo f;
        std::strncpy(f.name, n.c_str(), kGenericNameMax); // checked by the sink: never truncated
        f.kind = k;
        out.push_back(f);
      }
    }

    // After commit() and rebuild_particle_offsets(): copies the staged
    // properties into the grid's flat arrays, in grid order (particle k of the
    // grid is staged atom order[k]: the grid was empty). Frees the staging.
    // Returns the generic fields created.
    std::vector<GenericFieldInfo> fill_properties(GridT &grid)
    {
      std::vector<GenericFieldInfo> created;
      const std::size_t n = m_order.size();
      if (grid.number_of_particles() != n)
        return created; // not reached: commit() put exactly these atoms in the empty grid
      for (const StagedProperty &p : m_props)
      {
        const std::size_t first = created.size();
        fields_of(p.name, p.ncomp, created);
        const double *d = p.data.get();
        const std::size_t nc = p.ncomp;
        if (nc == 1 || nc == 3 || nc == 9)
        {
          void *data = generic_field_data(grid, created[first]);
          if (nc == 1)
          {
            double *out = static_cast<double *>(data);
            for (std::size_t k = 0; k < n; ++k)
              out[k] = d[m_order[k]];
          }
          else if (nc == 3)
          {
            Vec3d *out = static_cast<Vec3d *>(data);
            for (std::size_t k = 0; k < n; ++k)
            {
              const double *v = d + 3 * static_cast<std::size_t>(m_order[k]);
              out[k] = Vec3d{v[0], v[1], v[2]};
            }
          }
          else
          {
            Mat3d *out = static_cast<Mat3d *>(data);
            for (std::size_t k = 0; k < n; ++k)
            {
              const double *v = d + 9 * static_cast<std::size_t>(m_order[k]);
              Mat3d &m = out[k];
              m.m11 = v[0], m.m12 = v[1], m.m13 = v[2];
              m.m21 = v[3], m.m22 = v[4], m.m23 = v[5];
              m.m31 = v[6], m.m32 = v[7], m.m33 = v[8];
            }
          }
        }
        else
        {
          for (std::size_t c = 0; c < nc; ++c)
          {
            double *out = static_cast<double *>(generic_field_data(grid, created[first + c]));
            for (std::size_t k = 0; k < n; ++k)
              out[k] = d[nc * static_cast<std::size_t>(m_order[k]) + c];
          }
        }
      }
      m_props.clear();
      m_order.clear();
      m_order.shrink_to_fit();
      return created;
    }

    // Appends the staged atoms to the grid cells and frees the staging. The id of
    // atom i is file_ids[i] - 1 (1-based ids from the file) or id_offset + i.
    void commit(GridT &grid, std::uint64_t id_offset = 0, const std::int64_t *file_ids = nullptr)
    {
      using CellTuple = std::remove_cv_t<std::remove_reference_t<decltype(grid.cells()[0][0])>>;
      const std::size_t ncells = grid.number_of_cells();
      const std::size_t n = m_n;
      const Staged *staged = m_staged.get();
      const double *vel = m_vel.get();

      std::vector<std::uint32_t> start(ncells + 1, 0); // counting sort by cell
      for (std::size_t i = 0; i < n; ++i)
        ++start[staged[i].cell + 1];
      for (std::size_t c = 0; c < ncells; ++c)
        start[c + 1] += start[c];
      std::vector<std::uint32_t> order(n);
      {
        std::vector<std::uint32_t> next(start.begin(), start.end() - 1);
        for (std::size_t i = 0; i < n; ++i)
          order[next[staged[i].cell]++] = static_cast<std::uint32_t>(i);
      }

      auto cells = grid.cells();
      const auto &alloc = grid.cell_allocator();
      for (std::size_t c = 0; c < ncells; ++c)
      {
        const std::size_t count = start[c + 1] - start[c];
        if (count == 0)
          continue;
        auto &cell = cells[c];
        const std::size_t base = cell.size(); // append (the grid is empty when reading)
        cell.resize(base + count, alloc);
        for (std::size_t j = 0; j < count; ++j)
        {
          const std::uint32_t i = order[start[c] + j];
          const Staged &s = staged[i];
          const double *v = vel ? vel + 3 * i : kZero;
          const std::uint64_t id = file_ids ? static_cast<std::uint64_t>(file_ids[i] - 1) : id_offset + i;
          cell.set_tuple(base + j, CellTuple(ParticleTupleIO(s.r[0], s.r[1], s.r[2], v[0], v[1], v[2], id, s.type)));
        }
      }
      m_staged.reset();
      m_vel.reset();
      m_n = 0;
      if (!m_props.empty()) // kept for fill_properties()
        m_order = std::move(order);
    }

  private:
    struct StagedProperty
    {
      std::string name;
      std::size_t ncomp;
      std::unique_ptr<double[]> data;
    };
    struct Staged
    { // 32 B; no default member initializers (no zero-fill)
      double r[3];
      std::uint32_t cell;
      std::uint8_t type;
    };
    static_assert(sizeof(Staged) == 32);
    static constexpr double kZero[3] = {0, 0, 0};
    std::size_t m_n = 0;
    std::unique_ptr<Staged[]> m_staged;
    std::unique_ptr<double[]> m_vel;
    std::vector<StagedProperty> m_props;
    std::vector<std::uint32_t> m_order;
  };

  // File position -> grid-frame position and grid cell, shared by every sink.
  // begin() fits the deck's domain to the file cell (BoxTransform, ADR 0010) and
  // caches the grid geometry; place() is the per-atom work.
  class GridPlacer
  {
  public:
    aio::Status begin(const aio::atoms::FrameInfo &info, Domain &domain, ReadBoundsSelectionMode bounds_mode, bool file_has_pbc, const char *no_cell_hint)
    {
      if (info.natoms > std::numeric_limits<std::uint32_t>::max())
        return aio::make_error(aio::Errc::Unsupported, 0, 0, "more than 2^32 atoms on one rank");
      const auto &h = info.cell.h;
      const double det = h(0, 0) * (h(1, 1) * h(2, 2) - h(1, 2) * h(2, 1)) - h(0, 1) * (h(1, 0) * h(2, 2) - h(1, 2) * h(2, 0)) + h(0, 2) * (h(1, 0) * h(2, 1) - h(1, 1) * h(2, 0));
      if (!(std::abs(det) > 0.0))
        return aio::make_error(aio::Errc::Unsupported, 0, 0, "the file defines no simulation cell (%s)", no_cell_hint);
      if (file_has_pbc)
      {
        const bool dp[3] = {domain.periodic_boundary_x(), domain.periodic_boundary_y(), domain.periodic_boundary_z()};
        if (dp[0] != info.cell.pbc[0] || dp[1] != info.cell.pbc[1] || dp[2] != info.cell.pbc[2])
          lout << "Warning: file pbc = " << std::boolalpha << info.cell.pbc[0] << "," << info.cell.pbc[1] << "," << info.cell.pbc[2] << " differs from the domain's periodic flags; the domain's are used" << std::endl;
      }
      m_box.init(info.cell, domain, bounds_mode);
      m_domain = &domain;
      m_dims = domain.grid_dimension();
      m_grid_origin = domain.origin();
      m_cell = domain.cell_size();
      m_dims_d[0] = static_cast<double>(m_dims.i);
      m_dims_d[1] = static_cast<double>(m_dims.j);
      m_dims_d[2] = static_cast<double>(m_dims.k);
      m_periodic[0] = domain.periodic_boundary_x();
      m_periodic[1] = domain.periodic_boundary_y();
      m_periodic[2] = domain.periodic_boundary_z();
      return {};
    }

    // Grid-frame position r of the file position p and its cell index. False:
    // outside the domain on a non-periodic axis.
    AIO_INLINE bool place(const Vec3d &p, Vec3d &r, std::uint32_t &cell) const
    {
      // grid frame, relative to the domain origin (0 + x is exact: bit-identical when bmin = 0)
      r = m_box.inv_xform * (p - m_box.origin);
      if (m_periodic[0])
        r.x = wrap(r.x, m_box.grid_len.x);
      if (m_periodic[1])
        r.y = wrap(r.y, m_box.grid_len.y);
      if (m_periodic[2])
        r.z = wrap(r.z, m_box.grid_len.z);
      r = r + m_box.grid_origin;
      // Fast path: inside the grid, truncation == floor and domain_periodic_location()
      // would not change r (same arithmetic as make_ijk((r - origin) / cell_size)).
      const double qx = (r.x - m_grid_origin.x) / m_cell;
      const double qy = (r.y - m_grid_origin.y) / m_cell;
      const double qz = (r.z - m_grid_origin.z) / m_cell;
      IJK loc;
      if (AIO_LIKELY(qx >= 0.0 && qy >= 0.0 && qz >= 0.0 && qx < m_dims_d[0] && qy < m_dims_d[1] && qz < m_dims_d[2]))
        loc = IJK{static_cast<ssize_t>(qx), static_cast<ssize_t>(qy), static_cast<ssize_t>(qz)};
      else
      {
        loc = domain_periodic_location(*m_domain, r);
        if (loc.i < 0 || loc.j < 0 || loc.k < 0 || loc.i >= m_dims.i || loc.j >= m_dims.j || loc.k >= m_dims.k)
          return false;
      }
      cell = static_cast<std::uint32_t>(grid_ijk_to_index(m_dims, loc));
      return true;
    }

    const BoxTransform &box() const { return m_box; }

  private:
    BoxTransform m_box;
    const Domain *m_domain = nullptr;
    IJK m_dims{};
    Vec3d m_grid_origin{};
    double m_cell = 1.0;
    double m_dims_d[3] = {0, 0, 0};
    bool m_periodic[3] = {true, true, true};
  };

  // Generic per-atom properties without a grid slot (ADR 0011), shared by every
  // sink: selected ones are bound to GridFiller staging (then stored as generic
  // fields); the others stay unbound, so aio skips them without parsing.
  class PropertyFilter
  {
  public:
    PropertyFilter(PropertySelection sel, std::vector<std::string> reserved) : m_sel(std::move(sel)), m_reserved(std::move(reserved)) {}

    // aio's bind_property hook: called once per component, before the atoms.
    template <class FillerT> aio::atoms::FieldBinding bind(FillerT &filler, std::string_view name, std::size_t ncomp, std::size_t comp)
    {
      if (!m_sel.selects(name))
        return {};
      if (comp == 0)
      {
        m_seen.emplace_back(name);
        std::string why;
        for (const auto &[fname, kind] : generic_field_names(std::string(name), ncomp))
        {
          if (fname.size() > kGenericNameMax)
            why = "generic field name '" + fname + "' longer than " + std::to_string(kGenericNameMax) + " characters";
          else if (std::find(m_reserved.begin(), m_reserved.end(), fname) != m_reserved.end())
            why = "'" + fname + "' is a field of the grid";
        }
        if (!why.empty())
        {
          if (m_sel.named(name))
            m_error += (m_error.empty() ? "" : "; ") + ("property '" + std::string(name) + "': " + why);
          else
            m_skipped.push_back(std::string(name) + " (" + why + ")");
          m_rejected.emplace_back(name);
          return {};
        }
      }
      else if (std::find(m_rejected.begin(), m_rejected.end(), name) != m_rejected.end())
      {
        return {};
      }
      double *d = filler.property_storage(name, ncomp);
      if (!d)
        return {};
      return aio::atoms::FieldBinding::of(d, static_cast<std::ptrdiff_t>(ncomp), static_cast<std::ptrdiff_t>(comp));
    }

    // At the end of the read: the selection is valid and every named property was found.
    aio::Status check(std::string &msg, const char *missing_hint) const
    {
      if (!m_error.empty())
      {
        msg = m_error;
        return aio::make_error(aio::Errc::Sink, 0, 0, "invalid generic property selection");
      }
      std::string missing;
      for (const std::string &n : m_sel.names)
        if (std::find(m_seen.begin(), m_seen.end(), n) == m_seen.end())
          missing += (missing.empty() ? "" : ", ") + n;
      if (!missing.empty())
      {
        msg = std::string(missing_hint) + ": " + missing;
        return aio::make_error(aio::Errc::Sink, 0, 0, "selected properties not found");
      }
      return {};
    }

    const std::string &error() const { return m_error; }
    // Properties selected by "*" but not stored (name too long, grid field).
    const std::vector<std::string> &skipped() const { return m_skipped; }

  private:
    PropertySelection m_sel;
    std::vector<std::string> m_reserved; // reserved_field_names()
    std::vector<std::string> m_seen, m_rejected, m_skipped;
    std::string m_error;
  };

  // aio stream sink, shared by the text formats: type map and placement of each
  // atom; insertion is left to GridFiller (rank 0).
  template <class GridT> class GridStreamSink
  {
  public:
    GridStreamSink(Domain &domain, GridT &grid, const ReaderSettings &settings) : m_domain(domain), m_grid(grid), m_species(settings.species), m_bounds_mode(settings.bounds_mode), m_read_velocities(settings.read_velocities), m_file_has_pbc(settings.file_has_pbc), m_props(settings.properties, reserved_field_names(typename GridT::Fields{})) {}

    // Header known: box, domain and grid geometry, before any atom line.
    aio::Status begin(const aio::atoms::FrameInfo &info)
    {
      m_ntypes = info.ntypes;
      if (aio::Status st = m_place.begin(info, m_domain, m_bounds_mode, m_file_has_pbc, "extxyz: no Lattice"); !st.ok())
        return st;
      set_grid_geometry(m_grid, m_domain);
      m_filler.reset(info.natoms);
      return {};
    }

    aio::atoms::FieldMask wanted() const
    {
      using F = aio::atoms::Field;
      // Without velocities the reader skips the whole Velocities section.
      if (!m_read_velocities)
        return {F::Id, F::Type, F::X, F::Y, F::Z};
      return {F::Id, F::Type, F::X, F::Y, F::Z, F::VX, F::VY, F::VZ};
    }

    void set_type_label(std::int64_t type, std::string_view label) { m_types.set_label(type, label); }

    aio::atoms::FieldBinding bind_property(std::string_view name, std::size_t ncomp, std::size_t comp, aio::atoms::ScalarType) { return m_props.bind(m_filler, name, ncomp, comp); }

    // Before the first atom line: the type labels are complete, and we learn
    // whether velocities come on the atom lines (extxyz) rather than in a
    // section of their own (LAMMPS Velocities -> velocity()).
    void atom_fields(aio::atoms::FieldMask m)
    {
      if (!m_types.freeze(m_ntypes, m_species, m_msg))
        m_type_error = true;
      m_vel_in_atoms = m.has(aio::atoms::Field::VX) || m.has(aio::atoms::Field::VY) || m.has(aio::atoms::Field::VZ);
    }

    aio::Errc atom(std::size_t index, const aio::atoms::AtomValues &a)
    {
      if (m_type_error)
        return aio::Errc::Sink;
      if (!m_props.error().empty())
      {
        m_msg = m_props.error();
        return aio::Errc::Sink;
      }
      if (a.type < 1 || static_cast<std::size_t>(a.type) > m_types.size())
      {
        m_msg = "atom type " + std::to_string(a.type) + " outside [1, " + std::to_string(m_types.size()) + "]";
        return aio::Errc::Sink;
      }
      Vec3d r;
      std::uint32_t cell;
      if (!m_place.place(Vec3d{a.x, a.y, a.z}, r, cell))
      {
        m_msg = "atom outside the domain on a non-periodic axis";
        return aio::Errc::Sink;
      }
      m_filler.stage(index, cell, r, m_types[static_cast<std::size_t>(a.type - 1)]);
      if (m_vel_in_atoms)
        m_filler.stage_velocity(index, a.vx, a.vy, a.vz);
      ++m_count;
      return aio::Errc::Ok;
    }

    aio::Errc velocity(std::size_t index, double vx, double vy, double vz)
    {
      m_filler.stage_velocity(index, vx, vy, vz);
      return aio::Errc::Ok;
    }

    aio::Status end()
    {
      if (m_types.late_label())
        return aio::make_error(aio::Errc::Unsupported, 0, 0, "type labels after the Atoms section are not supported");
      return m_props.check(m_msg, "not in the file, string-valued, or read as a grid field");
    }

    // Puts the staged atoms into the grid (after a successful read). id = index.
    void commit() { m_filler.commit(m_grid); }

    // After commit() and rebuild_particle_offsets(): generic properties ->
    // grid flat arrays. Returns the fields created.
    std::vector<GenericFieldInfo> fill_properties() { return m_filler.fill_properties(m_grid); }
    const std::vector<std::string> &skipped_properties() const { return m_props.skipped(); }
    const std::string &message() const { return m_msg; }
    std::size_t count() const { return m_count; }
    const BoxTransform &box() const { return m_place.box(); }

  private:
    Domain &m_domain;
    GridT &m_grid;
    const ParticleSpecies *m_species;
    ReadBoundsSelectionMode m_bounds_mode;
    bool m_read_velocities = true;
    bool m_file_has_pbc = false;
    PropertyFilter m_props;
    GridPlacer m_place;
    TypeMap m_types;
    GridFiller<GridT> m_filler;
    std::size_t m_ntypes = 0, m_count = 0;
    bool m_type_error = false;
    bool m_vel_in_atoms = false;
    std::string m_msg;
  };

  // Result sent to every rank (from rank 0, or from the first failing rank).
  struct ReadOutcome
  {
    int ok = 1;
    char msg[sizeof(aio::Error::msg) + 512] = {};
  };

  // Every rank learns the outcome; fatal_error everywhere on failure, so no
  // rank waits forever. Parallel: the first failing rank's message.
  inline void agree_on_outcome(const char *op_name, MPI_Comm comm, bool parallel, ReadOutcome &out)
  {
    int rank = 0, nranks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);
    int from = 0;
    if (parallel)
    {
      from = out.ok ? nranks : rank;
      MPI_Allreduce(MPI_IN_PLACE, &from, 1, MPI_INT, MPI_MIN, comm);
    }
    if (from < nranks)
      MPI_Bcast(&out, sizeof(ReadOutcome), MPI_BYTE, from, comm);
    if (!out.ok)
      fatal_error() << op_name << ": " << out.msg << std::endl;
  }

  // Shared flow of every atom file reader. `make_sink()` builds the sink;
  // `read_fn(path, sink)` reads the file into it and returns an aio::Status.
  //   - serial (parallel = false): rank 0 reads, as the previous
  //     read_external_atom_format_v2; the other ranks get the geometry.
  //   - parallel: every rank reads its share (HDF5 MPIO) into its own grid; the
  //     deck's load balancing / migration redistributes the atoms, as after any read.
  // A sink may have `std::string global_check(MPI_Comm, bool parallel)` (empty: ok),
  // run after a successful read on the reading ranks (collective when parallel).
  template <class GridT, class MakeSink, class ReadFn> void read_into_grid(const char *op_name, MPI_Comm comm, Domain &domain, GridT &grid, const std::string &path, bool parallel, MakeSink &&make_sink, ReadFn &&read_fn)
  {
    lout << "======== " << op_name << " ========" << std::endl;
    int rank = 0, nranks = 1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nranks);

    if (grid.number_of_particles() != 0)
    {
      lerr << op_name << ": system already defined, skipping" << std::endl;
      return;
    }
    // The deck's domain (cell size, dims, bounds, periodic, expandable, xform) is
    // kept: BoxTransform fits it to the file cell (ADR 0010).

    const bool reads = parallel || rank == 0;
    ReadOutcome out;
    std::vector<GenericFieldInfo> generic; // created by the readers, then on every rank (serial)
    unsigned long long count = 0;
    bool uniform_scale = true;
    bool checked = false;
    double ms_read = 0, ms_check = 0, ms_grid = 0;
    {
      auto sink = make_sink();
      using clock = std::chrono::steady_clock;
      if (parallel)
        MPI_Barrier(comm); // times the read only (MPI's first collectives may set up connections)
      const auto t0 = clock::now();
      aio::Status st;
      if (reads)
        st = read_fn(path.c_str(), sink);
      const auto t1 = clock::now();
      if (reads && !st.ok())
      {
        out.ok = 0;
        const aio::Error &e = st.error();
        const char *sep = sink.message().empty() ? "" : " -- ";
        if (e.line > 0)
          std::snprintf(out.msg, sizeof out.msg, "%s:%llu: %s%s%s", path.c_str(), static_cast<unsigned long long>(e.line), e.msg, sep, sink.message().c_str());
        else
          std::snprintf(out.msg, sizeof out.msg, "%s: %s%s%s", path.c_str(), e.msg, sep, sink.message().c_str());
      }
      agree_on_outcome(op_name, comm, parallel, out);

      const auto t2 = clock::now();
      if constexpr (requires { sink.global_check(comm, parallel); })
      {
        checked = true;
        if (reads)
        {
          const std::string err = sink.global_check(comm, parallel);
          if (!err.empty())
          {
            out.ok = 0;
            std::snprintf(out.msg, sizeof out.msg, "%s: %s", path.c_str(), err.c_str());
          }
        }
        agree_on_outcome(op_name, comm, parallel, out);
      }

      const auto t3 = clock::now();
      if (reads)
      {
        sink.commit();
        grid.rebuild_particle_offsets(); // flat arrays are indexed by these offsets
        generic = sink.fill_properties();
        count = sink.count();
        uniform_scale = sink.box().uniform_scale;
        if (rank == 0)
          for (const std::string &sk : sink.skipped_properties())
            lout << "Skipped property = " << sk << std::endl;
      }
      ms_read = std::chrono::duration<double, std::milli>(t1 - t0).count();
      ms_check = std::chrono::duration<double, std::milli>(t3 - t2).count();
      ms_grid = std::chrono::duration<double, std::milli>(clock::now() - t3).count();
    }

    // Domain: rank 0's for everyone. Parallel readers computed it from the same file cell: it must match.
    const Domain mine = domain;
    MPI_Bcast(&domain, sizeof(Domain), MPI_BYTE, 0, comm);
    if (parallel)
    {
      const bool same = mine.bounds().bmin == domain.bounds().bmin && mine.bounds().bmax == domain.bounds().bmax && mine.cell_size() == domain.cell_size() && mine.grid_dimension() == domain.grid_dimension() && mine.xform() == domain.xform();
      int differs = same ? 0 : 1;
      MPI_Allreduce(MPI_IN_PLACE, &differs, 1, MPI_INT, MPI_MAX, comm);
      if (differs)
        fatal_error() << op_name << ": ranks computed different domains from the file cell" << std::endl;
    }
    else if (rank != 0)
    {
      set_grid_geometry(grid, domain);
    }
    grid.rebuild_particle_offsets();

    // Same generic fields on every rank (empty where there are no particles; parallel: each reader made them).
    int ngeneric = static_cast<int>(generic.size());
    MPI_Bcast(&ngeneric, 1, MPI_INT, 0, comm);
    if (ngeneric > 0 && !parallel)
    {
      generic.resize(static_cast<std::size_t>(ngeneric));
      MPI_Bcast(generic.data(), static_cast<int>(ngeneric * sizeof(GenericFieldInfo)), MPI_BYTE, 0, comm);
      if (rank != 0)
        for (const GenericFieldInfo &f : generic)
          generic_field_data(grid, f);
    }

    if (parallel)
    {
      MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
      double ms[3] = {ms_read, ms_check, ms_grid};
      MPI_Allreduce(MPI_IN_PLACE, ms, 3, MPI_DOUBLE, MPI_MAX, comm);
      ms_read = ms[0], ms_check = ms[1], ms_grid = ms[2];
    }
    lout << onika::format_string("Read %llu atoms in %.1f ms (read %.1f ms, ", count, ms_read + ms_check + ms_grid, ms_read) << (checked ? onika::format_string("id check %.1f ms, ", ms_check) : std::string()) << onika::format_string("grid insert %.1f ms; %s", ms_grid, parallel ? "MPIO on " : "rank 0") << (parallel ? std::to_string(nranks) + (nranks > 1 ? " ranks" : " rank") : std::string()) << ")" << std::endl;
    lout << "Uniform scale    = " << std::boolalpha << uniform_scale << std::endl;
    lout << "Particles        = " << count << std::endl;
    lout << "Domain XForm     = " << domain.xform() << std::endl;
    lout << "Domain bounds    = " << domain.bounds() << std::endl;
    lout << "Domain size      = " << bounds_size(domain.bounds()) << std::endl;
    lout << "Real size        = " << bounds_size(domain.bounds()) * Vec3d{domain.xform().m11, domain.xform().m22, domain.xform().m33} << std::endl;
    lout << "Cell size        = " << domain.cell_size() << std::endl;
    lout << "Grid dimensions  = " << domain.grid_dimension() << " (" << grid_cell_count(domain.grid_dimension()) << " cells)" << std::endl;
    lout << "Domain periodic  = " << std::boolalpha << domain.periodic_boundary_x() << "," << domain.periodic_boundary_y() << "," << domain.periodic_boundary_z() << std::endl;
    lout << "Domain expandable = " << std::boolalpha << domain.expandable() << std::endl;
    for (const GenericFieldInfo &f : generic)
      lout << "Generic field    = " << f.name << " (" << generic_kind_name(f.kind) << ")" << std::endl;

    if (!check_domain(domain))
    {
      lout << "domain = " << domain << std::endl;
      fatal_error() << op_name << ": invalid domain configuration" << std::endl;
    }
    lout << "===========================================" << std::endl << std::endl;
  }

  // The text readers: rank 0 reads through GridStreamSink.
  template <class GridT, class ReadFn> void read_atom_file(const char *op_name, MPI_Comm comm, Domain &domain, GridT &grid, const ReaderSettings &settings, const std::string &path, ReadFn &&read_fn)
  {
    read_into_grid(op_name, comm, domain, grid, path, false, [&] { return GridStreamSink<GridT>(domain, grid, settings); }, read_fn);
  }

} // namespace exaStamp::atom_io
