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

// read_hdf5: reads atoms from an HDF5 file described by a schema of dataset
// paths (aio ADR 0014-0016), serially on rank 0 or with MPI-IO on every rank.
#include <string>
#include <string_view>
#include <vector>

#include <mpi.h>

#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>

#include <exaStamp/io/hdf5_reader.h>

namespace exaStamp {

using namespace exanb;

template <class GridT, class = AssertGridHasFields<GridT, field::_rx, field::_ry, field::_rz, field::_vx, field::_vy,
                                                   field::_vz, field::_id, field::_type>>
class ReadHdf5Node : public OperatorNode {
  EXASTAMP_ATOM_READER_SLOTS();
  ADD_SLOT(atom_io::Hdf5Schema, schema, INPUT, REQUIRED,
           DocString{"Dataset paths: positions, velocities, ids, types, cell, origin, pbc, fields"});
  ADD_SLOT(bool, mpio, INPUT, true,
           DocString{"Every rank reads its share with MPI-IO (parallel HDF5 only; otherwise rank 0 reads, with a "
                     "warning)"});
  ADD_SLOT(bool, stream_sink, INPUT, false,
           DocString{"Debug / benchmark: read on rank 0 through the text readers' stream sink (buffered)"});

public:
  inline void execute() override final {
#if !AIO_WITH_HDF5
    fatal_error() << "read_hdf5: exaStamp was built without HDF5 (configure with WITH_HDF5=ON and a libhdf5 "
                     "installation, e.g. -DHDF5_ROOT=<prefix>)"
                  << std::endl;
#else
    const aio::formats::hdf5::Schema& s = schema->s;
    if (s.positions.empty())
      fatal_error() << "read_hdf5: schema.positions is required" << std::endl;
    if (s.cell.empty())
      fatal_error() << "read_hdf5: schema.cell is required" << std::endl;
    if (!*remap_atom && s.ids.empty())
      fatal_error() << "read_hdf5: remap_atom: false needs schema.ids (or set remap_atom: true)" << std::endl;

    atom_io::ReaderSettings settings = EXASTAMP_ATOM_READER_SETTINGS();
    std::vector<std::string> names;
    for (const auto& p : s.properties)
      names.push_back(p.name);
    settings.properties = atom_io::PropertySelection(names);
    settings.file_has_pbc = !s.pbc.empty();

    if (*stream_sink) {
      // GridStreamSink takes the type count and labels from the reader: the species names
      std::vector<std::string_view> labels;
      if (species.has_value())
        for (std::size_t k = 0; k < species->size(); ++k)
          labels.emplace_back(species->at(k).m_name);
      aio::formats::hdf5::Options opt;
      opt.by_id_if_present = !*remap_atom;
      opt.species = labels;
      atom_io::read_atom_file("read_hdf5", *mpi, *domain, *grid, settings, *file,
                              [&](const char* path, auto& sink) {
                                return aio::formats::hdf5::read_file(path, s, sink, opt);
                              });
      return;
    }
    atom_io::read_hdf5_into_grid(*mpi, *domain, *grid, settings, *file, s, *remap_atom, *mpio);
#endif
  }

  inline std::string documentation() const override final {
    return R"EOF(
Reads atoms from an HDF5 file. The schema gives the dataset paths; each
per-atom dataset has N rows.

  positions   one (N,3) dataset, or a list of three (N) datasets [x, y, z] (required)
  velocities  same forms (optional; read_velocities)
  ids         (N) integers 1..N (needed unless remap_atom: true, which numbers atoms in row order)
  types       (N) integers: type t is species t-1 (optional: every atom type 1)
  cell        3x3 (rows = lattice vectors a, b, c) or (3) diagonal (required)
  origin      (3) (optional), pbc (3) integers (optional: compared with the domain)
  fields      { name: path } per-atom datasets (N) or (N,k) stored as generic fields
              (1 column -> scalar, 3 -> Vec3, 9 -> Mat3), valid until particles first move

cell, origin and pbc are a dataset path or an attribute, written "object@name".

With a parallel HDF5 and mpio: true, every rank reads its share of the rows
(MPI-IO) into its own grid; the atoms then go to their owners as after any
read. Otherwise rank 0 reads everything (with a warning when HDF5 has no
parallel I/O and several ranks run).

Example:

  - read_hdf5:
      file: system.h5
      schema:
        positions: /atoms/pos
        velocities: /atoms/vel
        ids: /atoms/id
        types: /atoms/type
        cell: /box@h
        origin: /box/origin
        fields: { force: /atoms/force }
)EOF";
  }
};

template <class GridT> using ReadHdf5NodeTmpl = ReadHdf5Node<GridT>;

ONIKA_AUTORUN_INIT(read_hdf5) {
  OperatorNodeFactory::instance()->register_factory("read_hdf5", make_grid_variant_operator<ReadHdf5NodeTmpl>);
}

} // namespace exaStamp
