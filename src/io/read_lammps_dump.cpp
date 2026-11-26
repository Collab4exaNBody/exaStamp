// read_lammps_dump: reads the first frame of a LAMMPS text dump with aio (ADR 0012).
//
//   read_lammps_dump: { file: dump.lammpstrj }                  # .gz/.bz2/.xz detected from content
//   read_lammps_dump: { file: dump.lammpstrj, fields: ["*"] }   # every other per-atom column
//   read_lammps_dump: { file: dump.lammpstrj, fields: [c_pe] }  # these columns only
//
// BOX BOUNDS gives the cell (orthogonal, restricted or general triclinic).
// Positions come from the first complete column family of x > xs > xu > xsu
// (scaled ones are converted, unwrapped ones wrapped into the box). Types: the
// type column (type t is species t-1), else the element column (matched to
// the species names); the species slot is required. Velocities: vx vy vz
// (read_velocities). Atom ids: the id column if there is one (remap_atom:
// false, ids must be 1..N), else the order in the file.
//
// `fields` selects other per-atom columns (one scalar each) to store as exanb
// generic fields; see read_extxyz. Generic fields are valid until particles
// first move.
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

#include <exaStamp/io/atom_file_reader.h>

namespace exaStamp {

using namespace exanb;

template <class GridT, class = AssertGridHasFields<GridT, field::_rx, field::_ry, field::_rz, field::_vx, field::_vy,
                                                   field::_vz, field::_id, field::_type>>
class ReadLammpsDumpNode : public OperatorNode {
  using StringList = std::vector<std::string>;
  EXASTAMP_ATOM_READER_SLOTS();
  ADD_SLOT(StringList, fields, INPUT, StringList{},
           DocString{"Per-atom columns to store as generic scalar fields; \"*\" = all. Default: none"});

public:
  inline void execute() override final {
    if (!species.has_value() || species->empty())
      fatal_error() << "read_lammps_dump: species must be defined (dumps do not state the type count)" << std::endl;

    std::vector<std::string_view> names;
    for (std::size_t k = 0; k < species->size(); ++k)
      names.emplace_back(species->at(k).m_name);

    aio::formats::lammps_dump::Options opt;
    opt.by_id_if_present = !*remap_atom;
    opt.species = names;

    atom_io::ReaderSettings settings = EXASTAMP_ATOM_READER_SETTINGS();
    settings.properties = atom_io::PropertySelection(*fields);
    settings.file_has_pbc = true;
    atom_io::read_atom_file("read_lammps_dump", *mpi, *domain, *grid, settings, *file,
                            [&](const char* path, auto& sink) {
                              return aio::formats::lammps_dump::read_file(path, sink, opt);
                            });
  }

  inline std::string documentation() const override final {
    return R"EOF(
Reads the first frame of a LAMMPS text dump into the grid: positions (first
complete family of x > xs > xu > xsu), types (type column: type t is species
t-1; else element column, matched to the species names) and velocities
(read_velocities). Other per-atom columns listed in `fields` ("*" = all) are
stored as scalar generic fields, valid until particles first move.

Example:

  - read_lammps_dump:
      file: dump.lammpstrj
      fields: [ c_pe ]
)EOF";
  }
};

template <class GridT> using ReadLammpsDumpNodeTmpl = ReadLammpsDumpNode<GridT>;

ONIKA_AUTORUN_INIT(read_lammps_dump) {
  OperatorNodeFactory::instance()->register_factory("read_lammps_dump",
                                                    make_grid_variant_operator<ReadLammpsDumpNodeTmpl>);
}

} // namespace exaStamp
