// read_extxyz: reads the first frame of an extended XYZ file with aio (ADR 0011).
//
//   read_extxyz: { file: system.xyz }                         # .gz/.bz2/.xz detected from content
//   read_extxyz: { file: system.xyz, fields: ["*"] }          # every other per-atom property
//   read_extxyz: { file: system.xyz, fields: [force, eps] }   # these properties only
//
// Lattice gives the cell (required); pos, type (or species, or Z) and velo go
// to the grid (velocities: read_velocities). Atom ids: the id column if there
// is one (remap_atom: false), else the order in the file. Types come from the
// type column when there is one (type t is species t-1); otherwise species
// symbols (or the element of Z) are matched to the species names.
//
// `fields` selects other per-atom properties to store as exanb generic fields
// (grid flat arrays): 1 column -> scalar, 3 -> Vec3, 9 -> Mat3 (row-major),
// any other count -> scalars name_0..name_{n-1}. Integers and logicals are
// stored as doubles. Properties not selected are skipped without parsing.
// Generic fields are valid until particles first move (exanb does not
// migrate flat arrays): use them in setup_system (read -> project -> dump).
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
class ReadExtXYZNode : public OperatorNode {
  using StringList = std::vector<std::string>;
  EXASTAMP_ATOM_READER_SLOTS();
  ADD_SLOT(StringList, fields, INPUT, StringList{},
           DocString{"Per-atom properties to store as generic fields (scalar/Vec3/Mat3); \"*\" = all. Default: none"});

public:
  inline void execute() override final {
    if (!species.has_value() || species->empty())
      fatal_error() << "read_extxyz: species must be defined (atom types are matched to species names)" << std::endl;

    std::vector<std::string_view> names;
    for (std::size_t k = 0; k < species->size(); ++k)
      names.emplace_back(species->at(k).m_name);

    aio::formats::extxyz::Options opt;
    opt.by_id_if_present = !*remap_atom;
    opt.species = names;

    atom_io::ReaderSettings settings = EXASTAMP_ATOM_READER_SETTINGS();
    settings.properties = atom_io::PropertySelection(*fields);
    settings.file_has_pbc = true;
    atom_io::read_atom_file("read_extxyz", *mpi, *domain, *grid, settings, *file,
                            [&](const char* path, auto& sink) {
                              return aio::formats::extxyz::read_file(path, sink, opt);
                            });
  }

  inline std::string documentation() const override final {
    return R"EOF(
Reads the first frame of an extended XYZ file (Lattice required) into the grid:
positions, types (type column: type t is species t-1; else species / Z column,
matched to the species names) and
velocities (read_velocities). Other per-atom properties listed in `fields`
("*" = all) are stored as generic fields: 1 column -> scalar, 3 -> Vec3,
9 -> Mat3, other counts -> name_0..name_{n-1}. Generic fields are valid until
particles first move.

Example:

  - read_extxyz:
      file: system.xyz
      fields: [ force, stress ]
)EOF";
  }
};

template <class GridT> using ReadExtXYZNodeTmpl = ReadExtXYZNode<GridT>;

ONIKA_AUTORUN_INIT(read_extxyz) {
  OperatorNodeFactory::instance()->register_factory("read_extxyz", make_grid_variant_operator<ReadExtXYZNodeTmpl>);
}

} // namespace exaStamp
