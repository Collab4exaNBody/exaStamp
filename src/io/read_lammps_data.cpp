// read_lammps_data: reads a LAMMPS data file (write_data output) with aio.
// Replaces read_external_atom_format_v2 for LAMMPS files.
//
//   read_lammps_data: { file: system.lmp }     # .gz/.bz2/.xz detected from content
//
// Supported: every atom style of aio (comment after "Atoms", or atomic/full
// guessed from the column count), optional image flags, orthogonal and
// restricted triclinic boxes, Masses type labels, Velocities
// (read_velocities: false skips them; velocities are then zero).
// Only id, type, x y z, vx vy vz are converted: other columns are skipped.
#include <string>

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
class ReadLammpsDataNode : public OperatorNode {
  EXASTAMP_ATOM_READER_SLOTS();
  ADD_SLOT(std::string, atom_style, INPUT, "",
           DocString{"LAMMPS atom style if the file has no '# style' comment after 'Atoms'"});

public:
  inline void execute() override final {
    aio::formats::lammps_data::Options opt;
    opt.keep_file_order = *remap_atom;
    opt.atom_style = *atom_style;
    atom_io::read_atom_file("read_lammps_data", *mpi, *domain, *grid, EXASTAMP_ATOM_READER_SETTINGS(), *file,
                            [&](const char* path, auto& sink) {
                              return aio::formats::lammps_data::read_file(path, sink, opt);
                            });
  }
};

template <class GridT> using ReadLammpsDataNodeTmpl = ReadLammpsDataNode<GridT>;

ONIKA_AUTORUN_INIT(read_lammps_data) {
  OperatorNodeFactory::instance()->register_factory("read_lammps_data",
                                                    make_grid_variant_operator<ReadLammpsDataNodeTmpl>);
}

} // namespace exaStamp
