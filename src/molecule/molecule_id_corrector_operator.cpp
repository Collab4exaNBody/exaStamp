#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <memory>
#include <functional>
#include <mpi.h>

namespace exaStamp
{

  using IdMap = std::map<uint64_t, uint64_t>;

  // Recursive function
  // Add an id molecule to all atoms of the molecule
  template<typename CellsType>
  void correct_molecules_id(CellsType& cells,
                            const IdMap& id_map,
                            const uint64_t& id,
                            size_t cell, size_t pos)
  {
    //add new id to the molecule
    cells[cell][field::idmol][pos] = id;

    // look the first neighbours
    for(int i=0; i<4; ++i)
      {
        if(cells[cell][field::idmol][pos]!=id &&
           cells[cell][field::cmol][pos][i]!=std::numeric_limits<uint64_t>::max())
          {
            size_t new_cell;
            size_t new_pos;
            decode_cell_particle(id_map.at(cells[cell][field::cmol][pos][i]),new_cell, new_pos);

            correct_molecules_id(cells, id_map, id, new_cell, new_pos);
          }
      }
  }


  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_cmol, field::_idmol>
    >
  class MoleculeIdCorrectorNode : public OperatorNode
  {
    ADD_SLOT( GridT , grid   , INPUT_OUTPUT);
    ADD_SLOT( IdMap , id_map , INPUT);

  public:
    inline void execute ()  override final
    {
      GridT& grid = *(this->grid);
      IdMap&  id_map = *(this->id_map);

      auto cells = grid.cells();

      // Correct id from ghost
      // When the ghost are updated, idmol can be change from one proc, to another
      // We correct this by keeping the min id along the molecule
      size_t n_cells = grid.number_of_cells();
      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          for(size_t pos=0;pos<cells[cell_i].size();pos++)
            {
              uint64_t id_to_test = cells[cell_i][field::idmol][pos];

              for(int i=0; i<4; ++i)
                {
                  if(cells[cell_i][field::cmol][pos][i]!=id_to_test)
                    {
                      correct_molecules_id(cells, id_map,
                                           cells[cell_i][field::cmol][pos][i]<id_to_test?cells[cell_i][field::cmol][pos][i]:id_to_test
                                           , cell_i, pos);
                    }
                }

            }
        }


      // for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      //   {
      //     size_t n = cells[cell_i].size();
      //     for(size_t i=0;i<n;i++)
      //       {
      //         std::cout << "Atom number " << cells[cell_i][field::id][i] << " associate to the molecule " << cells[cell_i][field::idmol][i] << std::endl;
      //       }
      //   }


    }

  };

  template<class GridT> using MoleculeIdCorrectorNodeTmpl = MoleculeIdCorrectorNode<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory("correct_molecules_id",make_grid_variant_operator< MoleculeIdCorrectorNodeTmpl > );
  }

}
