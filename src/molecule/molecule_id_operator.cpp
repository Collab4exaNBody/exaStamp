#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <memory>
#include <numeric>
//#include <functional>
#include <mpi.h>

namespace exaStamp
{

  using IdMap = std::map<uint64_t, uint64_t>;

  // Recursive function
  // Add an id molecule to all atoms of the molecule
  template<typename CellsType>
  void set_molecules_id(CellsType& cells,
                        const IdMap& id_map,
                        const uint64_t& id,
                        size_t cell, size_t pos)
  {
    if(cells[cell][field::idmol][pos]==std::numeric_limits<uint64_t>::max())
      {
        //add new id to the molecule
        cells[cell][field::idmol][pos] = id;

        // look the first neighbours
        for(int i=0; i<4; ++i)
          {
            if(cells[cell][field::cmol][pos][i]!=std::numeric_limits<uint64_t>::max())
              {
                size_t new_cell;
                size_t new_pos;
                decode_cell_particle(id_map.at(cells[cell][field::cmol][pos][i]),new_cell, new_pos);

                set_molecules_id(cells, id_map, id, new_cell, new_pos);
              }
          }
      }
  }


  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_id, field::_cmol, field::_idmol>
    >
  class MoleculeIdNode : public OperatorNode
  {
    ADD_SLOT( GridT    , grid   , INPUT_OUTPUT);
    ADD_SLOT( IdMap    , id_map , INPUT);
    ADD_SLOT( MPI_Comm , mpi    , INPUT);

  public:
    inline void execute ()  override final
    {
      GridT& grid = *(this->grid);
      IdMap&  id_map = *(this->id_map);
      MPI_Comm comm = *mpi;

      auto cells = grid.cells();

      // Create table of max id to give an unique id on all molecules
      int nprocs;
      MPI_Comm_size(comm, &nprocs);
      uint64_t id=0, tab_id_max[nprocs];

      //set an id of molecules
      //we follow the connectivity list
      size_t n_cells = grid.number_of_cells();
      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          if(!grid.is_ghost_cell(cell_i))
            for(size_t i=0;i<cells[cell_i].size();i++)
              {
                //If the atom hasn't already an idmol
                if(cells[cell_i][field::idmol][i]==std::numeric_limits<uint64_t>::max())
                  {
                    id++;
                    set_molecules_id(cells, id_map, id, cell_i, i);
                  }
              }
        }

      MPI_Allgather( &id, 1, MPI_LONG_DOUBLE, &tab_id_max, 1,MPI_LONG_DOUBLE, comm);

      int rank;
      MPI_Comm_rank(comm, &rank);

      // gestion of MPI :
      // ids of molecules appear several time in function of the number of core
      // The following code ensure that all id are unique
      uint64_t sum = std::accumulate(tab_id_max, tab_id_max+rank, 0);
      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
            for(size_t i=0;i<cells[cell_i].size();i++)
              {
                cells[cell_i][field::idmol][i] += sum;
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

  template<class GridT> using MoleculeIdNodeTmpl = MoleculeIdNode<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("set_molecules_id", make_grid_variant_operator< MoleculeIdNodeTmpl > );
  }

}
