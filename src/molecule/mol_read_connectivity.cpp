#include <yaml-cpp/yaml.h>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/file_utils.h>

#include <exanb/core/log.h>

#include <memory>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_id, field::_cmol>
    >
  class MolReadConnectivityNode : public OperatorNode
  {
    ADD_SLOT( GridT       , grid              , INPUT_OUTPUT );
    ADD_SLOT( std::string , filename , INPUT );

  public:
    inline void execute ()  override final
    {
      GridT& grid = *(this->grid);
      std::string file_name = data_file_path( *filename );

      //read connectivity file
      YAML::Node node = YAML::LoadFile(file_name);
      YAML::Node cNode = node["connectivity"];

      if(!cNode.IsMap())
        {
          lerr << "Impossible to read the connectivity file "<< file_name << " : node connectivity need to be a map." << std::endl;
          abort();
        }


      size_t n_cells = grid.number_of_cells();
      auto cells = grid.cells();

#pragma omp parallel
      {
#pragma omp for
        for(size_t cell_i=0;cell_i<n_cells;cell_i++)
          {
            for(size_t i =0; i<cells[cell_i].size();++i)
              {
                const uint64_t p_id = cells[cell_i][field::id][i];

                //here, we initialiaze idmol for after
                //maybe it's better to do an operator for that
                cells[cell_i][field::idmol][i] = std::numeric_limits<uint64_t>::max();

                for(size_t pos=0; pos<cNode[(std::to_string(p_id))].size();++pos)
                  {
                    cells[cell_i][field::cmol][i].at(pos) = cNode[(std::to_string(p_id))][pos].as<uint64_t>();
                  }

                //cmol is an array of size 4, so if the connectivity is lower, we need to fill the end with max uint64_t = 18446744073709551615
                std::fill(cells[cell_i][field::cmol][i].begin()+cNode[(std::to_string(p_id))].size(), cells[cell_i][field::cmol][i].end(), std::numeric_limits<uint64_t>::max());
              }
          }
      }

      // for(size_t cell_i=0;cell_i<n_cells;cell_i++)
      //   {
      //     size_t n = cells[cell_i].size();
      //     for(size_t i=0;i<n;i++)
      //       {
      //         std::cout << "Atom number " << cells[cell_i][field::id][i] << " associate to atomes " << cells[cell_i][field::cmol][i].at(0) << " " << cells[cell_i][field::cmol][i].at(1) << " " << cells[cell_i][field::cmol][i].at(2) << " " << cells[cell_i][field::cmol][i].at(3) << " "  << std::endl;
      //       }
      //   }


    }
  };

    template<class GridT> using MolReadConnectivityNodeTmpl = MolReadConnectivityNode<GridT>;

    // === register factories ===
    CONSTRUCTOR_FUNCTION
    {
      OperatorNodeFactory::instance()->register_factory( "mol_read_connectivity", make_grid_variant_operator< MolReadConnectivityNodeTmpl >);
    }

  }
