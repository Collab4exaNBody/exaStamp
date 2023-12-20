#include <yaml-cpp/yaml.h>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/basic_types_yaml.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/molecule/id_map.h>
#include <exaStamp/molecule/mol_connectivity.h>

#include <exanb/core/log.h>

#include <memory>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_id, field::_cmol>
    >
  class ReorganizeMoleculeImpropers : public OperatorNode
  {
    //using ChemicalImpropers = std::vector< std::array<uint64_t,4> >;

    ADD_SLOT( ParticleSpecies          , species               , INPUT );
    ADD_SLOT( IdMap                    , id_map                , INPUT );
    ADD_SLOT( IdMapGhosts              , id_map_ghosts         , INPUT );
    ADD_SLOT( ChemicalImpropers        , chemical_impropers    , INPUT_OUTPUT );

  public:
    inline void execute ()  override final
    {
      ParticleSpecies& species     = *(this->species);
      IdMap& id_map                = *(this->id_map);
      IdMapGhosts& id_map_ghosts   = *(this->id_map_ghosts);
      ChemicalImpropers& impropers = *chemical_impropers;


      size_t       cell[4];
      size_t       pos [4];
      unsigned int type[4];

      for(size_t i =0; i<impropers.size();++i)
        {
          uint64_t a0 = impropers[i].at(0);
          uint64_t a1 = impropers[i].at(1);
          uint64_t a2 = impropers[i].at(2);
          uint64_t a3 = impropers[i].at(3);

          //lerr << "In src/molecule/reorganise_molecule_improper.cpp impropers considered between atom global ids "  << a0 << " " << a1 << " " << a2 << " " << a3 << std::endl << std::flush;

          if(id_map.find(a0) != id_map.end())
            {
              decode_cell_particle(id_map.at(a0), cell[0], pos[0], type[0]);
            }
          else if(id_map_ghosts.find(a0) != id_map_ghosts.end())
            {
              decode_cell_particle(id_map_ghosts.find(a0)->second, cell[0], pos[0], type[0]);
            }
          else
            {
              lerr << "Reorganize_improper_tatb : id " << a0 << " not in id_map or id_map_ghost. Maybe id_map is missing in the config file." << std::endl;
              abort();
            }

          if(id_map.find(a1) != id_map.end())
            {
              decode_cell_particle(id_map.at(a1), cell[1], pos[1], type[1]);
            }
          else if(id_map_ghosts.find(a1) != id_map_ghosts.end())
            {
              decode_cell_particle(id_map_ghosts.find(a1)->second, cell[1], pos[1], type[1]);
            }
          else
            {
              lerr << "Reorganize_improper_tatb : id " << a1 << " not in id_map or id_map_ghost. Maybe id_map is missing in the config file.";
                abort();
            }

          if(id_map.find(a2) != id_map.end())
            {
              decode_cell_particle(id_map.at(a2), cell[2], pos[2], type[2]);
            }
          else if(id_map_ghosts.find(a2) != id_map_ghosts.end())
            {
              decode_cell_particle(id_map_ghosts.find(a2)->second, cell[2], pos[2], type[2]);
            }
          else
            {
              lerr << "Reorganize_improper_tatb : id " << a2 << " not in id_map or id_map_ghost. Maybe id_map is missing in the config file." << std::endl;
                abort();
            }


          if(id_map.find(a3) != id_map.end())
            {
              decode_cell_particle(id_map.at(a3), cell[3], pos[3], type[3]);
            }
          else if(id_map_ghosts.find(a3) != id_map_ghosts.end())
            {
              decode_cell_particle(id_map_ghosts.find(a3)->second, cell[3], pos[3], type[3]);
            }
          else
            {
              lerr << "Reorganize_improper_tatb : id " << a3 << " not in id_map or id_map_ghost. Maybe id_map is missing in the config file.";
                abort();
            }





          if( (species[type[0]].name() == "Cno2") || (species[type[0]].name() == "CNO2") )
            {
              if( (species[type[3]].name() != "No2") && (species[type[3]].name() != "NO2") )
                {
                  if( (species[type[1]].name() == "No2") || (species[type[1]].name() == "NO2") )
                    {
                      std::array<uint64_t, 4> tmp{a0, a2, a3, a1};
                      impropers[i] = tmp;
                    }
                  else if( (species[type[2]].name() == "No2") || (species[type[2]].name() == "NO2") )
                    {
                      std::array<uint64_t, 4> tmp{a0, a1, a3, a2};
                      impropers[i] = tmp;
                    }
                }
              continue;
            }

          if( (species[type[0]].name() == "Cnh2") || (species[type[0]].name() == "CNH2") )
            {
              if( (species[type[3]].name() != "Nh2") && (species[type[3]].name() != "NH2") )
                {
                  if( (species[type[1]].name() == "Nh2") || (species[type[1]].name() == "NH2") )
                    {
                      std::array<uint64_t, 4> tmp{a0, a2, a3, a1};
                      impropers[i] = tmp;
                    }
                  else if( (species[type[2]].name() == "Nh2") || (species[type[2]].name() == "NH2") )
                    {
                      std::array<uint64_t, 4> tmp{a0, a1, a3, a2};
                      impropers[i] = tmp;
                    }
                }
              continue;
            }

          if( (species[type[0]].name() == "No2") || (species[type[0]].name() == "NO2") )
            {
              if( (species[type[3]].name() != "Cno2") && (species[type[3]].name() != "CNO2") )
                {
                  if( (species[type[1]].name() == "Cno2") || (species[type[1]].name() == "CNO2") )
                    {
                      std::array<uint64_t, 4> tmp{a0, a2, a3, a1};
                      impropers[i] = tmp;
                    }
                  else if( (species[type[2]].name() == "Cno2") || (species[type[2]].name() == "CNO2") )
                    {
                      std::array<uint64_t, 4> tmp{a0, a1, a3, a2};
                      impropers[i] = tmp;
                    }
                }
              continue;
            }


          if( (species[type[0]].name() == "Nh2") || (species[type[0]].name() == "NH2") )
            {
              //std::cout << "BEFORE : " << species[type[0]].name() << " " << species[type[1]].name() << " " << species[type[2]].name() << " " << species[type[3]].name() << std::endl;
              if( (species[type[3]].name() != "Cnh2") && (species[type[3]].name() != "CNH2") )
                {
                  if( (species[type[1]].name() == "Cnh2") || (species[type[1]].name() == "CNH2") )
                    {
                      std::array<uint64_t, 4> tmp{a0, a2, a3, a1};
                      impropers[i] = tmp;
                    }
                  else if( (species[type[2]].name() == "Cnh2") || (species[type[2]].name() == "CNH2") )
                    {
                      std::array<uint64_t, 4> tmp{a0, a1, a3, a2};
                      impropers[i] = tmp;
                    }
                }
              //std::cout << "AFTER : " << species[type[0]].name() << " " << species[type[1]].name() << " " << species[type[2]].name() << " " << species[type[3]].name() << std::endl;
              continue;
            }
        }

          // for(size_t i =0; i<impropers.size();++i)
      //   {
      //     uint64_t a0 = impropers[i].at(0);
      //     uint64_t a1 = impropers[i].at(1);
      //     uint64_t a2 = impropers[i].at(2);
      //     uint64_t a3 = impropers[i].at(3);

      //     decode_cell_particle(id_map.at(a0), cell[0], pos[0], type[0]);
      //     decode_cell_particle(id_map.at(a1), cell[1], pos[1], type[1]);
      //     decode_cell_particle(id_map.at(a2), cell[2], pos[2], type[2]);
      //     decode_cell_particle(id_map.at(a3), cell[3], pos[3], type[3]);

      //     std::cout << species[type[0]].name() << " " << species[type[1]].name() << " " << species[type[2]].name() << " " << species[type[3]].name() << std::endl;
      //   }
    }
  };

  template<class GridT> using ReorganizeMoleculeImpropersTmpl = ReorganizeMoleculeImpropers<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "reorganize_molecule_improper", make_grid_variant_operator< ReorganizeMoleculeImpropersTmpl > );
  }

}

