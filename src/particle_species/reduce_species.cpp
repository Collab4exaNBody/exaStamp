#include <exaStamp/particle_species/particle_specie_yaml.h>

#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/log.h>
#include <exanb/core/yaml_utils.h>
#include <exanb/core/particle_type_id.h>

#include <iostream>
#include <string>
#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  template<
      class GridT
    , class = AssertGridHasFields< GridT, field::_type >
    >
  struct ReduceSpecies : public OperatorNode
  {
    ADD_SLOT(MPI_Comm        , mpi     , INPUT , MPI_COMM_WORLD );
    ADD_SLOT(GridT           , grid    , INPUT_OUTPUT , DocString{"Particle grid"} );
    ADD_SLOT(ParticleSpecies , species , INPUT_OUTPUT );
    ADD_SLOT(ParticleTypeMap , particle_type_map , INPUT_OUTPUT );

    inline void execute () override final
    {
      unsigned long n_species_max = species->size();
      MPI_Allreduce(MPI_IN_PLACE,&n_species_max,1,MPI_UNSIGNED_LONG,MPI_MAX,*mpi);

      std::vector<unsigned long> specy_used( n_species_max , 0 );

      auto cells = grid->cells();
      size_t n_cells = grid->number_of_cells();
#     pragma omp parallel
      {
        #pragma omp for schedule(dynamic)
        for(size_t i=0;i<n_cells;i++)
        {
          size_t n_praticles = cells[i].size();
          const auto* __restrict__ types = cells[i][field::type];
          for(size_t j=0;j<n_praticles;j++)
          {
            assert( types[j] < specy_used.size() );
            specy_used[ types[j] ] = 1;
          }
        }
      }

      for(size_t i=0;i<species->size();i++)
      {
        if(specy_used[i]>0 && species->at(i).m_rigid_atom_count>1 )
        {
          for(size_t j=0;j<species->at(i).m_rigid_atom_count;j++)
          {
            int rmol_atom_type = species->at(i).m_rigid_atoms[j].m_atom_type;
            if( rmol_atom_type != -1 )
            {
              specy_used[rmol_atom_type] = 1;
            }
          }
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, specy_used.data(), n_species_max, MPI_UNSIGNED_LONG, MPI_MAX, *mpi);

      ParticleSpecies rspecies;
      for(size_t i=0;i<species->size();i++)
      {
        if(specy_used[i]!=0)
        {
          specy_used[i] = rspecies.size();
          rspecies.push_back( species->at(i) );
        }
        else
        {
          specy_used[i] = std::numeric_limits<unsigned long>::max(); // must not be used
        }
      }

      for(auto &sp:rspecies)
      {
        for(size_t j=0;j<sp.m_rigid_atom_count;j++)
        {
          int rmol_atom_type = sp.m_rigid_atoms[j].m_atom_type;
          if( rmol_atom_type != -1 )
          {
            sp.m_rigid_atoms[j].m_atom_type = specy_used[ rmol_atom_type ];
          }
        }
      }

      ldbg <<rspecies.size()<<" used among "<<species->size()<<std::endl;
      *species = rspecies;

#     pragma omp parallel
      {
        #pragma omp for schedule(dynamic)
        for(size_t i=0;i<n_cells;i++)
        {
          size_t n_praticles = cells[i].size();
          auto* __restrict__ types = cells[i][field::type];
          for(size_t j=0;j<n_praticles;j++)
          {
            types[j] = specy_used[types[j]];
          }
        }
      }
      
      // update particle type name map
      particle_type_map->clear();
      for(unsigned int a=0;a<species->size();a++)
      {
        (*particle_type_map) [ species->at(a).name() ] = a;
      }
    }

  };

  template<class GridT> using ReduceSpeciesTmpl = ReduceSpecies<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "reduce_species", make_grid_variant_operator<ReduceSpeciesTmpl> );
  }

}
