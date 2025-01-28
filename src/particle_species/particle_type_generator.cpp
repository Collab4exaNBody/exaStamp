#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <onika/log.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/parallel/random.h>

#include <random>

namespace exaStamp
{
  using namespace exanb;
  
  // ================== Thermodynamic state compute operator ======================
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_type>
    >
  class ParticleTypeGeneratorNode : public OperatorNode
  {  

    ADD_SLOT( GridT           , grid    , INPUT_OUTPUT );
    ADD_SLOT( ParticleSpecies , species , INPUT );

  public:
    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      const ParticleSpecies& species = *(this->species);
      
      if( species.empty() ) { return; }
      
      size_t n_species = species.size();
      auto cells = grid.cells();
      size_t n_cells = grid.number_of_cells();
      
#     pragma omp parallel
      {
        auto & re = rand::random_engine();
        std::uniform_int_distribution<> idist(0,n_species-1);
#       pragma omp for
        for(size_t i=0;i<n_cells;i++)
        {
          uint8_t* __restrict__ part_types = cells[i][field::type];
          size_t n_part = cells[i].size();
          for(size_t j=0;j<n_part;j++)
          {
            part_types[j] = idist( re );
          }
        }
      }

    }
  };
    
  template<class GridT> using ParticleTypeGeneratorNodeTmpl = ParticleTypeGeneratorNode<GridT>;
    
  // === register factories ===  
  ONIKA_AUTORUN_INIT(particle_type_generator)
  {
   OperatorNodeFactory::instance()->register_factory("particle_type_generator", make_grid_variant_operator< ParticleTypeGeneratorNodeTmpl >
    );
  }

}

