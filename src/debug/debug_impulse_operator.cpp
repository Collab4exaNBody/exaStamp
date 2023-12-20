#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exaStamp/particle_species/particle_specie.h>

#include <exanb/particle_neighbors/grid_particle_neighbors.h>

#include <memory>

namespace exaStamp
{
  using namespace exanb;

  // =====================================================================
  // ========================== TestWeight ========================
  // =====================================================================

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_vx , field::_vy , field::_vz , field::_type>
    >
  class TestImpulseNode : public OperatorNode
  {
    ADD_SLOT( GridT , grid , INPUT);
    ADD_SLOT(ParticleSpecies        , species       , INPUT );

  public:
    void execute() override final
    {
      const GridT& grid = *(this->grid);
      const ParticleSpecies species = *(this->species);

      auto cells = grid.cells();
      size_t n_cells = grid.number_of_cells();

      Vec3d total_impulse = {0,0,0};

      for(size_t cell_a=0;cell_a<n_cells;cell_a++)
        {
          if(! grid.is_ghost_cell(cell_a) )
            {
              size_t n_particles = cells[cell_a].size();

              for(size_t p_a=0;p_a<n_particles;p_a++)
                {
                  total_impulse += {species.at(cells[cell_a][field::type][p_a]).m_mass * cells[cell_a][field::vx][p_a],
                                   species.at(cells[cell_a][field::type][p_a]).m_mass * cells[cell_a][field::vy][p_a],
                                   species.at(cells[cell_a][field::type][p_a]).m_mass * cells[cell_a][field::vz][p_a]};
                }
            }
        }
      lout << "[DEBUG] total momentum of the system : [" << total_impulse << "]. Must be close of 0 or at least constant." << std::endl;
    }
  };

  template<class GridT> using TestImpulseNodeTmpl = TestImpulseNode<GridT>;
  
  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "debug_impulse", make_grid_variant_operator< TestImpulseNodeTmpl > );
  }

}
