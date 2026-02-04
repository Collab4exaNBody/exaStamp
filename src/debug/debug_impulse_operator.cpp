/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>
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
  ONIKA_AUTORUN_INIT(debug_impulse_operator)
  {
    OperatorNodeFactory::instance()->register_factory( "debug_impulse", make_grid_variant_operator< TestImpulseNodeTmpl > );
  }

}
