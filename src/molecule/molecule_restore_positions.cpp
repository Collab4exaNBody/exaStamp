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
#include <exanb/core/grid.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exaStamp/molecule/molecule_species.h>

namespace exaStamp
{

  template<typename GridT>
  class MoleculeRestorePositions : public OperatorNode
  {
    ADD_SLOT( GridT                   , grid                      , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( MoleculeFoldedPositions , molecule_folded_positions , INPUT_OUTPUT , REQUIRED , DocString{"positions saved by molecule_unfold"} );

  public:
    inline std::string documentation() const override final
    {
      return R"EOF(
Restores the particle positions saved by the last molecule_unfold, so that molecules can be unfolded for an output
in the middle of a simulation : molecule_unfold, write output, molecule_restore_positions. Particles must not have
moved between cells in between.
)EOF";
    }

    inline void execute () override final
    {
      auto cells = grid->cells();
      const size_t n_cells = grid->number_of_cells();
      const auto & fp = *molecule_folded_positions;
      size_t k = 0;
      for(size_t c=0;c<n_cells;c++) if( ! grid->is_ghost_cell(c) )
      {
        const size_t n = cells[c].size();
        if( k + n > fp.m_rx.size() ) fatal_error() << "molecule_restore_positions : particles changed since molecule_unfold" << std::endl;
        for(size_t p=0;p<n;p++,k++)
        {
          cells[c][field::rx][p] = fp.m_rx[k];
          cells[c][field::ry][p] = fp.m_ry[k];
          cells[c][field::rz][p] = fp.m_rz[k];
        }
      }
      if( k != fp.m_rx.size() ) fatal_error() << "molecule_restore_positions : particles changed since molecule_unfold" << std::endl;
    }
  };

  template<class GridT> using MoleculeRestorePositionsTmpl = MoleculeRestorePositions<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(molecule_restore_positions)
  {
   OperatorNodeFactory::instance()->register_factory( "molecule_restore_positions", make_grid_variant_operator<MoleculeRestorePositionsTmpl> );
  }

}
