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

#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid_fields.h>
#include <exanb/core/grid.h>

#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;

  // 
  template<class GridT , class FieldSubSetT = typename GridT::field_set_t>
  struct InitGridFlavorNode : public OperatorNode
  {
    ADD_SLOT(GridT, grid, INPUT_OUTPUT );

    inline InitGridFlavorNode()
    {
      set_profiling(false);
    }

    inline void execute () override final
    {
      if( grid->number_of_cells() == 0 )
      {
        grid->set_cell_allocator_for_fields( FieldSubSetT{} );
        grid->rebuild_particle_offsets();
      }
    }      
  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(grid_flavor)
  {
    // just an alias to multimat
	  OperatorNodeFactory::instance()->register_factory(
        "grid_flavor_minimal",
        make_compatible_operator< InitGridFlavorNode< GridFromFieldSet<MultiMatFieldSet> > >
        );

      OperatorNodeFactory::instance()->register_factory(
        "grid_flavor_multimat",
        make_compatible_operator< InitGridFlavorNode< GridFromFieldSet<MultiMatFieldSet> > >
        );

      OperatorNodeFactory::instance()->register_factory(
        "grid_flavor_full",
        make_compatible_operator< InitGridFlavorNode< GridFromFieldSet<MoleculeFieldSet> > >
        );

      OperatorNodeFactory::instance()->register_factory(
        "grid_flavor_full_mechanics",
        make_compatible_operator< InitGridFlavorNode< GridFromFieldSet<FullFieldMechSet> > >
        );
      
	OperatorNodeFactory::instance()->register_factory(
        "grid_flavor_multimat_mechanics",
        make_compatible_operator< InitGridFlavorNode< GridFromFieldSet<MultimatMechFieldSet> > >
        );
      
      OperatorNodeFactory::instance()->register_factory(
        "grid_flavor_rigidmol",
        make_compatible_operator< InitGridFlavorNode< GridFromFieldSet<RigidMoleculeFieldSet> > >
        );
  }

}

