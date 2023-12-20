#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/field_sets.h>
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
  CONSTRUCTOR_FUNCTION
  {
	  OperatorNodeFactory::instance()->register_factory(
        "grid_flavor_minimal",
        make_compatible_operator< InitGridFlavorNode< GridFromFieldSet<MinimalFieldSet> > >
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

