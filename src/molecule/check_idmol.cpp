#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/molecule/molecule_list.h>
#include <exaStamp/molecule/molecule_species.h>

#include <unordered_map>

#include <omp.h>

namespace exaStamp
{
  using namespace exanb;

  template< class GridT >
  class CheckMoleculeIdMol : public OperatorNode
  {
    ADD_SLOT( GridT    , grid   , INPUT_OUTPUT);
    ADD_SLOT( MoleculeSpeciesVector , molecules  , INPUT, REQUIRED, DocString{"Molecule descriptions"} );

  public:
    inline void execute ()  override final
    {
      if( ! grid->has_allocated_field(field::idmol) )
      {
        fatal_error() << "check_idmol can only work if particle field 'idmol' is present" << std::endl;
      }
    
      auto cells = grid->cells_accessor();
      auto field_idmol = grid->field_const_accessor( field::idmol );
      const size_t n_cells = grid->number_of_cells();
      
#     pragma omp parallel
      {              
#       pragma omp for schedule(guided)        
        for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          bool is_ghost = grid->is_ghost_cell( cell_i );
          const size_t n_particles = cells[cell_i].size();
          for(size_t i=0;i<n_particles;i++)
          {
            const auto idmol = cells[cell_i][field_idmol][i];
            if( idmol == std::numeric_limits<uint64_t>::max() )
            {
              fatal_error() << "Unitialized idmol @ Cell#"<<cell_i<<" / P#"<<i<<" , ghost="<<std::boolalpha<<is_ghost <<std::endl;
            }
            int mtype = molecule_type_from_id( idmol );
            if( mtype<0 || mtype >= molecules->m_molecules.size() )
            {
              fatal_error() << "Bad molecule type "<<mtype<<" , number of molecule species = "<<molecules->m_molecules.size()<<" , ghost="<<std::boolalpha<<is_ghost<<std::endl;
            }
            int place = molecule_place_from_id( idmol );
            if( place<0 || place >= molecules->m_molecules[mtype].m_nb_atoms )
            {
              fatal_error() << "Bad atom place "<<place<<", molecule atoms = "<<molecules->m_molecules[mtype].m_nb_atoms<<" , ghost="<<std::boolalpha<<is_ghost<<std::endl;
            }
          }
        }
      }
      
    }

  };

  template<class GridT> using CheckMoleculeIdMolTmpl = CheckMoleculeIdMol<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("check_idmol", make_grid_variant_operator< CheckMoleculeIdMolTmpl > );
  }

}
