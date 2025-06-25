#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>

#include <onika/log.h>
#include <exaStamp/molecule/molecule_species.h>

namespace exaStamp
{
  using namespace exanb;
  
  class MoleculeSpeciesSetup : public OperatorNode
  {    
    ADD_SLOT( MoleculeSpeciesVector , molecules , INPUT_OUTPUT , DocString{"Molecule descriptions"} );

  public:
    inline bool is_sink() const override final { return true; }
    
    inline void execute ()  override final
    {
      for(const auto& mol:molecules->m_molecules)
      {
        ldbg << "+ molecule "<< mol.name() << std::endl;
      }
    }

  };

  // === register factories ===
  ONIKA_AUTORUN_INIT(molecule_species)
  {
    OperatorNodeFactory::instance()->register_factory( "molecule_species", make_simple_operator< MoleculeSpeciesSetup > );
  }

}

