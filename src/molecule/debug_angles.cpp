#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>

#include <memory>
#include <mpi.h>

namespace exanb
{

  // =====================================================================
  // ========================== TestAngles        ========================
  // =====================================================================

  using ChemicalAngles    = std::vector< std::array<uint64_t,3> >;

  struct DebugAnglesNode : public OperatorNode
  {
    ADD_SLOT( MPI_Comm        , mpi             , INPUT );
    ADD_SLOT( ChemicalAngles  , chemical_angles    , INPUT );

    void execute() override final
    {

      ChemicalAngles angles = *chemical_angles;

      // MPI Initialization
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      lout << "nb of proc : " << np << std::endl;
      for(int p=0; p<np; ++p)
        {
          if(rank==p)
            for(const std::array<uint64_t,3>& angle : angles)
              std::cout << angle[0] << " " << angle[1] << " " << angle[2] << " " << std::endl;

          MPI_Barrier(*mpi);
        }
    }
  };

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "debug_angles", make_compatible_operator<DebugAnglesNode> );

  }

}
