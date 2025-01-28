#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <onika/log.h>

#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  class ApplicationVersionInfo : public OperatorNode
  {  
    ADD_SLOT( bool , show_details, INPUT , true , DocString{"Print additional details"}  );

  public:
    inline void execute () override final
    {
      lout << "exaStamp v"<< EXASTAMP_VERSION << " ("
#     ifndef NDEBUG
      << "debug"
#     else
      << "release"
#     endif
      <<")" << std::endl;
      if( *show_details )
      {
        int mpi_maj=0, mpi_min=0;
        MPI_Get_version( &mpi_maj, &mpi_min );
        lout << " , MPI v"<<mpi_maj<<'.'<<mpi_min << std::endl;
      }
    }
  };
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(version_info)
  {
   OperatorNodeFactory::instance()->register_factory( "version_info", make_simple_operator< ApplicationVersionInfo > );
  }

}

