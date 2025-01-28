#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/domain.h>

#include <mpi.h>

namespace exaStamp
{
  using namespace exanb;
  using namespace onika;

  class DomainVolume : public OperatorNode
  {  
    ADD_SLOT( MPI_Comm, mpi, INPUT, MPI_COMM_WORLD, DocString{"MPI communicator"} );
    ADD_SLOT( Domain  , domain, INPUT, REQUIRED  , DocString{"Simulation domain"});
    ADD_SLOT( double  , volume, INPUT_OUTPUT      , DocString{"Computed volume"}  );

  public:
    inline bool is_sink() const override final { return true; }

    inline void execute () override final
    {
       *volume = domain->bounds().bmax.x - domain->bounds().bmin.x;
       lout << "Volume = " << *volume << std::endl;
    }

    inline std::string documentation() const override final
    {
      return "This is the documentation of domain_volume component. (c) 2023 T. Carrard";
    }
  };
  
 // === register factories ===  
  ONIKA_AUTORUN_INIT(domain_volume)
  {
   OperatorNodeFactory::instance()->register_factory( "domain_volume", make_simple_operator< DomainVolume > );
  }

}

