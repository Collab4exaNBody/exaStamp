#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/domain.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include "kim_init.h"

namespace exaStamp
{

  using namespace exanb;

  class KIMInitOperator : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( std::string, model     , INPUT , REQUIRED );

  public:
    // Operator execution
    inline void execute () override final
    {

      std::cout << "KIM Initialization function" << std::endl;
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(kim_init)
  {  
    OperatorNodeFactory::instance()->register_factory( "kim_init" , make_simple_operator< KIMInitOperator > );
  }

}
