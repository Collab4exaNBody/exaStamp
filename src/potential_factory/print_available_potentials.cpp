#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>

#include <exaStamp/potential_factory/pair_potential_factory.h>

namespace exaStamp
{
  using namespace exanb;

  struct PrintAvailablePotentials : public OperatorNode
  {  
    ADD_SLOT( std::string , header , INPUT , std::string("")   , DocString{"Message header to print"} );

    // -----------------------------------------------
    // ----------- Operator documentation ------------
    inline std::string documentation() const override final
    {
      return R"EOF(
        Prints available pair potentials to the standard output.
        )EOF";
    }

    inline void execute () override final
    {
      lout << *header ;
      for( const auto & s : PairPotentialFactory::available_potentials() )
      {
        lout << s << std::endl;
      }
    }

  };
    
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "print_available_potentials", make_simple_operator<PrintAvailablePotentials> );
  }

}

