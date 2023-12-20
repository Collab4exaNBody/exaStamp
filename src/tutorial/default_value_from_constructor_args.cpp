#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/log.h>
#include <exanb/core/cpp_utils.h>

#include <vector>

namespace exaStamp
{
  using namespace exanb;

  static const std::vector<double> myvec = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };

  class DefaultValueFromCTorArgs : public OperatorNode
  {      
    // ========= I/O slots =======================
    ADD_SLOT( std::vector<double> , input1, INPUT, std::make_tuple(size_t(10),5.0) );
    ADD_SLOT( std::vector<double> , input2, INPUT, std::make_tuple(myvec.begin()+1,myvec.begin()+6) );

  public:
    // Operator execution
    inline void execute () override final
    {      
      lout<<"input1 : length = " << input1->size()<<", content =";
      for(const auto& x:*input1) { lout << " " << x; }
      lout << std::endl;

      lout<<"input2 : length = " << input2->size()<<", content =";
      for(const auto& x:*input2) { lout << " " << x; }
      lout << std::endl;
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {  
    OperatorNodeFactory::instance()->register_factory( "default_slot_value_from_ctor_args" , make_simple_operator< DefaultValueFromCTorArgs > );
  }

}


