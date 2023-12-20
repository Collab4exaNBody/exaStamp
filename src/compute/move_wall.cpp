#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/core/basic_types_stream.h>
#include <vector>
#include <iomanip>

namespace exaStamp
{
  using namespace exanb;

  class MoveWall : public OperatorNode
  {
   
    ADD_SLOT( Vec3d  , init_normal , INPUT , Vec3d{1.0,0.0,0.0} );
    ADD_SLOT( double , init_offset , INPUT , 0.0 );
    ADD_SLOT( double , init_cutoff , INPUT , REQUIRED );
    ADD_SLOT( double , init_epsilon , INPUT , make_quantity(1.0e-19,"J").convert() );

    ADD_SLOT( bool   , push_wall, INPUT , REQUIRED );
    ADD_SLOT( long   , timestep            , INPUT , REQUIRED );

    // outputs for walll
    ADD_SLOT( Vec3d  , normal , OUTPUT );
    ADD_SLOT( double , offset , OUTPUT );
    ADD_SLOT( double , cutoff , OUTPUT );
    ADD_SLOT( double , epsilon , OUTPUT );

  public:

    inline void execute () override final
    {
      *normal = *init_normal;
      *cutoff = *init_cutoff;
      *epsilon = *init_epsilon;
      if( *push_wall )
      {
        *offset = (*init_offset) + (*timestep) * 0.0001;
      }
      else
      {
        *offset = *init_offset;
      }
    }

  };
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "move_wall", make_simple_operator< MoveWall > );
  }

}

