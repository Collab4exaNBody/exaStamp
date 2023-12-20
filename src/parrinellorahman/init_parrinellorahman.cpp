#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/basic_types_stream.h>
#include <exaStamp/parrinellorahman/parrinellorahman.h>
#include <exanb/core/domain.h>
#include <exanb/core/physics_constants.h>

#include <iostream>
#include <string>

namespace exaStamp
{
  using namespace exanb;

  struct InitParrinelloRahmanNode : public OperatorNode
  {
    static constexpr Mat3d all_one_matrix { 1.,1.,1., 1.,1.,1., 1.,1.,1. };

    ADD_SLOT( double , Text     , INPUT , REQUIRED , DocString{"target temperature"} );
    ADD_SLOT( double , masseNVT , INPUT , REQUIRED);
    ADD_SLOT( double , Pext     , INPUT , REQUIRED);
    ADD_SLOT( double , masseNPT , INPUT , REQUIRED);
    ADD_SLOT( Mat3d  , hmask    , INPUT , all_one_matrix , DocString{"scale factors for H matrix, use to constrain deformation"} );
    ADD_SLOT( Mat3d  , hblend   , INPUT , make_identity_matrix() , DocString{"matrix that enable blending of diagonal components of hp"} );

    ADD_SLOT( Domain , domain   , INPUT , REQUIRED);

    ADD_SLOT( ParrinelloRahmanContext , parrinello_rahman_ctx , OUTPUT );

    inline void execute () override final
    {
      *parrinello_rahman_ctx = ParrinelloRahmanContext { { *Text, *masseNVT, *Pext, *masseNPT, *hmask, *hblend } };

      static const double conv_pressure = 1.e4 * legacy_constant::atomicMass * 1e30;                            // internal units to Pascal

      lout << "=== Parrinello-Rahman scheme ===" << std::endl;
      lout << "  Text     : " << parrinello_rahman_ctx->m_config.m_Text               << " K" << std::endl;
      lout << "  masseNVT : " << parrinello_rahman_ctx->m_config.m_masseNVT           << " Da"<< std::endl;
      lout << "  Pext     : " << parrinello_rahman_ctx->m_config.m_Pext*conv_pressure << " Pa"<< std::endl;
      lout << "  masseB   : " << parrinello_rahman_ctx->m_config.m_masseB             << " Da"<< std::endl;
      lout << "  Hmask    : " << parrinello_rahman_ctx->m_config.m_hmask              << " "  << std::endl;
      lout << "  Hblend   : " << parrinello_rahman_ctx->m_config.m_hblend             << " "  << std::endl;
      lout << "================================" <<std::endl<<std::endl;

      parrinello_rahman_ctx->h = multiply( domain->xform() , diag_matrix( domain->bounds_size() ) );
      parrinello_rahman_ctx->G = multiply( transpose(parrinello_rahman_ctx->h), parrinello_rahman_ctx->h );
      parrinello_rahman_ctx->updateMembers();
      parrinello_rahman_ctx->apply_mask();
      parrinello_rahman_ctx->print( ldbg );
    }
  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "init_parrinellorahman", make_compatible_operator< InitParrinelloRahmanNode > );
  }

}


