#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/domain.h>
#include <exanb/core/log.h>
#include <exanb/core/cpp_utils.h>
#include <exaStamp/potential/ewald/ewald.h>

namespace exaStamp
{

  using namespace exanb;

  class EwaldInitOperator : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( double     , epsilon     , INPUT , REQUIRED );
    ADD_SLOT( double     , alpha       , INPUT , REQUIRED );
    ADD_SLOT( double     , radius      , INPUT , REQUIRED );
    ADD_SLOT( long       , kmax        , INPUT , REQUIRED );
    ADD_SLOT( Domain     , domain      , INPUT , OPTIONAL );
    ADD_SLOT( EwaldParms , ewald_config, INPUT_OUTPUT );
    ADD_SLOT( double     , rcut        , OUTPUT );
    ADD_SLOT( double     , rcut_max    , INPUT_OUTPUT , 0.0 );

  public:
    // Operator execution
    inline void execute () override final
    {
      using ewald_constants::fpe0;
      using ewald_constants::epsilonZero;

      if( domain.has_value() )
      {
        auto & p = *ewald_config;
        
        if( ( *alpha > 0.0 && *alpha != p.alpha ) || p.volume==0.0 || *radius != p.radius || *epsilon != p.epsilon || ( *kmax > 0 && *kmax != p.kmax ) )
        {
          auto domainSize = domain->bounds_size();
          auto xform = domain->xform();
          if( ! is_diagonal( xform ) )
          {
            lerr << "Domain XForm is not diagonal, cannot compute domain box size" << std::endl;
            std::abort();
          }
          domainSize = xform * domainSize;
          if( ! ( domain->periodic_boundary_x() && domain->periodic_boundary_y() && domain->periodic_boundary_z() ) )
          {
            lerr << "Domain must be entierly periodic, cannot initialize ewald." << std::endl;
            std::abort();
          }

          ewald_init_parameters( *alpha , *radius , *epsilon , *kmax , domainSize , p , ldbg<<"" );
          if( p.volume > 0.0 )
          {
            lout << "====== Ewald configuration ======" << std::endl;      
            lout << "size    = "<< domainSize << std::endl;
            lout << "period. = "<< std::boolalpha << domain->periodic_boundary_x() 
                                << " , " << std::boolalpha << domain->periodic_boundary_y() 
                                << " , " << std::boolalpha << domain->periodic_boundary_z() << std::endl;                    
            lout << "alpha   = "<<p.alpha << std::endl;
            lout << "radius  = "<<p.radius << std::endl;
            lout << "epsilon = "<<p.epsilon << std::endl;
            lout << "kmax    = "<<p.kmax << std::endl;
            lout << "nk      = "<<p.nk << std::endl;
            lout << "nknz    = "<<p.nknz << std::endl;
            lout << "bt      = "<<p.bt << std::endl;
            lout << "gm      = "<<p.gm << std::endl;
            lout << "lm      = "<< p.lm.x <<','<<p.lm.y<<','<<p.lm.z << std::endl;
            lout << "volume  = "<< p.volume << std::endl;
            lout << "fpe0    = "<< fpe0 << std::endl;
            lout << "GnMax   = "<< p.GnMax << std::endl;
            lout << "=================================" << std::endl;
          }
        }
      }

      *rcut = *radius;
      *rcut_max = std::max( *rcut_max , *radius );
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {  
    OperatorNodeFactory::instance()->register_factory( "ewald_init" , make_simple_operator< EwaldInitOperator > );
  }

}


