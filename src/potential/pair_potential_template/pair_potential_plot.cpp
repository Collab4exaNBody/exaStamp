#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/log.h>

#include <fstream>

#include "potential.h"

#define CLASS_NAME USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_plot)

namespace exaStamp
{
  using namespace exanb;

  struct CLASS_NAME : public OperatorNode
  {
    ADD_SLOT( USTAMP_POTENTIAL_PARAMS , parameters , INPUT , REQUIRED );
    ADD_SLOT( double                  , rcut       , INPUT , REQUIRED );
    ADD_SLOT( long                    , samples    , INPUT , 1000 );
    ADD_SLOT( std::string             , file       , INPUT , std::string( std::string(USTAMP_STR(USTAMP_POTENTIAL_NAME)) + std::string(".csv") ) );
    ADD_SLOT( ParticleSpecies         , species    , INPUT , REQUIRED );
    ADD_SLOT( std::string             , type_a     , INPUT );
    ADD_SLOT( std::string             , type_b     , INPUT );
 
    inline void execute () override final
    {
      const USTAMP_POTENTIAL_PARAMS p = *parameters;
      const ParticleSpecies& species = *(this->species);

      size_t specy_index_a = 0;
      size_t specy_index_b = 0;
      if( species.size()!=1 && ! ( type_a.has_value() && type_b.has_value() ) )
      {
        lerr<<"Error: missing type_a and/or type_b inputs"<<std::endl;
        std::abort();
      }
  
      if( type_a.has_value() && type_b.has_value() )
      {
        std::string specy_name_a = *type_a;
        std::string specy_name_b = *type_b;
        for(size_t s=0;s<species.size();s++)
        {
          if( species[s].m_name == specy_name_a ) { specy_index_a = s; }
          if( species[s].m_name == specy_name_b ) { specy_index_b = s; }
        }
        ldbg << "specy_name_a="<<specy_name_a<< ", specy_index_a = "<<specy_index_a<<std::endl;
        ldbg << "specy_name_b="<<specy_name_b<< ", specy_index_b = "<<specy_index_b<<std::endl;
      }

      double min_e = std::numeric_limits<double>::max();
      double min_e_r = -1.0;

      // gather required parameters
      const PairPotentialParameters pair = { species[specy_index_a] , species[specy_index_b] };
      double rc = *rcut;
      
      size_t n = *samples;
      std::ofstream fout( *file );
      double step = rc / n;
      double small_step = step*1.e-2;

      for(size_t i=1;i<n;i++)
      {       
        double r = i * step;        
        double e = 0.0, de = 0.0;
        double e1 = 0.0;
        double e2 = 0.0;
        double approx_de = 0.0;
        USTAMP_POTENTIAL_COMPUTE( p, pair, r-small_step, e1, approx_de );
        USTAMP_POTENTIAL_COMPUTE( p, pair, r+small_step, e2, approx_de );
        USTAMP_POTENTIAL_COMPUTE( p, pair, r, e, de );
        if( e < min_e ) { min_e=e; min_e_r=r; }
        approx_de = (e2-e1) / (2.0*small_step);
        double approx_diff = approx_de - de;
        double ampl = std::max( std::abs(approx_de) , std::abs(de) );
        if( ampl > 0.0 ) { approx_diff /= ampl; }
        fout << r << " " << e << " " << de << " " << approx_de << " " << approx_diff << std::endl;
      }
      fout << std::flush;
      fout.close();
      
      ldbg << "Minimum Energy : " << min_e << " at r="<<min_e_r<<std::endl;
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( USTAMP_STR(CLASS_NAME) , make_compatible_operator< CLASS_NAME > );
  }

}

#undef CLASS_NAME
#undef CLASS_NAME
#undef OPERATOR_NAME_STR

