#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
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

  class EwaldPotentialEnergyShiftOperator : public OperatorNode
  {
    // ========= I/O slots =======================
    ADD_SLOT( EwaldParms , ewald_config            , INPUT , REQUIRED );
    ADD_SLOT( double     , sum_square_charge       , INPUT , REQUIRED );
    ADD_SLOT( bool       , trigger_thermo_state    , INPUT , OPTIONAL );
    ADD_SLOT( EwaldRho   , ewald_rho               , INPUT , REQUIRED );
    ADD_SLOT( double     , potential_energy_shift  , OUTPUT );

  public:
    // Operator execution
    inline void execute () override final
    {
      using ewald_constants::fpe0;
      // compute auto interaction
      // fpe0 = 4piEpsilon0

      bool need_energy = false;
      if( trigger_thermo_state.has_value() ) need_energy = *trigger_thermo_state ;

      if( need_energy )
      {
        ldbg << "------------------------------------------------"<<std::endl<<std::flush;
        ldbg << "Calculation of the auto-interaction contribution"<<std::endl<<std::flush;
        ldbg << "------------------------------------------------"<<std::endl<<std::flush;
        ldbg << "sum_square_charge = " << (*sum_square_charge) <<std::endl<<std::flush;
        ldbg << "alpha = " << ewald_config->alpha <<std::endl<<std::flush;

        *potential_energy_shift = - 1. / fpe0 * ewald_config->alpha / std::sqrt(M_PI) * (*sum_square_charge);

        // reciprocal enegy
        const size_t nk = ewald_rho->nk;
        double re = 0.;
  #     pragma omp parallel for schedule(static) reduction(+:re)
        for (size_t k=0; k<nk; ++k)
        {
          re +=  ewald_config->Gdata[k].Gc * complex_norm( ewald_rho->rho[k] ); // (totalRho_r[k]*totalRho_r[k] + totalRho_i[k]*totalRho_i[k]);
        }
        //re /= total_particles;

        ldbg << "potential_energy_shift = " << (*potential_energy_shift) <<std::endl;
        ldbg << "reciprocal enegy = " << re  <<std::endl;
        *potential_energy_shift += re;
        ldbg << "total energy shift = " << (*potential_energy_shift) <<std::endl;
      }
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {  
    OperatorNodeFactory::instance()->register_factory( "ewald_potential_energy_shift" , make_simple_operator< EwaldPotentialEnergyShiftOperator > );
  }

}


