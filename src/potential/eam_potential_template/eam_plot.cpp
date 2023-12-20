#include <exanb/core/basic_types.h>
#include <exanb/core/basic_types_operators.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/log.h>
#include <exanb/core/cpp_utils.h>

#include <iostream>

#include "potential.h"

#define EamPotentialName USTAMP_CONCAT(USTAMP_POTENTIAL_NAME,_plot)

#define EamPotentialOperatorName EamPotentialName
#define EamPotentialStr USTAMP_STR(EamPotentialName)

namespace exaStamp
{
  using namespace exanb;

  class EamPotentialOperatorName : public OperatorNode
  {  
    // ========= I/O slots =======================
    ADD_SLOT( USTAMP_POTENTIAL_PARMS, parameters        , INPUT );
    ADD_SLOT( double                , rcut              , INPUT );
    ADD_SLOT( long                  , samples    , INPUT , 1000 );
    ADD_SLOT( std::string           , file       , INPUT , std::string( EamPotentialStr ".csv" ) );
    
   public:
    // Operator execution
    inline void execute () override final
    {
      double phiCut = 0.;
      double rhoCut = 0.;
#     ifdef USTAMP_POTENTIAL_RCUT
      {
        double Rho = 0.;
        double dRho = 0.;
        double Phi = 0.;
        double dPhi = 0.;
        USTAMP_POTENTIAL_EAM_RHO( *parameters, *rcut, Rho, dRho );          
        USTAMP_POTENTIAL_EAM_PHI( *parameters, *rcut, Phi, dPhi );
        phiCut = Phi;
        rhoCut = Rho;
      }
#     endif

      double rc = *rcut;      
      size_t n = *samples;
      std::ofstream fout( *file );
      double step = rc / n;

      for(size_t i=1;i<n;i++)
      {       
        double r = i * step;        
        double Rho = 0.;
        double dRho = 0.;
        double Phi = 0.;
        double dPhi = 0.;
        double F = 0.;
        double dF = 0.;

        USTAMP_POTENTIAL_EAM_RHO( *parameters, r, Rho, dRho ); Rho -= rhoCut;          
        USTAMP_POTENTIAL_EAM_PHI( *parameters, r, Phi, dPhi ); Phi -= phiCut;
        USTAMP_POTENTIAL_EAM_EMB( *parameters, Rho, F, dF );
        fout << "r="<<r << " rho=" <<  Rho << " drho=" << dRho << " phi=" <<  Phi << " dphi=" << dPhi << " F=" << F << " dF=" << dF << std::endl;
      }
      fout << std::flush;
      fout.close();
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( EamPotentialStr , make_simple_operator< EamPotentialOperatorName > );
  }

}

