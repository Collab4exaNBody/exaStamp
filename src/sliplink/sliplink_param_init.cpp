/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exaStamp/sliplink/sliplink.h>

#include <cmath>

namespace exaStamp
{
  class SlipLinkParamInitOperator : public OperatorNode
  {
    ADD_SLOT(long               , number_of_chains    , INPUT , 1000 );
    ADD_SLOT(long               , beads_per_chain     , INPUT , 100 );
    ADD_SLOT(long               , number_of_sliplinks , INPUT , 1000 );
    ADD_SLOT(double             , beta                , INPUT , 1.0 );
    ADD_SLOT(double             , Ns                  , INPUT , 0.5 );
    ADD_SLOT(double             , kBT                 , INPUT , 1.0 );
    ADD_SLOT(double             , ksi                 , INPUT , 1.0 );
    ADD_SLOT(double             , h1                  , INPUT , 0.01 );
    ADD_SLOT(double             , sigma_rcut_scale    , INPUT , 5.0 );
    ADD_SLOT(double             , bond_max_dist       , INPUT_OUTPUT , 0.0 ); // molecule bond max distance
    ADD_SLOT(double             , dt                  , INPUT_OUTPUT , 0.0 );
    ADD_SLOT(SlipLinkParameters , sliplink_config     , OUTPUT );

  public:

    inline void execute () override final
    {
      SlipLinkParameters& conf = *sliplink_config;

      conf.number_of_chains = *number_of_chains;
      conf.beads_per_chain = *beads_per_chain;
      conf.number_of_sliplinks = *number_of_sliplinks;
      conf.beta = *beta;
      conf.Ns = *Ns;
      conf.kBT = *kBT;
      conf.ksi = *ksi;
      conf.sigma_rcut_scale = *sigma_rcut_scale;

      // timestep
      double tauzero = ( conf.ksi * conf.beta * conf.beta ) / ( 3.0 * M_PI * M_PI * conf.kBT );
      conf.h = (*h1) * tauzero;
      if( *dt == 0.0 ) { *dt = conf.h; }
      else { conf.h = *dt; }

      // computed parameters

      conf.ksis = 0.1 * conf.ksi;

      conf.sigma1 = std::sqrt( 2.0 * conf.h * conf.kBT / conf.ksi );      
      conf.sigma2 = conf.beta / std::sqrt(3.0);
      conf.sigma3 = std::sqrt(conf.Ns) * conf.beta / std::sqrt(3.0);
      conf.sigma4 = std::sqrt( ( 2.0 * conf.h * conf.kBT ) / conf.ksis );

      conf.cte2 = 1.0 / conf.ksi;
      conf.cte7 = 3.0 * conf.kBT / ( conf.beta * conf.beta );
      conf.cte8 = 3.0 * conf.kBT / ( conf.Ns * conf.beta * conf.beta );
      conf.cte4 = 3.0 * conf.kBT / ( conf.Ns * conf.beta * conf.beta * conf.ksis );

      conf.bond_max_dist = conf.sigma2 * conf.sigma_rcut_scale ;

      // generate some values to outputs slots later reused by other operators
      *bond_max_dist = std::max( *bond_max_dist , conf.bond_max_dist );

      lout << "============ SlipLinks ==============" << std::endl;
      lout << "Number of chains = "<<conf.number_of_chains<< std::endl;
      lout << "Beads per chain  = "<< conf.beads_per_chain << std::endl;
      lout << "Num. sliplinks   = "<< conf.number_of_sliplinks << std::endl;      
      lout << "Bond max. dist.  = "<< conf.bond_max_dist << std::endl;
      lout << "dt               = "<< *(dt) << std::endl;      
//#     ifndef NDEBUG
      lout << "beta             = "<< conf.beta << std::endl;      
      lout << "Ns               = "<< conf.Ns << std::endl;      
      lout << "kBT              = "<< conf.kBT << std::endl;      
      lout << "ksi              = "<< conf.ksi << std::endl;      
      lout << "h                = "<< conf.h << std::endl;      
      lout << "sigma_rcut_scale = "<< conf.sigma_rcut_scale << std::endl;      
      lout << "sigma1           = "<< conf.sigma1 << std::endl;      
      lout << "sigma2           = "<< conf.sigma2 << std::endl;      
      lout << "sigma3           = "<< conf.sigma3 << std::endl;      
      lout << "cte7             = "<< conf.cte7 << std::endl;      
      lout << "cte8             = "<< conf.cte8 << std::endl;      
      lout << "ksis             = "<< conf.ksis << std::endl;      
      lout << "cte4             = "<< conf.cte4 << std::endl;      
      lout << "sigma4           = "<< conf.sigma4 << std::endl;      
//#     endif
      lout << "=====================================" << std::endl;
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Written by Michelin & CEA/DIF
Credits to Claire Lemarchand, Ioanis Tanis, Thierry Carrard.

Initializes Slip-Link parameters and constants.
)EOF";
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(sliplink_param_init)
  {
    OperatorNodeFactory::instance()->register_factory( "sliplink_param_init", make_compatible_operator< SlipLinkParamInitOperator > );
  }

}

