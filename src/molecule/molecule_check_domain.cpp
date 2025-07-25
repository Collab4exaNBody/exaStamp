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
#include <exanb/core/domain.h>
#include <exanb/core/grid_algorithm.h>

#include <cmath>

namespace exaStamp
{
  using namespace exanb;

  class MoleculeCheckDomain : public OperatorNode
  {
    ADD_SLOT(Domain , domain        , INPUT , REQUIRED ); // for info printing only
    ADD_SLOT(double , bond_max_dist , INPUT , REQUIRED ); // molecule bond max distance
    ADD_SLOT(double , bond_max_dist_eps   , INPUT  , 1.e-3 ); // fraction of bond_max_dist.
    // a periodic copy of an original particle is guaranteed to be at least 2*(bond_max_dist+bond_max_eps) far away from original

  public:

    inline void execute () override final
    {    
      // domain must be at least sliplink_min_domain_size wide, so that we know how to duplicate beads/SL connectivity with periodic conditions
      double min_domain_size = ( *bond_max_dist * ( 1.0 + *bond_max_dist_eps ) ) * 2.0;
      assert( domain->bounds_size().x > min_domain_size && domain->bounds_size().y > min_domain_size && domain->bounds_size().z > min_domain_size );

      lout << "======= Molecule/Domain info ========" << std::endl;
      lout << "Domain bounds    = "<<domain->bounds()<<std::endl;
      lout << "Domain size      = "<<domain->bounds_size() <<std::endl;
      lout << "Cell size        = "<<domain->cell_size()<<std::endl;
      lout << "Grid dimensions  = "<<domain->grid_dimension()<<" ("<<grid_cell_count(domain->grid_dimension())<<" cells)"<< std::endl;
      lout << "Min. domain size = "<<min_domain_size<< std::endl;
      lout << "Bond max. dist.  = "<<(*bond_max_dist)<< std::endl;
      lout << "=====================================" << std::endl;
    }

    // -----------------------------------------------
    // -----------------------------------------------
    inline std::string documentation() const override final
    {
      return R"EOF(
Checks that domain size is greater than 2*bond_max_dist,
updates rcut_max to max(rcut_max,bond_max_dist).
)EOF";
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(molecule_check_domain)
  {
    OperatorNodeFactory::instance()->register_factory( "molecule_check_domain", make_compatible_operator< MoleculeCheckDomain > );
  }

}

