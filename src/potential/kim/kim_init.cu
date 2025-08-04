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

#include <onika/math/basic_types.h>
#include <onika/math/basic_types_operators.h>
#include <onika/math/basic_types_stream.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_factory.h>
#include <onika/scg/operator_slot.h>
#include <exanb/core/domain.h>
#include <onika/log.h>
#include <onika/cpp_utils.h>
#include "kim.h"
#define MY_ERROR(message)                                                \
  {                                                                      \
    std::cout << "* Error : \"" << message << "\" : " << __LINE__ << ":" \
              << __FILE__ << std::endl;                                  \
    exit(1);                                                             \
  }

namespace exaStamp
{

  using namespace exanb;
  using namespace KIM;
  
  class KIMInitOperator : public OperatorNode
  {
    // ========= I/O slots =======================
    //    ADD_SLOT( KIMParams,   parameters    , INPUT , REQUIRED );
    ADD_SLOT( std::string,   kim_model_name , INPUT_OUTPUT , REQUIRED );
    ADD_SLOT( double,        rcut_max       , INPUT_OUTPUT );    
    ADD_SLOT( KIMContext,    kim_ctx        , OUTPUT );
    
  public:
    // Operator execution
    inline void execute () override final
    {
      std::cout << "KIM Initialization function" << std::endl;

      KIM::Model * kim_local_model;
      int requestedUnitsAccepted;
      int error = KIM::Model::Create(KIM::NUMBERING::zeroBased,
                                     KIM::LENGTH_UNIT::A,
                                     KIM::ENERGY_UNIT::eV,
                                     KIM::CHARGE_UNIT::e,
                                     KIM::TEMPERATURE_UNIT::K,
                                     KIM::TIME_UNIT::ps,
                                     *kim_model_name,
                                     &requestedUnitsAccepted,
                                     &kim_local_model);
      if (error) { MY_ERROR("KIM::Model::Create()"); }
      
      // Check for compatibility with the model
      if (!requestedUnitsAccepted) { MY_ERROR("Must Adapt to model units"); }

      // Check that we know about all required routines
      int numberOfModelRoutineNames;
      KIM::MODEL_ROUTINE_NAME::GetNumberOfModelRoutineNames(&numberOfModelRoutineNames);
      
      for (int i = 0; i < numberOfModelRoutineNames; ++i)
        {
          KIM::ModelRoutineName modelRoutineName;
          int error
            = KIM::MODEL_ROUTINE_NAME::GetModelRoutineName(i, &modelRoutineName);
          if (error) { MY_ERROR("Unable to get ModelRoutineName."); }
          int present;
          int required;
          error = kim_local_model->IsRoutinePresent(modelRoutineName, &present, &required);
          if (error) { MY_ERROR("Unable to get routine present/required."); }

          std::cout << "Model routine name \"" << modelRoutineName.ToString()
                    << "\" has present = " << present
                    << " and required = " << required << "." << std::endl;

          if ((present == true) && (required == true))
            {
              using namespace KIM::MODEL_ROUTINE_NAME;
              if (!((modelRoutineName == Create)
                    || (modelRoutineName == ComputeArgumentsCreate)
                    || (modelRoutineName == Compute) || (modelRoutineName == Refresh)
                    || (modelRoutineName == ComputeArgumentsDestroy)
                    || (modelRoutineName == Destroy)))
                {
                  MY_ERROR("Unknown Routine \"" + modelRoutineName.ToString()
                           + "\" is required by model.");
                }
            }
        }

      // print model units
      KIM::LengthUnit lengthUnit;
      KIM::EnergyUnit energyUnit;
      KIM::ChargeUnit chargeUnit;
      KIM::TemperatureUnit temperatureUnit;
      KIM::TimeUnit timeUnit;

      kim_local_model->GetUnits(&lengthUnit, &energyUnit, &chargeUnit, &temperatureUnit, &timeUnit);

      std::cout << "LengthUnit is \"" << lengthUnit.ToString() << "\"" << std::endl
                << "EnergyUnit is \"" << energyUnit.ToString() << "\"" << std::endl
                << "ChargeUnit is \"" << chargeUnit.ToString() << "\"" << std::endl
                << "TemperatureUnit is \"" << temperatureUnit.ToString() << "\""
                << std::endl
                << "TimeUnit is \"" << timeUnit.ToString() << "\"" << std::endl;

      // check species
      int speciesIsSupported;
      int modelTaCode;
      error = kim_local_model->GetSpeciesSupportAndCode(KIM::SPECIES_NAME::Ta, &speciesIsSupported, &modelTaCode);
      if ((error) || (!speciesIsSupported))
        {
          MY_ERROR("Species Ta not supported");
        }

      KIM::ComputeArguments * computeArguments;
      error = kim_local_model->ComputeArgumentsCreate(&computeArguments);
      if (error) { MY_ERROR("Unable to create a ComputeArguments object."); }

      const double* cutoffs;
      double influencedistance;
      
      const int* modelWillNotRequestNeighborsOfNoncontributingParticles;
      int numberOfNeighborLists;
      kim_local_model->GetNeighborListPointers(&numberOfNeighborLists,
                                               &cutoffs,
                                               &modelWillNotRequestNeighborsOfNoncontributingParticles);
      double rcutmax = 0.;
      for (int i=0; i<numberOfNeighborLists;i++) {
        rcutmax = std::max(rcutmax, cutoffs[i]);
      }
      *rcut_max = rcutmax;
      std::cout << "RCUT MAX = " << *rcut_max << std::endl;
      kim_local_model->GetInfluenceDistance( &influencedistance);
      std::cout << "INF DIST = " << influencedistance << std::endl;
      
    }
    
    //    (*parameters).rcut = *rcut;    
  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(kim_init)
  {  
    OperatorNodeFactory::instance()->register_factory( "kim_init" , make_simple_operator< KIMInitOperator > );
  }

}
