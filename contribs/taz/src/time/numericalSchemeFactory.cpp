/// @file 
/// @brief Creation of a time-integration method


#include "io/input.hpp"

#include "time/timeIntegration.hpp"


/// @brief Constructor from an Input structure
/// @param [in] input Input structure
Configuration<NumericalScheme>::Configuration(const Input& input) {

  recovering = input.recovering;

  switch (input.scheme) {

  case (Input::NEWTON_VERLET_LEAPFROG) :
    type = NEWTON;
    subtype = VERLET_LEAPFROG;
    break;

  case (Input::NEWTON_VERLET_VELOCITY) :
    type = NEWTON;
    subtype = VERLET_VELOCITY;
    break;

  case (Input::LANGEVIN) :
    type = LANGEVIN;
    subtype = SPLITTING;
    frictionLangevin = input.frictionLangevin;
    temperatureThermostat = input.temperatureThermostat;
    break;
  }

}


/// @brief Create a NumericalScheme in function of a configuration structure
/// @param [in] configuration Configuration structure to convert
/// @return Created scheme
NumericalScheme* createScheme(Configuration<NumericalScheme>& configuration) {

  NumericalScheme* integrator = nullptr;

  switch (configuration.type) {

  case Configuration<NumericalScheme>::NEWTON :

    switch (configuration.subtype) {

    case Configuration<NumericalScheme>::VERLET_LEAPFROG : 
      integrator = new NewtonVerletLeapfrog(configuration.recovering);
      break;

    case Configuration<NumericalScheme>::VERLET_VELOCITY : 
      integrator = new NewtonVerletVelocity(configuration.recovering);
      break;

    case Configuration<NumericalScheme>::SPLITTING : 
      break;

    }
    break;

  case Configuration<NumericalScheme>::LANGEVIN :
    switch (configuration.subtype) {

    case Configuration<NumericalScheme>::VERLET_LEAPFROG : 
      break;

    case Configuration<NumericalScheme>::VERLET_VELOCITY : 
      break;

    case Configuration<NumericalScheme>::SPLITTING : 
      integrator = new LangevinSplitting(configuration.recovering,configuration.frictionLangevin, configuration.temperatureThermostat);
      break;
    }
    break;
    
  }

  return integrator;

}
