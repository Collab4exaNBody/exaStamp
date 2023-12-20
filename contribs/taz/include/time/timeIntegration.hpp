/// @file 
/// @brief Interface for a time-integrator and three implementations
///
/// @warning All implementations are in NVE or NVT ensemble. There is no other NPT
/// scheme available yet

#ifndef __SCHEMES_HPP_INCLUDED
#define __SCHEMES_HPP_INCLUDED


#include "utils/stampUnits.hpp"

#include "time/computeForceStrategy.hpp"


class DomainInterface;
class Input;


/// @brief Integrates particles trajectories (in phase space)
///
/// This object takes particles at a time \f$ t \f$ 
/// and compute their new state (position, velocity, ...) at time 
/// \f$ t+\Delta t \f$.
///
class NumericalScheme {
public:
  /// @brief Constructor
  /// @param [in] recovering Enable or disable the communication recovering
  NumericalScheme(bool recovering) {
    if (recovering)
      forceComputer = new ComputeForceRecovering();
    else
      forceComputer = new ComputeForceNoRecovering();
  }
  /// @brief Destructor
  virtual ~NumericalScheme() {
    if (forceComputer != nullptr) delete forceComputer;
  }

  /// @brief Do one step of the integration scheme for a domain
  /// @param [in,out] domain Domain where the step is done
  /// @param [in] time Time step size
  virtual void oneStep(DomainInterface* domain, double time) = 0;
  ComputeForceStrategy* forceComputer;  ///< strategy to compute forces
};


/// @brief Temporary structure to extract parameters from input in order to build an
/// integration scheme for the system
template <> struct Configuration<NumericalScheme> {

  // Parameters
  /// @brief Enumeration of the types of numerical schemes
  enum Type {
	  NEWTON, ///< Scheme based on Newton movement equations
	  LANGEVIN, ///< Fixed temperature scheme based on Langevin dynamics
  };
  /// @brief Enumeration of the subtypes of numerical schemes
  enum Subtype {
	  VERLET_LEAPFROG, ///< Leapfrog integration scheme
	  VERLET_VELOCITY, ///< Verlet Velocity integration scheme
	  SPLITTING, ///< Langevin Splitting integration scheme
  };

  /// @brief Default constructor
  Configuration() {}

  /// @brief Destructor (nothing to do)
  ~Configuration() {}

  Configuration(const Input& input);

  Type type;       ///< Type of Scheme to build
  Subtype subtype; ///< Subtype ...

  double frictionLangevin; ///< Friction term for a Langevin scheme
  double temperatureThermostat; ///< Target temperature for thermostats

  bool recovering; ///< Allow communication recovering
};


// Create a NumericalScheme from a configuration structure
NumericalScheme* createScheme(Configuration<NumericalScheme>& configuration);


/// @brief Leapfrog Verlet integration scheme
///
/// This class is used to compute particles' new state using the
/// Leapfrog Verlet algorithm
/// \f[ \left\{ \begin{array}{ccl} x^{t+^1/_2} &=&  x^t + \frac{\Delta t}{2} v^t \\ \mathtt{compute} && a^{t+^1/_2} \\ v^{t+1} &=& v^t + \Delta t a^{t+^1/_2} \\ x^{t+1} &=&  x^{t+^1/_2} + \frac{\Delta t}{2} v^{t+1} \end{array} \right. \f]
///
class NewtonVerletLeapfrog : public NumericalScheme {
public:
  /// @brief Default constructor
  ///
  /// @param [in] recovering Enable the communication recovering
  NewtonVerletLeapfrog(bool recovering) : NumericalScheme(recovering) {}
  virtual void oneStep(DomainInterface* domain, double time);
};


/// @brief Leapfrog Verlet integration scheme (version 2)
///
/// This class is used to compute particles' new state using the 
/// Velocity Verlet algorithm
/// \f[ \left\{ \begin{array}{ccl} x^{t+1} &=&  x^t + \Delta t v^t + \frac{\Delta t ^2}{2} a^t \\ v^{t+^1/_2} &=& v^t + \frac{\Delta t}{2} a^t \\ \mathtt{compute} && a^{t+1} \\ v^{t+1} &=& v^{t+^1/_2} + \frac{\Delta t}{2} a^{t+1} \end{array} \right. \f]
///
class NewtonVerletVelocity : public NumericalScheme {
public:
  /// @brief Default constructor
  ///
  /// @param [in] recovering Enable the communication recovering
  NewtonVerletVelocity(bool recovering) : NumericalScheme(recovering) {}
  virtual void oneStep(DomainInterface* domain, double time);
};


/// @brief Langevin integration scheme
///
/// This class is used to compute particles' new state using the 
/// Langevin thermostat
class LangevinSplitting : public NumericalScheme {

private:

  double gamma; ///< Friction parameter for the Langevin scheme
  double beta; ///< 1/kB*T

public:

  /// @brief Default constructor
  ///
  /// @param [in] recovering Enable the communication recovering
  /// @param [in] friction Friction term
  /// @param [in] temperature System temperature
  LangevinSplitting(bool recovering, double friction, double temperature)
    : NumericalScheme(recovering), gamma(friction), beta(1.00/(Stamp_Constant::boltzmann * temperature)) {}

  virtual void oneStep(DomainInterface* domain, double time);
};



#endif // __TIME_INTEGRATION_HPP_INCLUDED
