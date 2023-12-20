/// @file 
/// @brief Class to handle main time loop

#ifndef __TIME_HPP_INCLUDED
#define __TIME_HPP_INCLUDED


#include "io/input.hpp"


class TimeManager;


/// @brief Temporary structure gathering all data to initialize the physical time manager
template <> struct Configuration<TimeManager> {

  /// @brief Default constructor
  Configuration() {}
  /// @brief Destructor (nothing to do)
  ~Configuration() {}

  /// @brief Constructor from input data
  /// @param [in] input Input to convert
  Configuration(const Input& input)
    : start(0.),
      step(0),
      delta(input.delta),
      numberOfSteps(input.numberOfSteps),
      logRate(input.logRate),
      balanceRate(input.balanceRate),
      outputRate(input.outputRate),
      dumpRate(input.dumpRate) {
    if (input.legacyHeader!=nullptr) {
      start = input.legacyHeader->time;
      step  = input.legacyHeader->iterationNumber;
    }
    if (input.stopWallsTime > 0)
      stopWallsTime = input.stopWallsTime;
    else
      stopWallsTime = start+delta*numberOfSteps*2;
  }

  double start; ///< Start time
  uint   step; ///< Start step

  double delta; ///< Time step

  uint numberOfSteps; ///< Max number of steps
  int logRate; ///< Number of steps between each log writing
  int balanceRate; ///< Number of steps between each load balancing
  int outputRate; ///< Number of steps between each output writing
  int dumpRate; ///< Number of steps between each dump writing

  double stopWallsTime; ///< Time when walls should be stopped
};


/// @brief Manager for the physical time
class TimeManager {

public:

  /// @brief Constructor with a configuration
  /// @param [in] configuration Configuration to convert
  TimeManager(Configuration<TimeManager>& configuration) 
    : m_start(configuration.start), 
      m_final(m_start+configuration.delta*configuration.numberOfSteps),
      m_delta(configuration.delta),
      m_elapsed(m_start),
      m_currentStep(configuration.step), 
      m_numberOfSteps(configuration.numberOfSteps),
      m_logRate(configuration.logRate),
      m_balanceRate(configuration.balanceRate),
      m_dumpRate(configuration.dumpRate),
      m_outputRate(configuration.outputRate),
      m_stopWallsTime(configuration.stopWallsTime),
      m_outOfTime(false),
      m_runningWalls(true) {
    Global::seed += m_currentStep;
  }

  /// @brief Destructor (nothing to do)
  ~TimeManager() {}

  TimeManager& operator ++ ();
  TimeManager  operator ++ (int);

  /// @brief True if the calculation is finished
  ///
  /// Max number of steps exceeded or out of time
  bool isFinished() { return (m_currentStep>=m_numberOfSteps) || m_outOfTime; }

  /// @brief True if time for log writings
  bool log() { return m_currentStep % m_logRate == 0; }
  /// @brief True if time for output writings
  bool output() { return m_currentStep % m_outputRate == 0; }
  /// @brief True if time for load balancing
  bool balance() { return m_currentStep % m_balanceRate == 0; }
  /// @brief True if time for dump writings
  bool dump() { return m_currentStep % m_dumpRate == 0; }

  /// @brief Accessor to m_start
  double start() { return m_start; }
  /// @brief Accessor to m_final
  double final() { return m_final; }
  /// @brief Accessor to m_delta
  double delta() { return m_delta; }

  /// @brief Accessor to m_elapsed
  double elapsed() { return m_elapsed; }

  /// @brief Accessor to m_currentStep
  uint currentStep() { return m_currentStep; }
  /// @brief Accessor to m_numberOfSteps
  uint numberOfSteps() { return m_numberOfSteps; }

  /// @brief Accessor to m_logRate
  int logRate() { return m_logRate; }
  /// @brief Accessor to m_outputRate
  int outputRate() { return m_outputRate; }
  /// @brief Accessor to m_balanceRate
  int balanceRate() { return m_balanceRate; }
  /// @brief Accessor to m_dumpRate
  int dumpRate() { return m_dumpRate; }
  /// @brief Accessor to outOfTime
  bool outOfTime() { return m_outOfTime; }

  /// @brief True if time has reached m_stopWallsTime
  bool stopWalls() {
    if (m_runningWalls) {
      m_runningWalls = !(m_elapsed>=m_stopWallsTime);
      return !m_runningWalls;
    }
    return false;
  }
  
  /// @brief Set m_outOfTime at true to abort calculation
  ///
  ///
  void abort() { m_outOfTime=true; }

private:

  TimeManager& operator += (int t);

  double m_start; ///< Start time
  double m_final; ///< Final time
  double m_delta; ///< Time step
  double m_elapsed; ///< Elapsed time

  uint m_currentStep; ///< Number of the current step
  uint m_numberOfSteps; ///< Total number of steps

  int m_logRate; ///< Number of steps between each log writing
  int m_balanceRate; ///< Number of steps between each load balancing
  int m_dumpRate; ///< Number of steps between each dump writing
  int m_outputRate; ///< Number of steps between each output writing

  double m_stopWallsTime; ///< Time when walls should be stopped
  
  bool m_outOfTime; ///< Stop condition
  bool m_runningWalls; ///< False when walls have been stopped
  
};


/// @brief Increment operator
/// @param [in] t Time to add
inline TimeManager& TimeManager::operator += (int t) {
  m_currentStep+=t;
  m_elapsed += t*m_delta;
  return *this;
}

/// @brief Pre increment operator
///
///
inline TimeManager& TimeManager::operator ++ () {
  m_currentStep++;
  m_elapsed += m_delta;
  return *this;
}

/// @brief Post increment operator
///
///
inline TimeManager TimeManager::operator ++ (int) {
  TimeManager clone(*this);
  m_currentStep++;
  m_elapsed += m_delta;
  return clone;
}

#endif // __TIME_HPP_INCLUDED
