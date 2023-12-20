/// @file
/// @brief Definition of the Metrics class

#ifndef __METRICS_HPP_INCLUDED
#define __METRICS_HPP_INCLUDED



#include <sys/resource.h>

#include "parallel/mympi.hpp"
#include "parallel/commManager.hpp"

#include "utils/chrono.hpp"

#include "metricsDetails.hpp"



/// @brief A class to measure time and memory usage
///
/// Measure and store time used for the different program parts of the program :
/// initialization, computation, I/O and load balancing
class Metrics {
public:

  /// @brief Enumeration of the possible times to measure
  enum Time {
	  INIT,             ///< Initialization time
	  COMPUTE,          ///< Computation time
	  IO,               ///< I/O time
	  BALANCE,          ///< Load balancing time
    TIMELOOP,         ///< Measure time loop duration
    SIMULATIONSLEEP,  ///< Measure when the simulation is sleeping
    COPY,             ///< Measure the time to copy data
    ORCHESTRATORSLEEP,///< Measure when the orchestrator is sleeping
    ANALYTICS,        ///< Measure analytics duration
    NEIGHBOURS,       ///< Measure the time to build neighbours lists
    POTENTIAL,        ///< Measure the time to do potential kernel
    GHOST,            ///< Measure the time to MPI Comms
    REFINE,           ///< Measure the time to function linked to the refinement
    GENERIC           ///< For new timer
  };
  

  /// @brief Constructor
  Metrics() 
    : m_chrono_init(), m_chrono_compute(), m_chrono_io(), m_chrono_balance(), m_chrono_timeloop(),
      m_chrono_simulationsleep(), m_chrono_copy(), m_chrono_orchestratorsleep(), m_chrono_analytics(), 
      m_chrono_neighbours(), m_chrono_potential(), m_chrono_ghost(), m_chrono_refine(), m_chrono_generic(), 
      m_time_init(0.), m_time_compute(0.), m_time_io(0.), m_time_balance(0.), m_time_timeloop(0.), m_time_simulationsleep(0.), m_time_copy(0.), m_time_orchestratorsleep(0.), m_time_analytics(0.),
      m_time_neighbours(0), m_time_potential(0), m_time_ghost(0), m_time_refine(0), m_time_generic(0),
      data() {}

  /// @brief Destructor (nothing to do)
  ~Metrics() {}
  
  /// @brief Start chosen chrono
  /// @param [in] time Type of the measured time
  void tic (const Time& time) {
    switch (time) {
    case INIT   : m_chrono_init.tic();    break;
    case COMPUTE: m_chrono_compute.tic(); break;
    case IO     : m_chrono_io.tic();      break;
    case BALANCE: m_chrono_balance.tic(); break;
    case TIMELOOP: m_chrono_timeloop.tic(); break;
    case SIMULATIONSLEEP  : m_chrono_simulationsleep.tic()  ; break;
    case COPY             : m_chrono_copy.tic()             ; break;
    case ORCHESTRATORSLEEP: m_chrono_orchestratorsleep.tic(); break;
    case ANALYTICS        : m_chrono_analytics.tic()        ; break;
    case NEIGHBOURS : m_chrono_neighbours.tic() ; break;
    case POTENTIAL  : m_chrono_potential.tic()  ; break;
    case GHOST      : m_chrono_ghost.tic()      ; break;
    case REFINE     : m_chrono_refine.tic()     ; break;
    case GENERIC    : m_chrono_generic.tic()    ; break;
    }
  }
  /// @brief Stop chosen chrono
  /// @param [in] time Type of the measured time
  void toc (const Time& time) {
    switch (time) {
    case INIT   : m_chrono_init.toc();    break;
    case COMPUTE: m_chrono_compute.toc(); break;
    case IO     : m_chrono_io.toc();      break;   
    case BALANCE: m_chrono_balance.toc(); break;
    case TIMELOOP: m_chrono_timeloop.toc(); break;
    case SIMULATIONSLEEP  : m_chrono_simulationsleep.toc()  ; break;
    case COPY             : m_chrono_copy.toc()             ; break;
    case ORCHESTRATORSLEEP: m_chrono_orchestratorsleep.toc(); break;
    case ANALYTICS        : m_chrono_analytics.toc()        ; break;
    case NEIGHBOURS: m_chrono_neighbours.toc() ; break;
    case POTENTIAL : m_chrono_potential.toc()  ; break;
    case GHOST     : m_chrono_ghost.toc()      ; break;
    case REFINE    : m_chrono_refine.toc()     ; break;
    case GENERIC   : m_chrono_generic.toc()    ; break;
    }
  }
  /// @brief Store time measured by chosen chrono
  /// @param [in] time Type of the measured time
  void store (const Time& time) {
    switch (time) {
    case INIT   : m_time_init   +=m_chrono_init.elapsed();    break;
    case COMPUTE: m_time_compute+=m_chrono_compute.elapsed(); break;
    case IO     : m_time_io     +=m_chrono_io.elapsed();      break;
    case BALANCE: m_time_balance+=m_chrono_balance.elapsed(); break;
    case TIMELOOP         : m_time_timeloop     +=m_chrono_timeloop.elapsed()             ; break;
    case SIMULATIONSLEEP  : m_time_simulationsleep+=m_chrono_simulationsleep.elapsed()    ; break;
    case COPY             : m_time_copy+=m_chrono_copy.elapsed()                          ; break;
    case ORCHESTRATORSLEEP: m_time_orchestratorsleep+=m_chrono_orchestratorsleep.elapsed(); break;
    case ANALYTICS        : m_time_analytics+=m_chrono_analytics.elapsed()                ; break;
    case NEIGHBOURS : m_time_neighbours += m_chrono_neighbours.elapsed() ; data.fill_m_time_neighbours(m_chrono_neighbours.elapsed()) ; break;
    case POTENTIAL  : m_time_potential  += m_chrono_potential.elapsed()  ; data.fill_m_time_potential(m_chrono_potential.elapsed())   ; break;
    case GHOST      : m_time_ghost      += m_chrono_ghost.elapsed()      ; data.fill_m_time_ghost(m_chrono_ghost.elapsed())           ; break;
    case REFINE     : m_time_refine     += m_chrono_refine.elapsed()     ; data.fill_m_time_refine(m_chrono_refine.elapsed())         ; break;
    case GENERIC    : m_time_generic    += m_chrono_generic.elapsed()    ; data.fill_m_time_generic(m_chrono_generic.elapsed())       ; break;
    }
  }
  /// @brief Measure time on chosen chrono
  /// @param [in] time Type of the measured time
  /// @return Measured time
  double elapsed (const Time& time) {
    switch (time) {
    case INIT   : return m_chrono_init.elapsed();    break;
    case COMPUTE: return m_chrono_compute.elapsed(); break;
    case IO     : return m_chrono_io.elapsed();      break;
    case BALANCE: return m_chrono_balance.elapsed(); break;
    case TIMELOOP          : return m_chrono_timeloop.elapsed()         ; break;
    case SIMULATIONSLEEP   : return m_chrono_simulationsleep.elapsed()  ; break;
    case COPY              : return m_chrono_copy.elapsed()             ; break;
    case ORCHESTRATORSLEEP : return m_chrono_orchestratorsleep.elapsed(); break;
    case ANALYTICS  : return m_chrono_analytics.elapsed()  ; break;
    case NEIGHBOURS : return m_chrono_neighbours.elapsed() ; break;
    case POTENTIAL  : return m_chrono_potential.elapsed()  ; break;
    case GHOST      : return m_chrono_ghost.elapsed()      ; break;
    case REFINE     : return m_chrono_refine.elapsed()     ; break;
    case GENERIC    : return m_chrono_generic.elapsed()    ; break;
    }
    return 0;
  }
  /// @brief Return total time for chosen program part
  /// @param [in] time Type of the measured time
  /// @return Total time
  double total (const Time& time) {
    switch (time) {
    case INIT   : return m_time_init;    break;
    case COMPUTE: return m_time_compute; break;
    case IO     : return m_time_io;      break;
    case BALANCE: return m_time_balance; break;
    case TIMELOOP          : return m_time_timeloop          ; break;
    case SIMULATIONSLEEP   : return m_time_simulationsleep   ; break;
    case COPY              : return m_time_copy              ; break;
    case ORCHESTRATORSLEEP : return m_time_orchestratorsleep ; break;
    case ANALYTICS:   return m_time_analytics; break;
    case NEIGHBOURS : return m_time_neighbours ; break;
    case POTENTIAL :  return m_time_potential  ; break;
    case GHOST :      return m_time_ghost      ; break;
    case REFINE :     return m_time_refine     ; break;
    case GENERIC :    return m_time_generic    ; break;
    }
    return 0;
  }
  /// @brief Reset total time for chosen program part
  /// @param [in] time Type of the measured time
  void reset (const Time& time) {
    switch (time) {
    case INIT   : m_time_init=0.;    break;
    case COMPUTE: m_time_compute=0.; break;
    case IO     : m_time_io=0.;      break;
    case BALANCE: m_time_balance=0.; break;
    case TIMELOOP          : m_time_timeloop=0.         ; break;
    case SIMULATIONSLEEP   : m_time_simulationsleep=0.  ; break;
    case COPY              : m_time_copy=0.             ; break;
    case ORCHESTRATORSLEEP : m_time_orchestratorsleep=0.; break;
    case ANALYTICS         : m_time_analytics=0.; break;
    case NEIGHBOURS : m_time_neighbours = 0.; break;
    case POTENTIAL  : m_time_potential  = 0.; break;
    case GHOST      : m_time_ghost      = 0.; break;
    case REFINE     : m_time_refine     = 0.; break;
    case GENERIC    : m_time_generic    = 0.; break;
    }
  }

  /// @brief Get total memory usage
  /// @param [in] comm The communication manager
  /// @return Memory usage (kB)
  double memStats(CommManager& comm) {

    rusage _rusage;
    Array<double> localRssUsg(1, 0.);
    Array<double> totalRssUsg(1, 0.);
    
    getrusage(RUSAGE_SELF, &_rusage);
    localRssUsg[0] = (double) _rusage.ru_maxrss;
    comm.reduce(localRssUsg, totalRssUsg, MPI_SUM);

    return totalRssUsg[0];

  }

private:

  Chrono m_chrono_init;              ///< A chrono to measure time used for initialization
  Chrono m_chrono_compute;           ///< A chrono to measure time used for computation
  Chrono m_chrono_io;                ///< A chrono to measure time used for I/O
  Chrono m_chrono_balance;           ///< A chrono to measure time used for load balancing
  Chrono m_chrono_timeloop;          ///< A chrono to measure time used for time loop
  Chrono m_chrono_simulationsleep;   ///< A chrono to measure time used for simulation sleeping
  Chrono m_chrono_copy;              ///< A chrono to measure time used for data copy
  Chrono m_chrono_orchestratorsleep; ///< A chrono to measure time used for orchestrator sleeping
  Chrono m_chrono_analytics;         ///< A chrono to measure time used for analytics
  Chrono m_chrono_neighbours;        ///< A chrono to measure time used for neighbours
  Chrono m_chrono_potential;         ///< A chrono to measure time used for potential
  Chrono m_chrono_ghost;             ///< A chrono to measure time used for ghost
  Chrono m_chrono_refine;            ///< A chrono to measure time used for refine
  Chrono m_chrono_generic;           ///< A chrono to measure time used for generic

  double m_time_init;                ///< Store total time used for initialization
  double m_time_compute;             ///< Store total time used for computation
  double m_time_io;                  ///< Store total time used for I/O
  double m_time_balance;             ///< Store total time used for load balancing
  double m_time_timeloop;            ///< Store total time used for I/O
  double m_time_simulationsleep;     ///< Store total time used for simulation sleeping
  double m_time_copy;                ///< Store total time used for data copy
  double m_time_orchestratorsleep;   ///< Store total time used for orchestrator sleeping
  double m_time_analytics;           ///< Store total time used for analytics
  double m_time_neighbours;          ///< Store total time used for neighbours
  double m_time_potential;           ///< Store total time used for potential
  double m_time_ghost;               ///< Store total time used for ghost
  double m_time_refine;              ///< Store total time used for refine
  double m_time_generic;             ///< Store total time used for generic

  public :
  metricsDetails data;               ///< Each optional timers are stored for each iteration

};



#endif //  __METRICS_HPP_INCLUDED
