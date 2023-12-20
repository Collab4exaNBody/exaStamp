/// @file
/// @brief Definition of load balancing related classes and functions

#ifndef __BALANCE_HPP_INCLUDED
#define __BALANCE_HPP_INCLUDED


#include <string>

#include "parallel/mympi.hpp"

#if __use_lib_zoltan
#include "zoltan.h"
#endif


class DomainInterface;


/// @brief Namespace containing load balancing related functions
namespace LBS {


#if !__use_lib_zoltan
	/// @brief If Zoltan is not used, replace the Zoltan structure by a char to avoid errors
  typedef char Zoltan_Struct;
#endif
  
  
  /// @brief Enumeration of the load balancing methods
  enum Method {
    NONE, 						///< No load balancing
    BLOCK,						///< Block load balancing (designed for testing, cannot handle the double weights of ExaStamp)
    RANDOM,						///< Random load balancing (for real; do not use)
    REC_COORD_BSC,		///< Recursive coordinate bisection (cut system in half according to the weights until the number of process is reached)
    REC_INERT_BSC,		///< Recursive inertial bisection (same as Recursive Coordinate Bisection but without axis constraint)
    SPACE_FILL_CURVE, ///< Hilbert space filling curve
    PHG,							///< Parallel hypergraph partitioning (balancing method that consider edges)
    PARMETIS,					///< Parallel graph partitioning (balancing method that consider edges)
    SCOTCH						///< Scotch graph partitioning (balancing method that consider edges)
  };


  /// @brief Enumeration of load balancing approaches (supplementary option for the graph methods)
  enum GraphPartition {
    PARTITION,		///< Static partition (used to initialize the partition)
    REPARTITION,	///< Dynamique partition (change the partition with some memory of the previous state
    REFINE				///< Fast refining of the partition
  };


  /// @brief Enumeration of the possible states for the load balancer
  enum State {
    VOID,					///< No load balancing this step
    BALANCE,			///< Load balancing done
    NOTHING_TO_DO	///< Load balancing done, no change
  };


  /// @brief Translate an load balancer state into a string
  /// @param [in] state State to translate
  /// @return Resulting string
  inline std::string str(State state) {
    switch (state) {
    case VOID:          return "";          break;
    case BALANCE:       return "balancing"; break;
    case NOTHING_TO_DO: return "unchanged"; break;
    default:            return "error";     break;
    }
    return "";
  }


  void init(int argc, char** argv);


  /// @brief Tool to balance the workload between the nodes
  class LoadBalancer {

  public:

    LoadBalancer(MPI_Comm comm);
    ~LoadBalancer();

    void destroy();

    void setParameters();
    void setMethod(Method method_, GraphPartition graphPartition=REPARTITION);
    void setCallbackQueryFunctions(int numberOfDomains, DomainInterface** domains); 

    /// @brief Accessor to the load balancing method
    Method method() { return m_method; }
    /// @brief Accessor to the load balancer main structure
    Zoltan_Struct* zoltanStruct() { return m_zoltanStruct; }

  private:

    Method m_method; ///< Load balancing method
    Zoltan_Struct* m_zoltanStruct; ///< Load balancer main structure

  };

  
}

#endif // __BALANCE_HPP_INCLUDED
