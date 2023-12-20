/// @file
/// @brief Implementations for the load balancer


#include <iostream>

#include "parallel/balance.hpp"

#include "domain/domainInterface.hpp"


/// @brief Initialize the load balancer
/// @param [in] argc Program's arguments number.
/// @param [in] argv Program's arguments.
void LBS::init(int argc, char** argv) {

#if __use_lib_zoltan
  float version;
  Zoltan_Initialize(argc, argv, &version);
#endif

}


/// @brief Constructor
/// @param [in] comm MPI communicator for all nodes
LBS::LoadBalancer::LoadBalancer(MPI_Comm comm) 
      : m_method(NONE), m_zoltanStruct(nullptr) {

#if __use_lib_zoltan
      m_zoltanStruct = Zoltan_Create(comm);
#endif

}


/// @brief Destructor (nothing to do)
LBS::LoadBalancer::~LoadBalancer() {
}


/// @brief Ill named destructor for the load balancer
void LBS::LoadBalancer::destroy() {

#if __use_lib_zoltan
      Zoltan_Destroy(&m_zoltanStruct);
#endif

      m_zoltanStruct = nullptr;

}


/// @brief Set the parameters for the load balancer
void LBS::LoadBalancer::setParameters() {

#if __use_lib_zoltan

  Zoltan_Set_Param(m_zoltanStruct, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param(m_zoltanStruct, "NUM_LID_ENTRIES", "1");
  Zoltan_Set_Param(m_zoltanStruct, "DEBUG_LEVEL", "0");

  Zoltan_Set_Param(m_zoltanStruct, "IMBALANCE_TOL", "1.05");

  Zoltan_Set_Param(m_zoltanStruct, "OBJ_WEIGHT_DIM", "2"); 
  Zoltan_Set_Param(m_zoltanStruct, "EDGE_WEIGHT_DIM", "1"); 

  Zoltan_Set_Param(m_zoltanStruct, "RETURN_LISTS", "ALL");

  Zoltan_Set_Param(m_zoltanStruct, "REDUCE_DIMENSIONS", "1");
  // Zoltan_Set_Param(m_zoltanStruct, "DEGENERATE_RATIO", "10");
  // Zoltan_Set_Param(m_zoltanStruct, "RCB_MULTICRITERIA_NORM", "2");

  Zoltan_Set_Param(m_zoltanStruct, "KEEP_CUTS", "1");

#endif

}


/// @brief Set the method for the load balancer
/// @param [in] method_ Load balancing method
/// @param [in] graphPartition Load balancing approach
void LBS::LoadBalancer::setMethod(Method method_, GraphPartition graphPartition) {

  m_method = method_;

#if __use_lib_zoltan

  bool graph = false;

  switch (m_method) {

  case NONE:
    Zoltan_Set_Param(m_zoltanStruct, "LB_METHOD", "NONE");
    break;
    
  case BLOCK:
    Zoltan_Set_Param(m_zoltanStruct, "LB_METHOD", "BLOCK");
    break;

  case RANDOM:
    Zoltan_Set_Param(m_zoltanStruct, "LB_METHOD", "RANDOM");
    break;

  case REC_COORD_BSC:
    Zoltan_Set_Param(m_zoltanStruct, "LB_METHOD", "RCB");
    break;

  case REC_INERT_BSC:
    Zoltan_Set_Param(m_zoltanStruct, "LB_METHOD", "RIB");
    break;

  case SPACE_FILL_CURVE:
    Zoltan_Set_Param(m_zoltanStruct, "LB_METHOD", "HSFC");
    break;

  // In the graph cases set GRAPH_PACKAGE instead of LB_METHOD
  case PHG:
    graph = true;
    Zoltan_Set_Param(m_zoltanStruct, "GRAPH_PACKAGE", "PHG");
    break;

#ifdef __use_lib_metis
  case PARMETIS:
    graph = true;
    Zoltan_Set_Param(m_zoltanStruct, "GRAPH_PACKAGE", "PARMETIS");
    break;
#endif

#ifdef __use_lib_scotch
  case SCOTCH:
    graph = true;
    Zoltan_Set_Param(m_zoltanStruct, "GRAPH_PACKAGE", "SCOTCH");
    break;
#endif

  }

  // In the graph cases, set LB_METHOD to GRAPH and set the load balancing approach
  if (graph) {
    
    Zoltan_Set_Param(m_zoltanStruct, "LB_METHOD", "GRAPH");

    switch (graphPartition) {

    case PARTITION :
      Zoltan_Set_Param(m_zoltanStruct, "LB_APPROACH", "PARTITION");
      break;
    case REPARTITION :
      Zoltan_Set_Param(m_zoltanStruct, "LB_APPROACH", "REPARTITION");
      break;
    case REFINE :
      Zoltan_Set_Param(m_zoltanStruct, "LB_APPROACH", "REFINE");
      break;

    }

  }

#endif

}


/// @brief Set the functions used for the load balancing
/// @param [in] numberOfDomains Number of domain on the node
/// @param [in] domains Pointer to the domains on the node
void LBS::LoadBalancer::setCallbackQueryFunctions(int numberOfDomains, DomainInterface** domains) {

#if __use_lib_zoltan

  if (numberOfDomains!=1) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'LBS::setCallbackQueryFunctions(Zoltan_Struct*, int, DomainInterface*)' : Only one domain per node is currently handled. STOP." 
	     << std::endl;
    exit(-1);
  }

  domains[0]->setCallbackQueryFunctions(this);

#endif

}
