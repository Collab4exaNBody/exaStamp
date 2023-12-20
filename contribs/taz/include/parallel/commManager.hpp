/// @file 
/// @brief Definition of class CommManager

#ifndef __COMM_MANAGER_HPP_INCLUDED
#define __COMM_MANAGER_HPP_INCLUDED


#include <map>
#include <vector>

#include "globals.hpp"

#include "parallel/mympi.hpp"
#include "parallel/types/MPI_typesUtils.hpp"

#include "utils/array/array.hpp"


// Some forward declarations
class DomainInterface;
class ISession;
template <class T> class Session;


/// @brief Class to handle communications between nodes
///
/// There is only one CommManager per node and all the communications
/// between go through it or its sessions
class CommManager {

public:
  CommManager(MPI_Comm comm_);

  ~CommManager();

  void setDomainToNode(const std::vector< std::vector<int> >& domainsPerNode);

  const MPI_Comm& getCommunicator() const;
  int getNumberOfNodes() const;
  int getRank() const;

  int getNodeOfDomain(int domain);
  int getNumberOfDomainsOnNode();

  void updateCommInfo();
  void barrier() const;

  template <class T> void broadcast (T& inout, int root=Global::masterNode);
  template <class T> void broadcast (Array<T>& inout, int root=Global::masterNode);
  template <class T> void gather    (Array<T>& in, Array<T>& out, int dest=Global::masterNode);
  template <class T> void gatherV   (Array<T>& in, Array<T>& out, Array<int>& count, Array<int>& disp, int dest=Global::masterNode);
  template <class T> void allGather (Array<T>& in, Array<T>& out);
  template <class T> void allGatherV(Array<T>& in, Array<T>& out, Array<int>& count, Array<int>& disp);
  template <class T> void reduce    (Array<T>& in, Array<T>& out, MPI_Op operation, int dest=Global::masterNode);
  template <class T> void reduce    (T& in, T& out, MPI_Op operation, int dest=Global::masterNode);
  template <class T> void allReduce (Array<T>& in, Array<T>& out, MPI_Op operation);

  template <class T> Session<T>* getSession(int id);
  void closeSession(int id);

private:

  MPI_Comm comm;      ///< Communicator for all nodes
  int numberOfNodes;  ///< Number of nodes on the communicator
  int rank;           ///< Rank of the current node

  std::vector<DomainInterface*> domains;	///< Vector of pointers to the domains (not used as such and full of nullptr)

  Array<int> domainToNode;  ///< Index of the node it belongs to for each domain

  std::map<int,ISession*> sessions;  ///< Map of all active sessions of communications

};


/// @brief Broadcast a variable from one node to all
/// @tparam T Type of the broadcasted variable
/// @param [in,out] inout Broadcasted variable
/// @param [in] root Index of the broadcasting node, default=master node
template <class T> 
void CommManager::broadcast(T& inout, int root) {
  MPI_Bcast(&inout, 1, MPI__Type_get<T>(), root, comm);
}


/// @brief Broadcast an array from one node to all
/// @tparam T Type of the broadcasted array elements
/// @param [in,out] inout Broadcasted array
/// @param [in] root Index of the broadcasting node, default=master node
template <class T>
void CommManager::broadcast(Array<T>& inout, int root) {
  MPI_Bcast(inout.data(), inout.size(), MPI__Type_get<T>(), root, comm);
}


/// @brief Gather elements of arrays of the same size from all nodes to one
/// @tparam T Type of the gathered elements
/// @param [in] in Array on each node
/// @param [out] out Gathered array
/// @param [in] dest Index of the gathering node, default=master node
template <class T> 
void CommManager::gather(Array<T>& in, Array<T>& out, int dest) {
  MPI_Gather(in.data(), in.size(), MPI__Type_get<T>(), out.data(), in.size(), MPI__Type_get<T>(), dest, comm);
}


/// @brief Gather elements of arrays of the varying sizes from all nodes to one
/// @tparam T Type of the gathered elements
/// @param [in] in Array on each node
/// @param [out] out Gathered array
/// @param [in] count Number of elements on each node
/// @param [in] disp Where to put the elements of each node in the final array
/// @param [in] dest Index of the gathering node, default=master node
template <class T> 
void CommManager::gatherV(Array<T>& in, Array<T>& out, Array<int>& count, Array<int>& disp, int dest) {
  MPI_Gatherv(in.data(), in.size(), MPI__Type_get<T>(), out.data(), count.data(), disp.data(), MPI__Type_get<T>(), dest, comm);
}


/// @brief Gather elements of arrays of the same size from all nodes to all
/// @tparam T Type of the gathered elements
/// @param [in] in Array on each node
/// @param [out] out Gathered array
template <class T> 
void CommManager::allGather(Array<T>& in, Array<T>& out) {
  MPI_Allgather(in.data(), in.size(), MPI__Type_get<T>(), out.data(), in.size(), MPI__Type_get<T>(), comm);
}


/// @brief Gather elements of arrays of the varying sizes from all nodes to all
/// @tparam T Type of the gathered elements
/// @param [in] in Array on each node
/// @param [out] out Gathered array
/// @param [in] count Number of elements on each node
/// @param [in] disp Where to put the elements of each node in the final array
template <class T> 
void CommManager::allGatherV(Array<T>& in, Array<T>& out, Array<int>& count, Array<int>& disp) {
  MPI_Allgatherv(in.data(), in.size(), MPI__Type_get<T>(), out.data(), count.data(), disp.data(), MPI__Type_get<T>(), comm);
}


/// @brief Reduction from all nodes to one
/// @tparam T Type of the reduced element
/// @param [in] in Data to reduce (in an array)
/// @param [out] out Reduced data (in an array)
/// @param [in] operation Reduction operation
/// @param [in] dest Index of the reducing node, default=master node
template <class T>
void CommManager::reduce(Array<T>& in, Array<T>& out, MPI_Op operation, int dest) {
  MPI_Reduce(in.data(), out.data(), in.size(), MPI__Type_get<T>(), operation, dest, comm);
}



/// @brief Reduction from all nodes to one
/// @tparam T Type of the reduced element
/// @param [in] in Data to reduce (in an array)
/// @param [out] out Reduced data (in an array)
/// @param [in] operation Reduction operation
/// @param [in] dest Index of the reducing node, default=master node
template <class T>
void CommManager::reduce(T& in, T& out, MPI_Op operation, int dest) {
  MPI_Reduce(&in, &out, 1, MPI__Type_get<T>(), operation, dest, comm);
}


/// @brief Reduction from all nodes to all
/// @tparam T Type of the reduced element
/// @param [in] in Data to reduce (in an array)
/// @param [out] out Reduced data (in an array)
/// @param [in] operation Reduction operation
template <class T>
void CommManager::allReduce(Array<T>& in, Array<T>& out, MPI_Op operation) {
  MPI_Allreduce(in.data(), out.data(), in.size(), MPI__Type_get<T>(), operation, comm);
}


/// @brief Get a session of communications
/// @tparam Type of the data to exchange
/// @param [in] id Index of the session
/// @return Session corresponding to specified index if it exists, or newly created one
template <class T> 
Session<T>* CommManager::getSession(int id) {

  if (sessions.find(id)==sessions.end()) 
    sessions[id] = new Session<T>(this);

  return static_cast<Session<T>*>(sessions[id]);
}

#endif // __COMM_MANAGER_HPP_INCLUDED
