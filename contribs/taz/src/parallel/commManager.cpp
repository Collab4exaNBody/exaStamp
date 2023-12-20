/// @file 
/// @brief Implementations for the communication manager


#include <vector>

#include "domain/domainInfo.hpp"

#include "parallel/commManager.hpp"
#include "parallel/communications/session.hpp"


/// @brief Constructor
/// @param [in] comm MPI communicator
CommManager::CommManager(MPI_Comm comm) 
  : comm(comm), numberOfNodes(0), rank(-1),
    domains(), domainToNode(),
    sessions() {
	// Get rank and number of nodes from the communicator
  this->updateCommInfo();
}


/// @brief Destructor
CommManager::~CommManager() {

	// Delete all sessions
  for (auto& elem : sessions) {
    if (elem.second!=nullptr) 
      delete elem.second;
  }

  // And clear the sessions map
  sessions.clear();

}


/// @brief Link the domains to their nodes
/// @param [in] domainsPerNode Index of the domains on each node
void CommManager::setDomainToNode(const std::vector< std::vector<int> >& domainsPerNode) {

  domainToNode = Array<int>(Global::domainInfo.getNumberOfDomains());

  for (int node=0; node<numberOfNodes; node++) {
    for (auto domain : domainsPerNode[node]) 
      domainToNode[domain] = node;
  }

  domains.resize(domainsPerNode[rank].size(), nullptr);

}


/// @brief Accessor to the MPI communicator
const MPI_Comm& CommManager::getCommunicator() const {
  return comm;
}


/// @brief Accessor to the number of nodes
int CommManager::getNumberOfNodes() const {
  return numberOfNodes;
}


/// @brief Accessor to the rank of the current node
int CommManager::getRank() const {
  return rank;
}


/// @brief Get the node where the specified domain stands
/// @param [in] domain Domain index
/// @return node index
int CommManager::getNodeOfDomain(int domain) {
  return domainToNode[domain];
}


// @bief Get the number of domains on current node
/// @return Number of domains
int CommManager::getNumberOfDomainsOnNode() {
  return domains.size();
}


/// @brief Update the communication manager with data from the MPI communicator
void CommManager::updateCommInfo() {
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &numberOfNodes);
}


/// @brief Synchronize all nodes
void CommManager::barrier() const {
  MPI_Barrier(comm);
}


/// @brief Close (and desallocate) a session of communications
/// @param [in] id Session index
void CommManager::closeSession(int id) {

  delete sessions[id];
  sessions[id] = nullptr;
  sessions.erase(id);

}
