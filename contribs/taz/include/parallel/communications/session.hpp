/// @file 
/// @brief Declaration of the session class

#ifndef __SESSSION_HPP_INCLUDED
#define __SESSSION_HPP_INCLUDED


#include <map>
#include <set>

#include "parallel/commManager.hpp"
#include "parallel/communications/communication.hpp"
#include "parallel/communications/message.hpp"


/// @brief Base for the class Session
class ISession {

public:

  /// @brief Enumeration of the types of communication handled by the program
  ///
  /// If you want to add a new type of communication you must add a type to that enumeration
  /// and declare the associated message center
  enum Type {
  	PARTICLE_EXCHANGE, ///< Exchange of moving particles at each step after the position update
  	GHOST_UPDATE, ///< Update of the ghost at each step before the force calculation
  	EAM_EMB, ///< Exchange of the EAM data at each step in the middle of the force calculation
	UPDATE_DENSITIES, ///< Exchange of the local densities
	GHOST_UPDATE_VELOCITIES, ///< Update the ghosts' velocities
	GHOST_UPDATE_FLUCTUATION, ///< Update the ghosts' fluctuation forces
  	PARTICLE_SYNC, ///< Exchange of the moved particles after the load balancing
  	CELL_SYNC, ///< Exchange of the moved cells neighbors after the load balancing
  	GHOST_OWNER ///< Exchange of the owner of the ghost cells after the load balancing
  };

  /// @brief Constructor
  /// @param [in] commManager_ Pointer to the communication manager
  ISession(CommManager* commManager_)
    : commManager(commManager_), 
      numberOfDomainsOnNode(commManager->getNumberOfDomainsOnNode()) {
  }

  /// @brief Destructor (nothing to do)
  virtual ~ISession() {}

protected:
  
  CommManager* commManager; ///< Pointer to the communication manager
  uint numberOfDomainsOnNode; ///< Number of domains on the node

};


/// @brief Manage the point-to-point communications
/// @tparam T Type of data exchanged on the session
template <class T> class Session : public ISession {

public:

  Session(CommManager* commManager_);
  ~Session();

  void pushMessages(MessageSend<T>* s);
  void collectMessages(MessageRecv<T>* r);

protected:

  void execute();

  void packMessages();
  void exchangeMessages();

  void unpackMessages();

  void initCommunicationMap(MPI_Request*& req, MPI_Status*& stat);
  void fillSendBuffers();
  void allocRecvBuffers();

  void processRecvBuffers();
  void writeRecvMessages();

  void wait();

  void cleanRequestAndStatus();


private:


  // map domains <-> messages
  std::map<int, MessageSend<T>*> toSend; ///< Map of pointers to the message to send for each domain of the node
  std::map<int, MessageRecv<T>*> toRecv; ///< Map of pointers to the message to receive for each domain of the node

  // map node to communicates to reformated messages
  std::map<int, Communication<T>* > communications; ///< Map of the communication for each recipient node

  uint numberOfEffectiveCommunications; ///< Count the MPI communications
  MPI_Request* request; ///< MPI requests for all the communications
  MPI_Status* status; ///< MPI statuses for all the communications

};


/// @brief Constructor
/// @tparam T Type of data exchanged on the session
/// @param [in] commManager_ Communication manager
template <class T> Session<T>::Session(CommManager* commManager_)
  : ISession(commManager_),
    request(nullptr), 
    status(nullptr) {
}


/// @brief Destructor
/// @tparam T Type of data exchanged on the session
template <class T> Session<T>::~Session() {

	// Unlink the messages received
  for (auto& elem : toRecv) elem.second = nullptr;
  // Delete the communication
  for (auto& elem : communications) delete elem.second;

  //Clear the maps
  toSend.clear();
  toRecv.clear();
  communications.clear();

  // Clean the requests and statuses
  this->cleanRequestAndStatus();

}


/// @brief Send a message
///
/// Link the message to the session and send the message when all the messages of the node are there
/// @tparam T Type of data exchanged on the session
/// @param [in] s Pointer to where the new message is
template <class T> void Session<T>::pushMessages(MessageSend<T>* s) {

	// Get the domain of the new message
  int sendDomainIndex = s->getDomainIndex(); 

  // Add the pointer to the map if it's a new domain
  if (toSend.find(sendDomainIndex)==toSend.end()) 
    toSend[sendDomainIndex] = s;

  // If there is a message for all the domains of the node execute the session (send the messages)
  if (toSend.size()==numberOfDomainsOnNode)
    this->execute();
  
}


/// @brief Collect a message
///
/// Link the message to the session and collect the message when all the messages of the node are there
/// @tparam T Type of data exchanged on the session
/// @param [in] r Pointer to the where message to receive is
template <class T> void Session<T>::collectMessages(MessageRecv<T>* r) {

	// Get the domain of the collected message
	int recvDomainIndex = r->getDomainIndex();

  // Add the pointer to the map if it's a new domain
  if (toRecv.find(recvDomainIndex)==toRecv.end()) 
    toRecv[recvDomainIndex] = r;

  // If there is a message for all the domains of the node collect the data
  if (toRecv.size()==numberOfDomainsOnNode) 
    this->unpackMessages();

}


/// @brief Setup and start the communications
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::execute() {

	// Setup everything for the communications
  this->packMessages();

  // Exhange the data
  this->exchangeMessages();

}


/// @brief Setup communications from the messages
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::packMessages() {

	// Init the communications
  this->initCommunicationMap(request, status);

  // Put the messages data into the communications
  this->fillSendBuffers();

  // Wait for all exchanges to be finished (exchanges of domains data and content sizes)
  this->wait();

  // Clear current requests and statuses
  this->cleanRequestAndStatus();

  // Configure the incoming communications
  this->allocRecvBuffers();

}


/// @brief Finish the communications and rebuild messages
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::unpackMessages() {

	// Configure the messages  to receive
  this->processRecvBuffers();

  // Wait for all exchanges to be finished
  this->wait();

  // Clear current requests and statuses
  this->cleanRequestAndStatus();


  this->writeRecvMessages();

}


/// @brief Initialize the communications
/// @tparam T Type of data exchanged on the session
/// @param [out] req Pointer to the MPI requests for all communications
/// @param [out] stat r to the MPI statuses for all communications
template <class T> void Session<T>::initCommunicationMap(MPI_Request*& req, MPI_Status*& stat) {

	// Recipient node data
  std::map<int, int> totalSizePerNode; // Number of domains per node
  std::map<int, std::set<int> > destDomains; // Index of the domains per node
  
  // Loop on all messages sent by domains : create node map and compute total sent size per node
  //
  // For each message to send (it is an iterator on toSend, first is a sending domain and second is the message for this domain)
  for (auto it=toSend.begin(); it!=toSend.end(); ++it) {
  	// For each recipient domain (jt is an iterator on the previous message, first is a recipient domain and second is the data for this domain)
    for (auto jt=it->second->begin(); jt!=it->second->end(); ++jt) {
      
    	// Get the node of the recipient domain
      int destNode = commManager->getNodeOfDomain(jt->first);

      // If it's a new domain
      if (totalSizePerNode.find(destNode)==totalSizePerNode.end()) {
      	// Get the size of the data associated to the recipient domain as the size for the recipient node
      	totalSizePerNode[destNode] = jt->second.size();
      	// Add the recipient domain to the domain of the recipient node
      	destDomains[destNode].insert(jt->first);
      }
      // Else
      else {
      	// Add the size of the data associated to the recipient domain to the size for the recipient node
      	totalSizePerNode[destNode] += jt->second.size();
      	// Add the recipient domain to the domain of the recipient node
      	destDomains[destNode].insert(jt->first);
      }
      
    }
  }

  // Loop on node map : create communications structures
  //
  // For each node where data will be send
  for (auto it=totalSizePerNode.begin(); it!=totalSizePerNode.end(); ++it) {
    int destNode = it->first;
    // Build a new communication with the recipient domain and the total data size for that node
    Communication<T>* tmp = new Communication<T>(commManager->getCommunicator(), commManager->getRank(),destDomains[destNode], it->second);
    // And put it into the communications map
    communications.insert(std::pair<int, Communication<T>*>(destNode, tmp)); 
  }

  // Loop on all messages sent by domains : compute size per recipient domain
  //
  // For each message to send (it is an iterator on toSend, first is a sending domain and second is the message for this domain)
  for (auto it=toSend.begin(); it!=toSend.end(); ++it) {
  	// For each recipient domain (jt is an iterator on the previous message, first is a recipient domain and second is the data for this domain)
  	for (auto jt=it->second->begin(); jt!=it->second->end(); ++jt) {
      int destNode = commManager->getNodeOfDomain(jt->first);
      // Add the size of the message data to the size of the communication content for this domain
      communications[destNode]->getSendContent()->getSize(jt->first) += jt->second.size();
    }
  }

  // There are two MPI communications (one send and one receive) for each recipient node except if the recipient
  // node is the send node in which case there is non communication
  numberOfEffectiveCommunications = 2*destDomains.size();
  if (destDomains.find(commManager->getRank())!=destDomains.end()) numberOfEffectiveCommunications -= 2;
  // Get one request and one status for each communication
  req = new MPI_Request [numberOfEffectiveCommunications];
  stat = new MPI_Status [numberOfEffectiveCommunications];

  int count=0;
  // For each communication
  for (auto it=communications.begin(); it!=communications.end(); ++it) {
  	// Set the position for each domain in the node exchanged content
    it->second->getSendContent()->setStart();
    // Exchange size and domains data to setup the communication
    it->second->sendRecvSizes(it->first, &req[2*count], numberOfDomainsOnNode);
    count++;
  }

}


/// @brief Put the messages data into the communication content
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::fillSendBuffers() {

  // Loop on all messages sent by domains : copy data do buffer
  //
  //   it : iterator on toSend (first is a sending domain | second is the message to be send)
  //   jt : iterator on it->second (first is the dest domain | second is the container containing data)
  //

  // For each message to send (it is an iterator on toSend, first is a sending domain and second is the message for this domain)
  for (auto it : toSend) {
  	// For each recipient domain (jt is an iterator on the previous message, first is a recipient domain and second is the data for this domain)
    for (auto jt=it.second->begin(); jt!=it.second->end(); ++jt) {

    	// Get the node for that domain
      int destNode = commManager->getNodeOfDomain(jt->first);
      // Get a pointer to the content of the communication for that node
      Content<T>* tmp = communications[destNode]->getSendContent();
      // Get the data of the message into the content of the communication
      tmp->pack(jt->first, jt->second.data(), jt->second.size());

    }
  }

  // Clear the messages send
  for (auto& message : toSend) {
    message.second->clear();
    message.second = nullptr;
  }
  // And the map that point to them
  toSend.clear();

}


/// @brief Allocate and configure the communications incoming contents from the received size data
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::allocRecvBuffers() {

	for (auto it : communications)
    it.second->allocRecvBuffers(it.first);

}


/// @brief Exchange the data
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::exchangeMessages() {

	// Set new requests and statuses
  request = new MPI_Request [numberOfEffectiveCommunications];
  status = new MPI_Status [numberOfEffectiveCommunications];

  // Exchange the data for each communication
  int count=0;
  for (auto it : communications) {
    it.second->exchangeMessages(it.first, &request[2*count]);
    count++;
  }

}


/// @brief Setup messages to receive
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::processRecvBuffers() {

  std::map<int, int> totalSizePerDomain; // Size of the received message for each domain

  // For each message to received (it is an iterator on toRecv, first is a recipient domain and second is the message for this domain)
  for (auto& elem : toRecv)
  	// Set a 0 size to the message received
    totalSizePerDomain.insert(std::pair<int,int>(elem.first, 0));

  // Get total recv size per domain
 
  // For each communication
  for (auto& it : communications) {

  	// Get a pointer to the content
    Content<T>* tmp = it.second->getRecvContent();

    // For each domain in that content
    for (int i=0; i<tmp->getNumberOfDomains(); ++i)
    	// Add the size for the domain in that content to the total size for the domain
      totalSizePerDomain[tmp->getDomainI(i)] += tmp->getSizeI(i);

  }

  // Set the size for each message to receive
  for (auto& it : toRecv)
    it.second->setSize(totalSizePerDomain[it.first]);

}


/// @brief Get the data from the communications contents to the received messages
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::writeRecvMessages() {

  std::map<int, int> nextPerDomain; // Position to write in each message
  for (auto& elem : toRecv) 
    nextPerDomain.insert(std::pair<int,int>(elem.first, 0));

  // For each communication
  for (auto& it : communications) {

  	// Get a pointer to the content
    Content<T>* tmp = it.second->getRecvContent();

    // For each domain in this content
    for (int i=0; i<tmp->getNumberOfDomains(); ++i) {
    	// Get the domain index
      int destDomain = tmp->getDomainI(i);
      // And the data size
      int destSize = tmp->getSizeI(i);
      // Get data from the content to the message for this domain
      tmp->unpack(toRecv[destDomain]->getData() + nextPerDomain[destDomain], i);
      // Update the writing position
      nextPerDomain[destDomain] += destSize;
    }

  }

}


/// @brief Wait for all the communication to finish
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::wait() {
  MPI_Waitall(numberOfEffectiveCommunications, request, status);
}


/// @brief Delete current requests and statuses
/// @tparam T Type of data exchanged on the session
template <class T> void Session<T>::cleanRequestAndStatus() {

  if (request != nullptr) delete [] request;
  if (status  != nullptr) delete [] status;

  request = nullptr;
  status  = nullptr;

}

#endif // __SESSSION_HPP_INCLUDED
