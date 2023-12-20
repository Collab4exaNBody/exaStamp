#ifdef GHOST_V2
class msgAMR {

  private :

  std::vector<uint8_t> m_ti;
  std::vector<uint64_t> m_id;
  std::vector<double> m_rx;
  std::vector<double> m_ry;
  std::vector<double> m_rz;

  public :

  inline int size() {return m_ti.size():}

  inline uint8_t ti(int i) { return m_ti[i]; }
  inline uint64_t id(int i) { return m_id[i]; }
  inline double rx(int i) { return m_rx[i]; }
  inline double ry(int i) { return m_ry[i]; }
  inline double rz(int i) { return m_rz[i]; }

  inline uint8_t* ti() { return m_ti.data(); }
  inline uint64_t* id() { return m_id.data(); }
  inline double* rx() { return m_rx.data(); }
  inline double* ry() { return m_ry.data(); }
  inline double* rz() { return m_rz.data(); }


  void resize(int newSize)
  {
    m_ti.resize(newSize);
    m_id.resize(newSize);
    m_rx.resize(newSize);
    m_ry.resize(newSize);
    m_rz.resize(newSize);
  }

  void clear()
  {
    m_ti.clear();
    m_id.clear();
    m_rx.clear();
    m_ry.clear();
    m_rz.clear();
  }
}


/// @brief Class to store what data is send to which node
/// @tparam Type of data in the message
template <class T> class msgSendAMR{

	/// @brief Shortcut for a map that index vectors of data by integer indexes
  typedef std::map<int, msgAMR> map_t;
  /// @brief Shortcut for an iterator on map_t
  typedef typename map_t::iterator iterator;
  /// @brief Shortcut for a constant iterator on map_t
  typedef typename map_t::const_iterator const_iterator;

public:

  msgSendAMR(int index, const Array<uint>& neighbors);
  ~msgSendAMR();

  msgSendAMR(const msgSendAMR &msg);

  int getDomainIndex() const ;
  
  uint getTotalSize();

  uint getSize(int domain);

  msgAMR* getData(int domain);

  void push(int domain, const T& object);

  void clear();

  iterator begin();
  iterator end();

  const_iterator cbegin() const;
  const_iterator cend() const;

protected:

  int domainIndex; ///< Index of the current domain
  map_t data; ///< Map of the data to send indexed by the domain to send to

};


class msgRecvAMR {

public:

  msgRecvAMR (int index);
  ~msgRecvAMR ();

  template <class U> msgRecvAMR (const msgRecvAMR & msg);

  int getDomainIndex() const;

  uint getSize();

  msgAMR* getData();

  void setSize(uint size);

  void clear();

protected:

  int domainIndex; ///< Index of the current domain
  msgAMR data; ///< Data received

};



/// @file 
/// @brief Definition of the message center
#include "parallel/communications/session.hpp"


#define TMPL_ MessageCenterAMR  MessageCenterAM


/// @brief Tool that handle the point-to-point node communications for a specific data type
///
/// You must create new  MessageCenterAMR each time you want to exchange a new data type
/// @tparam T Data type exchanged in this center
class  MessageCenterAMR {

public:

	/// @brief Constructor
	/// @param [in] comm_ Pointer to the communication manager
	/// @param [in] type_ Unique identifier of this communication manager in the program
   MessageCenterAMR(CommManager* comm_, ISession::Type type_)
    : comm(comm_), type(type_), 
      toSend(nullptr), toRecv(nullptr) {}

  /// @brief Destructor
	///
  /// Delete all remaining messages
  virtual ~ MessageCenterAMR() {
    if (toSend != nullptr) delete toSend;
    if (toRecv != nullptr) delete toRecv;
  }

  void init(uint domainIndex, const Array<uint>& neighbors);

  void init(const  MessageCenterAMR &  MessageCenterAMR);

  void send();
  void collect();
  void clean();

  msgSendAMR * getmsgSendAMR ();
  msgRecvAMR * getmsgRecvAMR ();
  
  const msgSendAMR & getmsgSendAMR () const ;
  const msgRecvAMR & getmsgRecvAMR () const ;
  
protected:

  CommManager* comm; ///< Pojnter to the communication manager
  ISessionAMR::Type type; ///< Unique identifier of this communication manager in the program

  msgSendAMR * toSend; ///< Message to send
  msgRecvAMR * toRecv; ///< Message received

};


/// @brief Initialize the messages
/// @param [in] domainIndex index of the domain
/// @param [in] neighbors Neighbors of the domain
inline void TMPL_ MessageCenterAMR::init(uint domainIndex, const Array<uint>& neighbors) {
  toSend = new msgSendAMR (domainIndex, neighbors);
  toRecv = new msgRecvAMR (domainIndex);
}


/// @brief Initialize the messages by copying those of another message center
/// @tparam U Data type exchanged in the other message center
/// @param [in]  MessageCenterAMR Other message center
inline void TMPL_ MessageCenterAMR::init(const  MessageCenterAMR&  MessageCenterAMR) {
  toSend = new msgSendAMR ( MessageCenterAMR.getmsgSendAMR ());
  toRecv = new msgRecvAMR ( MessageCenterAMR.getmsgRecvAMR ());
}


/// @brief Send the message to send
inline void TMPL_ MessageCenterAMR::send() {
	// Get the session corresponding to the identifier of this message center
  SessionAMR* tmp = comm->getSession(type);
  // Send the message to send on the session
  tmp->pushMessages(toSend);
}


/// @brief Collect data into the message received
inline void TMPL_ MessageCenterAMR::collect() {
	// Get the session corresponding to the identifier of this message center
  SessionAMR* tmp = comm->getSession(type);
  // Collect the data
  tmp->collectMessages(toRecv);
}


/// @brief End the session and clear the message received
inline void TMPL_ MessageCenterAMR::clean() {
  comm->closeSession(type);
  toRecv->clear();
}


/// @brief Accessor to the message to send
/// @return Pointer to the message to send
inline msgSendAMR * TMPL_ MessageCenterAMR::getmsgSendAMR () {
  return toSend;
}


/// @brief Accessor to the message received
inline msgRecvAMR * TMPL_ MessageCenterAMR::getmsgRecvAMR () {
  return toRecv;
}


/// @brief Constant accessor to the message to send
/// @return Message to send
inline const msgSendAMR & TMPL_ MessageCenterAMR::getmsgSendAMR () const {
  return *toSend;
}


/// @brief Constant accessor to the message received
/// @return Message received
inline const msgRecvAMR & TMPL_ MessageCenterAMR::getmsgRecvAMR () const {
  return *toRecv;
}


#undef TMPL_ MessageCenterAMR


// ICI --- 

#include <map>
#include <set>

#include "parallel/commManager.hpp"
#include "parallel/communications/communication.hpp"
#include "parallel/communications/message.hpp"


/// @brief Manage the point-to-point communications
/// @tparam T Type of data exchanged on the session
class SessionAMR : public ISession {

public:

  SessionAMR(CommManager* commManager_);
  ~SessionAMR();

  void pushMessages(msgSendAMR* s);
  void collectMessages(msgRecvAMR* r);

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
  std::map<int, msgSendAMR*> toSend; ///< Map of pointers to the message to send for each domain of the node
  std::map<int, msgRecvAMR*> toRecv; ///< Map of pointers to the message to receive for each domain of the node

  // map node to communicates to reformated messages
  std::map<int, Communication<T>* > communications; ///< Map of the communication for each recipient node

  uint numberOfEffectiveCommunications; ///< Count the MPI communications
  MPI_Request* request; ///< MPI requests for all the communications
  MPI_Status* status; ///< MPI statuses for all the communications

};


SessionAMR::SessionAMR(CommManager* commManager_)
  : ISession(commManager_),
    request(nullptr), 
    status(nullptr) {
}


SessionAMR::~SessionAMR() {

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


void SessionAMR::pushMessages(msgSendAMR* s) {

	// Get the domain of the new message
  int sendDomainIndex = s->getDomainIndex(); 

  // Add the pointer to the map if it's a new domain
  if (toSend.find(sendDomainIndex)==toSend.end()) 
    toSend[sendDomainIndex] = s;

  // If there is a message for all the domains of the node execute the session (send the messages)
  if (toSend.size()==numberOfDomainsOnNode)
    this->execute();
  
}


void SessionAMR::collectMessages(msgRecvAMR * r) {

	// Get the domain of the collected message
	int recvDomainIndex = r->getDomainIndex();

  // Add the pointer to the map if it's a new domain
  if (toRecv.find(recvDomainIndex)==toRecv.end()) 
    toRecv[recvDomainIndex] = r;

  // If there is a message for all the domains of the node collect the data
  if (toRecv.size()==numberOfDomainsOnNode) 
    this->unpackMessages();

}


void SessionAMR::execute() {

	// Setup everything for the communications
  this->packMessages();

  // Exhange the data
  this->exchangeMessages();

}


/// @brief Setup communications from the messages
/// @tparam T Type of data exchanged on the session
void SessionAMR::packMessages() {

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
void SessionAMR::unpackMessages() {

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
void SessionAMR::initCommunicationMap(MPI_Request*& req, MPI_Status*& stat) {

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
void SessionAMR::allocRecvBuffers() {

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
#endif


