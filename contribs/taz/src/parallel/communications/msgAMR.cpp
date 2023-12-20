#ifdef GHOST_V2

/// @file 
/// @brief Implementations for the messages

/// @brief Constructor
/// @tparam T Type of data in the message
/// @param [in] index Index of the domain
/// @param [in] neighbors Indexes of the neighbors domains
 msgSendAMR::msgSendAMR(int index, const Array<uint>& neighbors)
  : domainIndex(index), data() {

	// For each neighbor domain
  for (uint i=0; i<neighbors.size(); ++i) {

  	// If it's not the current domain, its index is positive and it's not in the data map
    if (neighbors[i]==(uint)domainIndex || (int)neighbors[i]<0 || data.find(neighbors[i])!=data.end()) 
      continue;

    // Add an empty data vector for this index in the data map
    data[neighbors[i]] = msgAMR();

  }
  
}


/// @brief Destructor
///
/// Clear the data map
/// @tparam T Type of data in the message
msgSendAMR::~msgSendAMR() {
  data.clear();
}


/// @brief Copy constructor
/// @tparam T Type of data in the constructed message
/// @tparam U Type of data in the copied message
/// @param [in] msg Message to copy
msgSendAMR::msgSendAMR(const msgSendAMR & msg)
  : domainIndex(msg.getDomainIndex()), data() {

	// For each domain in the map of the other message
  for (auto it=msg.cbegin(); it!=msg.cend(); ++it) {

  	// Add a empty data vector for this index in the data map
#ifndef __MIC__
    data.emplace((*it).first, msgAMR());
#else
    data.insert(std::make_pair((*it).first, msgAMR()));
#endif

  }
  
}


/// @brief Acessor to the domain index
/// @tparam T Type of data in the message
inline int msgSendAMR::getDomainIndex() const {
  return domainIndex;
}


/// @brief Get the total size
/// @tparam T Type of data in the message
/// @return Total number of objects to send
inline uint msgSendAMR::getTotalSize() {

  uint n=0;

  // For each vector in the data map
  for (auto& elem : data) 

  	// Add its size to the count
    n += elem.second.size();

  // Return that count
  return n;

}


/// @brief Get the size of the message to send to a specified domain
/// @tparam T Type of data in the message
/// @param [in] domain Index of the recipient domain
/// @return Number of objects to send to that domain
inline uint msgSendAMR::getSize(int domain) {

	// Find the specified domain in the data map
  auto it = data.find(domain);

  // Return the size of the data vector if there is some data
  if (it!=data.end()) 
    return it->second.size();

  // Or return 0
  else 
    return 0;

}


/// @brief Get a pointer to the data for the specified domain
/// @tparam T Type of data in the message
/// @param [in] domain Index of the recipient domain
/// @return Pointer to the data to send to that domain
inline msgAMR* msgSendAMR::getData(int domain) {

	// Find the specified domain in the data map
  auto it = data.find(domain);

  // Return a pointer to the data if there is some
  if (it!=data.end()) 
    return it->second.data();

  // Or return nullptr
  else 
    return nullptr;

}


/// @brief Add an object to the message for a domain
/// @tparam T Type of data in the message
/// @param [in] domain Index of the recipient domain
/// @param [in] object Object to add
inline void msgSendAMR::push(int domain, const msgAMR& object) {

	// Find the specified domain in the data map
  auto it = data.find(domain);

  // If the domain was found, add the object to the corresponding data vector
  if (it!=data.end()) 
    it->second.push_back(object);

}


/// @brief Clear the message
/// @tparam Type of data in the message
inline void msgSendAMR::clear() {
  // For each vector in the data map
  for (auto& elem : data) 
  	// Clear that vector
    elem.second.clear();
}


/// @brief Go to the start of the message data
/// @tparam T Type of data in the message
/// @return Iterator on the start message data
inline typename msgSendAMR::iterator msgSendAMR::begin() {
  return data.begin();
}


/// @brief Go to the end of the message data
/// @tparam T Type of data in the message
/// @return Iterator on the end message data
inline typename msgSendAMR::iterator msgSendAMR::end() {
  return data.end();
}


/// @brief Go to the start of the message data
/// @tparam T Type of data in the message
/// @return Constant iterator on the start message data
inline typename msgSendAMR::const_iterator msgSendAMR::cbegin() const {
  return data.cbegin();
}


/// @brief Go to the end of the message data
/// @tparam T Type of data in the message
/// @return Constant iterator on the end message data
inline typename msgSendAMR::const_iterator msgSendAMR::cend() const {
  return data.cend();
}





//
//  MESSAGE RECV
//


/// @brief Constructor
/// @tparam Type of data in the message
/// @param [in] index Index of the current domain
msgRecvAMR ::msgRecvAMR (int index)
  : domainIndex(index), data() {
}


/// @brief Destructor (nothing to do)
msgRecvAMR ::~msgRecvAMR () {}


/// @brief Copy constructor
/// @tparam T Type of data in the constructed message
/// @tparam U Type of data in the copied message
/// @param [in] msg Copied message
msgRecvAMR ::msgRecvAMR (const msgRecvAMR & msg)
  : domainIndex(msg.getDomainIndex()), data() {
}


/// @brief Accessor to the domain index
inline int msgRecvAMR ::getDomainIndex() const {
  return domainIndex;
}


/// @brief Get the size of the message data
/// @tparam T Type of data in the message
/// @return Size of the data vector
inline uint msgRecvAMR ::getSize() {
  return data.size();
}


/// @brief Get a pointer to the message data
/// @tparam T Type of data in the message
/// @return Data pointer of the data vector
inline msgAMR* msgRecvAMR ::getData() {
  return &data;
}


/// @brief Resize the message data
/// @tparam T Type of data in the message
/// @param [in] size New size
inline void msgRecvAMR ::setSize(uint size) {
  data.resize(size);
}


/// @brief Clear the message data
/// @tparam T Type of data in the message
inline void msgRecvAMR ::clear() {
  data.clear();
}

#endif
