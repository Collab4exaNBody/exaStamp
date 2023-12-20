/// @file 
/// @brief Implementations for the messages

/// @brief Constructor
/// @tparam T Type of data in the message
/// @param [in] index Index of the domain
/// @param [in] neighbors Indexes of the neighbors domains
template <class T>  MessageSend<T>::MessageSend(int index, const Array<uint>& neighbors)
  : domainIndex(index), data() {

	// For each neighbor domain
  for (uint i=0; i<neighbors.size(); ++i) {

  	// If it's not the current domain, its index is positive and it's not in the data map
    if (neighbors[i]==(uint)domainIndex || (int)neighbors[i]<0 || data.find(neighbors[i])!=data.end()) 
      continue;

    // Add an empty data vector for this index in the data map
    data[neighbors[i]] = std::vector<T>();

  }
  
}


/// @brief Destructor
///
/// Clear the data map
/// @tparam T Type of data in the message
template <class T> MessageSend<T>::~MessageSend() {
  data.clear();
}


/// @brief Copy constructor
/// @tparam T Type of data in the constructed message
/// @tparam U Type of data in the copied message
/// @param [in] msg Message to copy
template <class T> template <class U> MessageSend<T>::MessageSend(const MessageSend<U>& msg)
  : domainIndex(msg.getDomainIndex()), data() {

	// For each domain in the map of the other message
  for (auto it=msg.cbegin(); it!=msg.cend(); ++it) {

  	// Add a empty data vector for this index in the data map
#ifndef __MIC__
    data.emplace((*it).first, std::vector<T>());
#else
    data.insert(std::make_pair((*it).first, std::vector<T>()));
#endif

  }
  
}


/// @brief Acessor to the domain index
/// @tparam T Type of data in the message
template <class T> inline int MessageSend<T>::getDomainIndex() const {
  return domainIndex;
}


/// @brief Get the total size
/// @tparam T Type of data in the message
/// @return Total number of objects to send
template <class T> inline uint MessageSend<T>::getTotalSize() {

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
template <class T> inline uint MessageSend<T>::getSize(int domain) {

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
template <class T> inline T* MessageSend<T>::getData(int domain) {

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
template <class T> inline void MessageSend<T>::push(int domain, const T& object) {

	// Find the specified domain in the data map
  auto it = data.find(domain);

  // If the domain was found, add the object to the corresponding data vector
  if (it!=data.end()) 
    it->second.push_back(object);

}

/// @brief Add an object to the message for a domain
/// @tparam T Type of data in the message
/// @param [in] domain Index of the recipient domain
/// @param [in] object Object to add
template <class T> inline void MessageSend<T>::push(int domain, const T* object, const int n) {

	// Find the specified domain in the data map
  auto it = data.find(domain);

  // If the domain was found, add the object to the corresponding data vector
  if (it!=data.end()) 
  {
    int size = it->second.size();
    it->second.resize(n+size);
    std::copy(object, object+n, it->second.data()+size);
  }
    it->second.push_back(object);

}


/// @brief Clear the message
/// @tparam Type of data in the message
template <class T> inline void MessageSend<T>::clear() {
  // For each vector in the data map
  for (auto& elem : data) 
  	// Clear that vector
    elem.second.clear();
}


/// @brief Go to the start of the message data
/// @tparam T Type of data in the message
/// @return Iterator on the start message data
template <class T> inline typename MessageSend<T>::iterator MessageSend<T>::begin() {
  return data.begin();
}


/// @brief Go to the end of the message data
/// @tparam T Type of data in the message
/// @return Iterator on the end message data
template <class T>
inline typename MessageSend<T>::iterator MessageSend<T>::end() {
  return data.end();
}


/// @brief Go to the start of the message data
/// @tparam T Type of data in the message
/// @return Constant iterator on the start message data
template <class T>
inline typename MessageSend<T>::const_iterator MessageSend<T>::cbegin() const {
  return data.cbegin();
}


/// @brief Go to the end of the message data
/// @tparam T Type of data in the message
/// @return Constant iterator on the end message data
template <class T>
inline typename MessageSend<T>::const_iterator MessageSend<T>::cend() const {
  return data.cend();
}


//
//  MESSAGE RECV
//


/// @brief Constructor
/// @tparam Type of data in the message
/// @param [in] index Index of the current domain
template <class T> MessageRecv<T>::MessageRecv(int index)
  : domainIndex(index), data(0) {
}


/// @brief Destructor (nothing to do)
template <class T> MessageRecv<T>::~MessageRecv() {}


/// @brief Copy constructor
/// @tparam T Type of data in the constructed message
/// @tparam U Type of data in the copied message
/// @param [in] msg Copied message
template <class T> template <class U> MessageRecv<T>::MessageRecv(const MessageRecv<U>& msg)
  : domainIndex(msg.getDomainIndex()), data(0) {
}


/// @brief Accessor to the domain index
template <class T> inline int MessageRecv<T>::getDomainIndex() const {
  return domainIndex;
}


/// @brief Get the size of the message data
/// @tparam T Type of data in the message
/// @return Size of the data vector
template <class T> inline uint MessageRecv<T>::getSize() {
  return data.size();
}


/// @brief Get a pointer to the message data
/// @tparam T Type of data in the message
/// @return Data pointer of the data vector
template <class T> inline T* MessageRecv<T>::getData() {
  return data.data();
}


/// @brief Resize the message data
/// @tparam T Type of data in the message
/// @param [in] size New size
template <class T> inline void MessageRecv<T>::setSize(uint size) {
  data.resize(size);
}


/// @brief Clear the message data
/// @tparam T Type of data in the message
template <class T> inline void MessageRecv<T>::clear() {
  data.clear();
}
