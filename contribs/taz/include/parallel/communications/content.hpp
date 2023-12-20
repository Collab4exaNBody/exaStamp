/// @file 
/// @brief Definition and implementation of (template) class Content

#ifndef __CONTENT_HPP_INCLUDED
#define __CONTENT_HPP_INCLUDED


#include <set>


/// @brief Class to store the communicated data in an easy to exchange form
/// (in current node to one other node communications)
/// @tparam T Type of data exchanged
template <class T> class Content {

public:

  Content();
  ~Content();

  Content(const MPI_Comm& comm_, const std::set<int>& domains, int totalSize_);
  Content(Content<T>* sendContent);
  Content(const MPI_Comm& comm_, int nbDomains);

  void pack(int domainIndex, T* buff, int count);
  void unpack(T* buff, int i);

  int getNumberOfDomains();

  int& getStart(int domainIndex);
  int& getSize(int domainIndex);

  int& getDomainI(int i);
  int& getStartI(int i);
  int& getSizeI(int i);

  int& getNext(int domainIndex);
  T& getData(int i);

  void setStart();
  void setDataBuffer();

  int findIndex(int domainIndex);

  int getIndexesBufferSize();
  int* getIndexesBuffer();

  int getDataBufferSize();
  T* getDataBuffer();

protected:

  const MPI_Comm& comm;  ///< MPI communicator for communications

  int* __domainStartSize;  ///< One big array for the domains, start indexes of messages and sizes of messages

  int numberOfDomains;  ///< Number of domain on the recipient node
  int totalSize;        ///< Total size of field data

  int* domain;  ///< Indexes for all domains exchanging messages (pointer to the right position in __domainStartSize)
  int* start;   ///< Start position (in data) of message for each domain (pointer to the right position in __domainStartSize)
  int* size;    ///< Size of message for each domain (pointer to the right position in __domainStartSize)

  int* next;  ///< Temporary array to store the data size already added to the content for each recipient domain
  T* data;    ///< Array containing the objects to send to all the domains of the other node (in one big bloc)

};


/// @brief Default constructor
/// @tparam T Type of data exchanged
template <class T> Content<T>::Content()
  : comm(MPI_COMM_WORLD),
    __domainStartSize(nullptr),
    numberOfDomains(0), totalSize(0),
    domain(nullptr), start(nullptr), size(nullptr), 
    next(nullptr), data(nullptr) {
}


/// @brief Constructor with arguments
///
/// Used for the send contents only (the total size of the data must be known)
/// @tparam T Type of data exchanged
/// @param [in] comm_ MPI communicator used for communications
/// @param [in] domains Indexes of the domains of interest in the recipient node
/// @param [in] totalSize_ Total number of object to send to the recipient node
template <class T> Content<T>::Content(const MPI_Comm& comm_, const std::set<int>& domains, int totalSize_)
  : comm(comm_),
    __domainStartSize(nullptr),
    numberOfDomains(domains.size()), totalSize(totalSize_),
    domain(nullptr), start(nullptr), size(nullptr), 
     next(nullptr), data(nullptr) {

	// If there are domains to send data to
  if (numberOfDomains>0) {

  	// Allocate the big domains/indexes/sizes array
    __domainStartSize = new int [this->getIndexesBufferSize()];

    // Set the three pointers to the previous array
    domain = &__domainStartSize[0]; 
    start  = &__domainStartSize[1*numberOfDomains]; 
    size   = &__domainStartSize[2*numberOfDomains]; 

    // Allocate the temporary positions array
    next = new int [numberOfDomains];
    // Allocate the data array
    data = new T [totalSize];

    // Fill domains and initialize positions and sizes to 0
    int count = 0;
    for (auto elem : domains) {
      domain[count] = elem;
      start [count] = 0;
      size  [count] = 0;
      next  [count] = 0;
      count++;
    }
    
  }

}


/// @brief Copy constructor
///
/// The purpose of that constructor is to build a received content from
/// the corresponding send content. I think this won't work in the case
/// of multidomains nodes since the number of domains in the send content
/// of a communication is the number of domains on the other node while
/// the number of domains in the received content should be that of the
/// send content one the other node, that is the number of domains on the
/// current node.
/// @tparam T Type of data exchanged
/// @param [in] sendContent Content to copy
template <class T> Content<T>::Content(Content<T>* sendContent)
  : comm(sendContent->comm),
    __domainStartSize(nullptr),
    numberOfDomains(sendContent->numberOfDomains), totalSize(0), 
    domain(nullptr), start(nullptr), size(nullptr), 
    next(nullptr), data(nullptr) {
  
  if (numberOfDomains>0) {

    __domainStartSize = new int [this->getIndexesBufferSize()];

    domain = &__domainStartSize[0]; 
    start = &__domainStartSize[1*numberOfDomains]; 
    size = &__domainStartSize[2*numberOfDomains]; 

    next = new int [numberOfDomains];
    for (int i=0; i<numberOfDomains; ++i) next[i] = 0;

  }

}


/// @brief Constructor with arguments for when the total size is unknown
///
/// This constructor should be used to initialize the received content of the
/// communication by using the number of domains on the current node as argument.
///
/// Since the total size is unknown and the size data will be truncated
/// later, only the allocation of the size array is done (no allocation of the
/// data array, no initialization of the size array
/// @tparam T Type of data exchanged
/// @param [in] comm_ MPI communicator used for communications
/// @param [in] nbDomains Indexes of the domains of interest in the recipient node
template <class T> Content<T>::Content(const MPI_Comm& comm_, const int nbDomains) :
		comm(comm_),
    __domainStartSize(nullptr),
    numberOfDomains(nbDomains),
    totalSize(0),
    domain(nullptr),
    start(nullptr),
    size(nullptr),
    next(nullptr),
    data(nullptr)
    {
		// Allocate the big domains/indexes/sizes array
		__domainStartSize = new int [this->getIndexesBufferSize()];

		// Set the three pointers to the previous array
    domain = &__domainStartSize[0];
    start  = &__domainStartSize[1*numberOfDomains];
    size   = &__domainStartSize[2*numberOfDomains];

    // Allocate the temporary positions array
    next = new int [numberOfDomains];

    }


/// @brief Destructor
///
/// Delete all allocated arrays
/// @tparam T Type of data exchanged
template <class T> Content<T>::~Content() {
  if (__domainStartSize!=nullptr) delete [] __domainStartSize;
  if (next!=nullptr) delete [] next;
  if (data!=nullptr) delete [] data;
}


/// @brief Pack a some data in the data array
/// @tparam T Type of data exchanged
/// @param [in] domainIndex Recipient domain for that data
/// @param [in] buff Data
/// @param [in] count Number of objects to pack
template <class T> void Content<T>::pack(int domainIndex, T* buff, int count) {
	// Copy the data at the position for the recipient domain (getStart(domainIndex))
	// plus the size of what has already been added for that domain (getNext(domainIndex))
  memcpy(data + getStart(domainIndex) + getNext(domainIndex), buff, count*sizeof(T));
  // Add the size of the data to the size what has already been added for that domain
  getNext(domainIndex) += count;
}


/// @brief Unpack the data for a domain
/// @tparam T Type of data exchanged
/// @param [out] buff Where the unpack data is put
/// @param [in] i Index for the domain in the domain array
template <class T> void Content<T>::unpack(T* buff, int i) {
  memcpy(buff, &data[start[i]], size[i]*sizeof(T));
}


/// @brief Access to number of domains on the recipient node
template <class T> int Content<T>::getNumberOfDomains() {
  return numberOfDomains;
}


/// @brief Get start position of the data for a domain
/// @tparam T Type of data exchanged
/// @param [in] domainIndex Domain index
/// @return Position where the data for the domain start
template <class T> int& Content<T>::getStart(int domainIndex) {
  return start[findIndex(domainIndex)];
}


/// @brief Get size of the data for a domain
/// @tparam T Type of data exchanged
/// @param [in] domainIndex Domain index
/// @return Number of object for the domain
template <class T> int& Content<T>::getSize(int domainIndex) {
  return size[findIndex(domainIndex)];
}


/// @brief Get the index of a domain
/// @tparam T Type of data exchanged
/// @param [in] i Index for the domain in the domain array
/// @return Domain index
template <class T> int& Content<T>::getDomainI(int i) {
  return domain[i];
}


/// @brief Get start position of the data for a domain
/// @tparam T Type of data exchanged
/// @param [in] i Index for the domain in the domain array
/// @return Position where the data for the domain start
template <class T> int& Content<T>::getStartI(int i) {
  return start[i];
}


/// @brief Get size of the data for a domain
/// @tparam T Type of data exchanged
/// @param [in] i Index for the domain in the domain array
/// @return Number of object for the domain
template <class T> int& Content<T>::getSizeI(int i) {
  return size[i];
}


/// @brief Get the position where new data must be added for a domain
/// @tparam T Type of data exchanged
/// @param [in] domainIndex Domain index
/// @return Size of the data already added to the content for the domain
template <class T> int& Content<T>::getNext(int domainIndex) {
  return next[findIndex(domainIndex)];
}


/// @brief Get an object from the data (not used)
/// @tparam T Type of data exchanged
/// @param [in] i Position of the object to get
/// @return Reference the the object
template <class T> T& Content<T>::getData(int i) {
  return data[i];
}


/// @brief Set the start positions knowing the sizes
/// @tparam T Type of data exchanged
template <class T> void Content<T>::setStart() {
  start[0] = 0;
  for (int i=1; i<numberOfDomains; ++i) start[i] = start[i-1] + size[i-1];
}


/// @brief Recover the total size of the data and allocate the data array
/// @tparam T Type of data exchanged
template <class T> void Content<T>::setDataBuffer() {
  totalSize = 0;
  for (int i=0; i<numberOfDomains; ++i) totalSize += size[i];
  data = new T [totalSize];
}


/// @brief Get the index for a domain in the domain array
/// @tparam T Type of data exchanged
/// @param [in] domainIndex Domain index
/// @return Index for the domain in the domain array or -1 if the domain was not found
template <class T> int Content<T>::findIndex(int domainIndex) {

  for (int i=0; i<numberOfDomains; ++i) {
    if (domain[i]==domainIndex) return i;
  }

  return -1;

}


/// @brief Get the size of the big domains/starts/sizes array
/// @tparam T Type of data exchanged
/// @return Size
template <class T> int Content<T>::getIndexesBufferSize() {
  return 3*numberOfDomains;
}


/// @brief Acessor to the big domains/starts/sizes array
/// @tparam T Type of data exchanged
template <class T> int* Content<T>::getIndexesBuffer() {
  return __domainStartSize;
}


/// @brief Acessor the the total size of the data
/// @tparam T Type of data exchanged
template <class T> int Content<T>::getDataBufferSize() {
  return totalSize;
}


/// @brief Accessor to the data
/// @tparam T Type of data exchanged
 template <class T>
T* Content<T>::getDataBuffer() {
  return data;
}

#endif // __CONTENT_HPP_INCLUDED
