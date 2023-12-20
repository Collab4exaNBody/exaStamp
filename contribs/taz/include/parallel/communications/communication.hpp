/// @file 
/// @brief Definition and implementation of (template) class Communication

#ifndef __COMMUNICATION_HPP_INCLUDED
#define __COMMUNICATION_HPP_INCLUDED


#include <set>

#include "parallel/mympi.hpp"
#include "parallel/communications/content.hpp"


/// @brief Class that handle the data exchange between this node and another node
/// @tparam T Type of data exchanged
template <class T> class Communication {

public:

  Communication();
  ~Communication();

  Communication(const MPI_Comm& comm_, int rank_, const std::set<int>& domains, int totalSize);

  Content<T>* getSendContent();
  Content<T>* getRecvContent();

  void sendRecvSizes(int destNode, MPI_Request* req, int nbLocalDomains);

  void allocRecvBuffers(int destNode);
  void exchangeMessages(int destNode, MPI_Request* req);


protected:

  const MPI_Comm& comm;  ///< MPI communicator used for communications
  int rank; ///< Rank of the node

  Content<T>* toSend;  ///< Content to send
  Content<T>* toRecv;  ///< Content to receive

  int sizeToRecv;  ///< Total size of the content to receive (has to be allocated before being received)

};


/// @brief Default constructor
/// @tparam T Type of data exchanged
template <class T> Communication<T>::Communication()
  : comm(MPI_COMM_WORLD), rank(MPI_COMM_NULL), 
    toSend(nullptr), toRecv(nullptr), sizeToRecv() {
}



/// @brief Constructor with arguments
/// @tparam T Type of data exchanged
/// @param [in] comm_ MPI communicator
/// @param [in] rank_ Node rank
/// @param [in] domains Set of the domains for which there will be communications on the other node
/// @param [in] totalSize Total size of the data to send to the other node
template <class T> Communication<T>::Communication(const MPI_Comm& comm_, int rank_, const std::set<int>& domains, int totalSize)
  : comm(comm_), rank(rank_), toSend(nullptr), toRecv(nullptr), sizeToRecv(0) {

	// Create a content for this communication
  toSend = new Content<T>(comm, domains, totalSize);

}


/// @brief Destructor
/// @tparam T Type of data exchanged
template <class T> Communication<T>:: ~Communication() {
	// Delete leaving and incoming contents
  if (toSend!=nullptr) delete toSend;
  if (toRecv!=nullptr) delete toRecv;
}


/// @brief Send and get the sizes of what will be exchanged in order to allocate received buffers
/// @tparam T Type of data exchanged
/// @param [in] destNode Index of the other node
/// @param [in,out] req Pointer to the MPI requests for this communication
/// @param [in] nbLocalDomains Number of domains on the current node
template <class T> void Communication<T>::sendRecvSizes(int destNode, MPI_Request* req, int nbLocalDomains) {

  int tag = 123; // Tag for this communications

  // If the recipient node is not the sending node
  if (destNode!=rank) {

  	// Old version : create the received content buffer from the send content buffer (probably won't work)
    // toRecv = new Content<T>(toSend);
  	// New version : create the received content buffer from the number of domains on the current node
  	toRecv = new Content<T>(comm,nbLocalDomains);

    // Start the communication of the sizes for the received content
    MPI_Irecv(toRecv->getIndexesBuffer(), toRecv->getIndexesBufferSize(), 
	      MPI__Type_get<int>(), destNode, tag, comm, req);
    // Start the communication of the sizes for the send content
    MPI_Isend(toSend->getIndexesBuffer(), toSend->getIndexesBufferSize(), 
	      MPI__Type_get<int>(), destNode, tag, comm, req+1);

  }
  // Else no MPI communication, just copy the send content into the received content
  else {
    toRecv = toSend;
  }

}


/// @brief Pointer accessor to the content to send
/// @tparam T Type of data exchanged
template <class T> Content<T>* Communication<T>::getSendContent() {
  return toSend;
}


/// @brief Pointer accessor to the content to be received
/// @tparam T Type of data exchanged
template <class T> Content<T>* Communication<T>::getRecvContent() {
  return toRecv;
}


/// @brief Allocate received content
/// @tparam T Type of data exchanged
/// @param [in] destNode Index of the other node
template <class T> void Communication<T>::allocRecvBuffers(int destNode) {
  if (destNode!=rank) toRecv->setDataBuffer();
}


/// @brief Exchange the data
/// @tparam T Type of data exchanged
/// @param [in] destNode Index of the other node
/// @param [in,out] req Pointer to the MPI requests for this communication
template <class T> void Communication<T>::exchangeMessages(int destNode, MPI_Request* req) {

  int tag = 456; // Tag for this communications

  // If the recipient node is not the sending node
  if (destNode!=rank) {

  	// Start the communication to get the data
    MPI_Irecv(toRecv->getDataBuffer(), toRecv->getDataBufferSize(), 
    	      MPI__Type_get<T>(), destNode, tag, comm, req);
    // Start the communication to give the data
    MPI_Isend(toSend->getDataBuffer(), toSend->getDataBufferSize(), 
    	      MPI__Type_get<T>(), destNode, tag, comm, req+1);

  }

}

#endif // __COMMUNICATION_HPP_INCLUDED
