/// @file 
/// @brief Definition of the message center

#ifndef __MESSAGE_CENTER_HPP_INCLUDED
#define __MESSAGE_CENTER_HPP_INCLUDED


#include "parallel/communications/session.hpp"


/// @brief Shortcut for a template that depend on the type to exchange
#define TMPLMC template <class T>
/// @brief Shortcut for the templated message center
#define TMPL_MessageCenter MessageCenter<T>


/// @brief Tool that handle the point-to-point node communications for a specific data type
///
/// You must create new MessageCenter each time you want to exchange a new data type
/// @tparam T Data type exchanged in this center
TMPLMC class MessageCenter {

public:

	/// @brief Constructor
	/// @param [in] comm_ Pointer to the communication manager
	/// @param [in] type_ Unique identifier of this communication manager in the program
  MessageCenter(CommManager* comm_, ISession::Type type_)
    : comm(comm_), type(type_), 
      toSend(nullptr), toRecv(nullptr) {}

  /// @brief Destructor
	///
  /// Delete all remaining messages
  virtual ~MessageCenter() {
    if (toSend != nullptr) delete toSend;
    if (toRecv != nullptr) delete toRecv;
  }

  void init(uint domainIndex, const Array<uint>& neighbors);

  template <class U>
  void init(const MessageCenter<U>& messageCenter);

  void send();
  void collect();
  void clean();

  MessageSend<T>* getMessageSend();
  MessageRecv<T>* getMessageRecv();
  
  const MessageSend<T>& getMessageSend() const ;
  const MessageRecv<T>& getMessageRecv() const ;
  
protected:

  CommManager* comm; ///< Pojnter to the communication manager
  ISession::Type type; ///< Unique identifier of this communication manager in the program

  MessageSend<T>* toSend; ///< Message to send
  MessageRecv<T>* toRecv; ///< Message received

};


/// @brief Initialize the messages
/// @param [in] domainIndex index of the domain
/// @param [in] neighbors Neighbors of the domain
TMPLMC inline void TMPL_MessageCenter::init(uint domainIndex, const Array<uint>& neighbors) {
  toSend = new MessageSend<T>(domainIndex, neighbors);
  toRecv = new MessageRecv<T>(domainIndex);
}


/// @brief Initialize the messages by copying those of another message center
/// @tparam U Data type exchanged in the other message center
/// @param [in] messageCenter Other message center
TMPLMC template <class U> inline void TMPL_MessageCenter::init(const MessageCenter<U>& messageCenter) {
  toSend = new MessageSend<T>(messageCenter.getMessageSend());
  toRecv = new MessageRecv<T>(messageCenter.getMessageRecv());
}


/// @brief Send the message to send
TMPLMC inline void TMPL_MessageCenter::send() {
	// Get the session corresponding to the identifier of this message center
  Session<T>* tmp = comm->getSession<T>(type);
  // Send the message to send on the session
  tmp->pushMessages(toSend);
}


/// @brief Collect data into the message received
TMPLMC inline void TMPL_MessageCenter::collect() {
	// Get the session corresponding to the identifier of this message center
  Session<T>* tmp = comm->getSession<T>(type);
  // Collect the data
  tmp->collectMessages(toRecv);
}


/// @brief End the session and clear the message received
TMPLMC inline void TMPL_MessageCenter::clean() {
  comm->closeSession(type);
  toRecv->clear();
}


/// @brief Accessor to the message to send
/// @return Pointer to the message to send
TMPLMC inline MessageSend<T>* TMPL_MessageCenter::getMessageSend() {
  return toSend;
}


/// @brief Accessor to the message received
/// @return Pointer to the message received
TMPLMC inline MessageRecv<T>* TMPL_MessageCenter::getMessageRecv() {
  return toRecv;
}


/// @brief Constant accessor to the message to send
/// @return Message to send
TMPLMC inline const MessageSend<T>& TMPL_MessageCenter::getMessageSend() const {
  return *toSend;
}


/// @brief Constant accessor to the message received
/// @return Message received
TMPLMC inline const MessageRecv<T>& TMPL_MessageCenter::getMessageRecv() const {
  return *toRecv;
}


#undef TMPLMC
#undef TMPL_MessageCenter

#endif // __MESSAGE_CENTER_HPP_INCLUDED
