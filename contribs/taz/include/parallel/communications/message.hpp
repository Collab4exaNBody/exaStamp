/// @file 
/// @brief Declaration of the message classes

#ifndef __MESSAGE_HPP_INCLUDED
#define __MESSAGE_HPP_INCLUDED


#include <map>
#include <vector>

#include "utils/array/array.hpp"


/// @brief Class to store what data is send to which node
/// @tparam Type of data in the message
template <class T> class MessageSend {

	/// @brief Shortcut for a map that index vectors of data by integer indexes
  typedef std::map<int, std::vector<T> > map_t;
  /// @brief Shortcut for an iterator on map_t
  typedef typename map_t::iterator iterator;
  /// @brief Shortcut for a constant iterator on map_t
  typedef typename map_t::const_iterator const_iterator;

public:

  MessageSend(int index, const Array<uint>& neighbors);
  ~MessageSend();

  template <class U>  MessageSend(const MessageSend<U>& msg);

  int getDomainIndex() const ;
  
  uint getTotalSize();

  uint getSize(int domain);

  T* getData(int domain);

  void push(int domain, const T& object);
  void push(int domain, const T* object, const int n);

  void clear();

  iterator begin();
  iterator end();

  const_iterator cbegin() const;
  const_iterator cend() const;

public:

  int domainIndex; ///< Index of the current domain
  map_t data; ///< Map of the data to send indexed by the domain to send to

};


/// @brief Class to store the received data
/// @tparam Type of data in the message
template <class T> class MessageRecv {

public:

  MessageRecv(int index);
  ~MessageRecv();

  template <class U> MessageRecv(const MessageRecv<U>& msg);

  int getDomainIndex() const;

  uint getSize();

  T* getData();

  void setSize(uint size);

  void clear();

protected:

  int domainIndex; ///< Index of the current domain
  std::vector<T> data; ///< Data received

};


#include "parallel/communications/message.hxx"

#endif // __MESSAGE_HPP_INCLUDED
