/// @file 
/// @brief Definition of an iterator class for Array

#include <iterator>


/// @brief Iterator on Array<T, align>
/// @tparam T Encased class
/// @tparam align Alignment
template <class T, size_t align>
class Array<T, align>::iterator : public std::iterator<std::random_access_iterator_tag, T> {

public:

  iterator();
  ~iterator();

  iterator(Array<T, align>* ref, const uint idx);
  iterator(const iterator& it);

  iterator& operator = (const iterator& it);

  bool operator == (const iterator& it) const;
  bool operator != (const iterator& it) const;

  T& operator * () const;
  T* operator -> () const;

  iterator& operator += (const int n);
  iterator& operator -= (const int n);

  iterator& operator ++ ();
  iterator& operator -- ();

  iterator operator + (const int n) const;
  iterator operator - (const int n) const;

  iterator operator ++ (int);
  iterator operator -- (int);

  iterator operator - (const iterator& it) const;

  bool operator < (const iterator& it) const;
  bool operator > (const iterator& it) const;

  bool operator <= (const iterator& it) const;
  bool operator >= (const iterator& it) const;

  /// @brief Used only for debug
  /// @return Iterator position
  uint pos() const { return m_idx; }

private:

  Array<T, align>* m_ref; ///< Pointer to the Array
  uint m_idx; ///< Position of the iterator

};


/// @brief Default constructor
template <class T, size_t align>
inline Array<T, align>::iterator::iterator() : m_ref(nullptr), m_idx(0) {
}


/// @brief Destructor (nothing to do)
template <class T, size_t align>
inline Array<T, align>::iterator::~iterator() {
}


/// @brief Constructor from an Array and a position
/// @param [in,out] ref Pointer to the Array
/// @param [in] idx Start position
template <class T, size_t align>
inline Array<T, align>::iterator::iterator(Array<T, align>* ref, const uint idx) 
  : m_ref(ref), m_idx(idx) {
}


/// @brief Copy constructor
/// @param [in] it Another iterator
template <class T, size_t align>
inline Array<T, align>::iterator::iterator(const iterator& it)
    : m_ref(it.m_ref), m_idx(it.m_idx) {
}


/// @brief Assignment operator
/// @param [in] it Iterator to assign
template <class T, size_t align>
inline typename Array<T, align>::iterator& Array<T, align>::iterator::operator = (const iterator& it) {
  m_ref = it.m_ref;
  m_idx = it.m_idx;
  return *this;
}


/// @brief Equality operator
/// @param [in] it Iterator to test
template <class T, size_t align>
inline bool Array<T, align>::iterator::operator == (const iterator& it) const {
  return m_ref==it.m_ref && m_idx==it.m_idx;
}


/// @brief Inequality operator
/// @param [in] it Iterator to test
template <class T, size_t align>
inline bool Array<T, align>::iterator::operator != (const iterator& it) const {
  return !(*this==it);
}


/// @brief Indirection operator
///
/// Access to content
/// @tparam T Encased class
/// @return Content pointed by the iterator
template <class T, size_t align>
inline T& Array<T, align>::iterator::operator * () const {
  return m_ref->m_block[m_idx];
}

/// @brief Dereference operator
///
/// Access to address
/// @tparam T Encased class
/// @return Address pointed by the iterator
template <class T, size_t align>
inline T* Array<T, align>::iterator::operator -> () const {
  return m_ref->m_block + m_idx;
}



/// @brief Addition assignment operator
///
/// Shift iterator n times to the right
/// @param [in] n Number of shifts
template <class T, size_t align>
inline typename Array<T, align>::iterator& Array<T, align>::iterator::operator += (const int n) {
  m_idx += n;
  return *this;
}


/// @brief Substraction assignment operator
///
/// Shift iterator n times to the left
/// @param [in] n Number of shifts
template <class T, size_t align>
inline typename Array<T, align>::iterator& Array<T, align>::iterator::operator -= (const int n) {
  m_idx -= n;
  return *this;
}


/// @brief Pre-increment operator
///
/// Right shift 1 element
template <class T, size_t align>
inline typename Array<T, align>::iterator& Array<T, align>::iterator::operator ++ () {
  return *this+=1;
}


/// @brief Pre-decrement operator
///
/// Left shift 1 element
template <class T, size_t align>
inline typename Array<T, align>::iterator& Array<T, align>::iterator::operator -- () {
  return *this-=1;
}


/// @brief Addition operator (integer)
/// @param [in] n Integer to add
/// @return Right shifted iterator
template <class T, size_t align>
inline typename Array<T, align>::iterator Array<T, align>::iterator::operator + (const int n) const {
  return iterator(*this)+=n;
}


/// @brief Subtraction operator (integer)
/// @param [in] n Integer to subtract
/// @return Left shifted iterator
template <class T, size_t align>
inline typename Array<T, align>::iterator Array<T, align>::iterator::operator - (const int n) const {
  return iterator(*this)-=n;
}


/// @brief Post-increment operator
///
/// Right shift 1 element
template <class T, size_t align>
inline typename Array<T, align>::iterator Array<T, align>::iterator::operator ++ (int) {
  iterator clone (*this);
  ++*this;
  return clone;
}


/// @brief Post-decrement operator
///
/// Left shift 1 element
template <class T, size_t align>
inline typename Array<T, align>::iterator Array<T, align>::iterator::operator -- (int) {
  iterator clone (*this);
  --*this;
  return clone;
}


/// @brief Subtraction operator (iterator)
/// @param [in] it Iterator to subtract
/// @return Copy of the first operator left shifted of the position of the second operator
template <class T, size_t align>
inline typename Array<T, align>::iterator Array<T, align>::iterator::operator - (const iterator& it) const {
  return iterator(*this)-it.m_idx;
}


/// @brief Comparison operator (<)
///
/// Compare positions
/// @param [in] it Iterator to compare
template <class T, size_t align>
inline bool Array<T, align>::iterator::operator < (const iterator& it) const {
  return m_idx < it.m_idx;
}


/// @brief Comparison operator (>)
///
/// Compare positions
/// @param [in] it Iterator to compare
template <class T, size_t align>
inline bool Array<T, align>::iterator::operator > (const iterator& it) const {
  return m_idx > it.m_idx;
}


/// @brief Comparison operator (<=)
///
/// Compare positions
/// @param [in] it Iterator to compare
template <class T, size_t align>
inline bool Array<T, align>::iterator::operator <= (const iterator& it) const {
  return m_idx <= it.m_idx;
}


/// @brief Comparison operator (>=)
///
/// Compare positions
/// @param [in] it Iterator to compare
template <class T, size_t align>
inline bool Array<T, align>::iterator::operator >= (const iterator& it) const {
  return m_idx >= it.m_idx;
}


/// @brief Get begining of the Array
/// @return Iterator to the first element
template <class T, size_t align>
inline typename Array<T, align>::iterator Array<T, align>::begin() {
  return iterator(this, 0);
}


/// @brief Get end of the Array
/// @return Iterator to the first element out of range
template <class T, size_t align>
inline typename Array<T, align>::iterator Array<T, align>::end() {
  return iterator(this, m_size);
}
