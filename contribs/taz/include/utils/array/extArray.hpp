/// @file 
/// @brief Class ExtArray

#ifndef __EXT_ARRAY_HPP_INCLUDED
#define __EXT_ARRAY_HPP_INCLUDED


#include "utils/array/array.hpp"
#include "utils/auxMath.hpp"


// Check bounds in access operator (only in debug mode)
#ifdef  __Array_check_bounds
#define __ExtArray_check_bounds
#endif


/// @brief An array class with slightly more complicated assignment methods
///
/// Size constrained and vectorization friendly
/// @tparam T Encased class
/// @tparam CHUNK Imposed divisor of the size, default 8
/// @tparam align Alignment, default 0
template <class T, uint CHUNK=8, size_t align=0> class ExtArray : public Array<T, align> {

public:

  ExtArray();
  ExtArray(const uint n);
  ExtArray(const uint n, const T& value);
  ExtArray(const ExtArray<T>& array);

  virtual ~ExtArray();

  /// @brief Get the size
  /// @return Number of elements
  virtual inline const uint& size() const { return m_numberOfElements; }

  /// @brief Get the capacity
  /// @return size of the memory allocation
  inline uint capacity() const { return this->m_size; }

  ExtArray<T>& operator = (const ExtArray<T>& array);

  virtual typename Array<T, align>::iterator begin();
  virtual typename Array<T, align>::iterator end();

  /// @brief Accessor to the last element
  /// @tparam T Encased class
  /// @return Reference of the last element
  inline T& back() { return this->m_block[m_numberOfElements-1]; }

  /// @brief Constant accessor to the last element
  /// @tparam T Encased class
  /// @return Constant reference of the last element
  inline const T& back() const { return this->m_block[m_numberOfElements-1]; }

  void push_back(const T& object);
  void pop_back(uint size=1);

  void pop(uint index, uint size=1);

  void assign(const uint n, const T& value);

  /// @brief Remove all elements
  /// @warning This method won't free the allocated memory
  inline void clear() { m_numberOfElements = 0; }

  void reserve(const uint n);
  void resize(const uint n);
  void resize(const uint n, const T& val);

  void shrink_to_fit();

  template <class U> void map(const U& indexes);

#ifdef __ExtArray_check_bounds

  virtual       T& operator [] (uint index);
  virtual const T& operator [] (uint index) const;

#endif // __ExtArray_check_bounds

private:

  uint m_numberOfElements;  ///< Actual number of elements in the ExtArray

};


/// @brief Default constructor
///
///
template <class T, uint CHUNK, size_t align> 
inline ExtArray<T, CHUNK, align>::ExtArray()
  : Array<T, align>(), 
    m_numberOfElements(0) {}


/// @brief Constructor from size
/// @tparam CHUNK Imposed divisor of the size
/// @param [in] n Size
template <class T, uint CHUNK, size_t align> 
inline ExtArray<T, CHUNK, align>::ExtArray(const uint n)
  : Array<T, align>(roundToChunk<CHUNK>(n)), 
    m_numberOfElements(n) {}


/// @brief Constructor from size and value
/// @tparam T Encased class
/// @tparam CHUNK Imposed divisor of the size
/// @param [in] n Size
/// @param [in] value Value for all the elements
template <class T, uint CHUNK, size_t align> 
inline ExtArray<T, CHUNK, align>::ExtArray(const uint n, const T& value)
  : Array<T, align>(roundToChunk<CHUNK>(n)),
    m_numberOfElements(n) {

  for (auto& elem : *this) elem = value; 

}


/// @brief Copy constructor
/// @warning : No check on the size
/// @param [in] array Another ExtArray
template <class T, uint CHUNK, size_t align> 
inline ExtArray<T, CHUNK, align>::ExtArray(const ExtArray<T>& array)
  : Array<T, align>(array),
    m_numberOfElements(array.m_numberOfElements) {
}


/// @brief Destructor
///
///
template <class T, uint CHUNK, size_t align> 
inline ExtArray<T, CHUNK, align>::~ExtArray() {
  m_numberOfElements = 0;
}


/// @brief Assignment operator
/// @tparam T Encased class
/// @param [in] array ExtArray to copy
template <class T, uint CHUNK, size_t align> 
ExtArray<T>& ExtArray<T, CHUNK, align>::operator = (const ExtArray<T>& array) {

	// Copy the capacity
  this->m_size = array.m_size;
  // Free the pre-existent array
  this->allocator.free(this->m_block);
  this->m_block = nullptr;
  // Copy the size
  m_numberOfElements = array.m_numberOfElements;
  
  if (this->m_size>0) {
  	// Allocate m_size elements
    this->m_block = this->allocator.alloc(this->m_size);
    // Copy m_numberOfElements elements of the other array into this one
    auxMemCpy(this->m_block, array.m_block, m_numberOfElements);
  }

  return *this;

}


/// @brief Get begining of the ExtArray
/// @return Iterator to the first element
template <class T, uint CHUNK, size_t align>
inline typename Array<T, align>::iterator ExtArray<T, CHUNK, align>::begin() {
  return typename Array<T, align>::iterator(this, 0);
}


/// @brief Get end of the ExtArray
/// @return Iterator to the first element out of range
template <class T, uint CHUNK, size_t align>
inline typename Array<T, align>::iterator ExtArray<T, CHUNK, align>::end() {
  return typename Array<T, align>::iterator(this, m_numberOfElements);
}


/// @brief Add an element at the end
/// @tparam T Encased class
/// @param [in] object Element to add
template <class T, uint CHUNK, size_t align> 
inline void ExtArray<T, CHUNK, align>::push_back(const T& object) {

	// If capacity is enough to add one element
  if (m_numberOfElements < this->m_size) {
  	// Do it
    this->m_block[m_numberOfElements++] = object;
  }
  // Else
  else {
    T tmp = object; // Why this copy ?
    // Reserve more memory
    this->reserve(this->m_size + 1);
    // Add the element
    this->m_block[m_numberOfElements++] = tmp;

  }

}


/// @brief Remove last element
///
///
template <class T, uint CHUNK, size_t align> 
inline void ExtArray<T, CHUNK, align>::pop_back(uint size) {

#ifdef __ExtArray_check_bounds
  if (m_numberOfElements<size) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'ExtArray::pop_back()' : array not that full ! Size = "
	     << size
	     << std::endl;
  }
#endif // __ExtArray_check_bounds

  m_numberOfElements -= size;
}


/// @brief Remove specified elements
/// @warning Does not retain the order of the elements, use with caution
/// @param [in] index Index of the removed element
/// @param [in] size Number of elements to remove, default 1
template <class T, uint CHUNK, size_t align> 
inline void ExtArray<T, CHUNK, align>::pop(uint index, uint size) {

#ifdef __ExtArray_check_bounds
  if (m_numberOfElements<size) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'ExtArray::pop(uint)' : array not that full ! Size = "
	     << size
	     << std::endl;
  }
  if (index>=m_numberOfElements) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'ExtArray::pop(uint)' : index out of bounds (" 
	     << index << "/" << m_numberOfElements << ") !" 
	     << std::endl;
  }
#endif // __ExtArray_check_bounds
  m_numberOfElements-=size;
  for(uint i(0); i<size; ++i) this->m_block[index+i] = this->m_block[m_numberOfElements+i];

}


/// @brief Assign size and value to array
/// @tparam T Encased class
/// @param [in] n Size
/// @param [in] value Value for all the elements
template <class T, uint CHUNK, size_t align> 
inline void ExtArray<T, CHUNK, align>::assign(const uint n, const T& value) {
  this->clear();
  this->resize(n, value);
}


/// @brief Reserve space
/// @tparam CHUNK Imposed divisor of the size
/// @param [in] n Desired capacity
template <class T, uint CHUNK, size_t align> 
void ExtArray<T, CHUNK, align>::reserve(const uint n) {

	// Get new capacity
  const uint nn = roundToChunk<CHUNK>(n);

  // If more than actual
  if (nn > this->m_size) {

  	// Allocate an array for the new capacity
    T* tmp = this->allocator.alloc(nn);

    // Copy the elements in the new array
    auxMemCpy(tmp, this->m_block, this->m_numberOfElements);

    // Free the old one
    this->allocator.free(this->m_block);

    // Reassign array and capacity
    this->m_block = tmp;
    this->m_size = nn;

  }

}


/// @brief Resize array
/// @param [in] n New size
template <class T, uint CHUNK, size_t align> 
void ExtArray<T, CHUNK, align>::resize(const uint n) {

	// If size is decreasing
  if (n < m_numberOfElements) {

  	// Just resize
    m_numberOfElements = n;

  }
  // If size is increasing
  else if (n>m_numberOfElements) {

  	// Reserve more space
    this->reserve(n);
    // And resize
    m_numberOfElements = n;

  }

}


/// @brief Resize with a value
/// @tparam T Encased class
/// @param [in] n New size
/// @param [in] val New value for all the new elements
template <class T, uint CHUNK, size_t align> 
void ExtArray<T, CHUNK, align>::resize(const uint n, const T& val) {

	// If size is decreasing
  if (n < m_numberOfElements) {

  	// Just resize
    m_numberOfElements = n;

  }
  // If size is increasing
  else if (n>m_numberOfElements) {

  	// Reserve more space
    this->reserve(n);
    // Assign value to new elements
    for (uint i=m_numberOfElements; i<n; ++i)
      this->m_block[i] = val;
    // And resize
    m_numberOfElements = n;

  }

}


/// @brief Shrink array so that the total size is exactly the number of elements (rouded to the chunk)
/// @tparam CHUNK Imposed divisor of the size
template <class T, uint CHUNK, size_t align> 
void ExtArray<T, CHUNK, align>::shrink_to_fit() {

	// Get the new size
  const uint n = roundToChunk<CHUNK>(m_numberOfElements);

  // If size is decreasing
  if (this->m_size > n) {

  	// Allocate an array for the new size
    T* tmp = this->allocator.alloc(n);
    // Copy the elements in the new array
    auxMemCpy(tmp, this->m_block, this->m_size);

    // Reassign capacity
    this->m_size = n;
    // Free the old one
    this->allocator.free(this->m_block);
    // Reassign array
    this->m_block = tmp;

  }

}


/// @brief Rearrange elements according to indexes given in an array
/// @tparam T Encased class
/// @tparam U Class of the array containing indexes (must have a constant access operator by [], no check)
/// @param [in] indexes Indexes (must have a minimal size of m_numberOfElements, no check)
template <class T, uint CHUNK, size_t align> template <class U>
void ExtArray<T, CHUNK, align>::map(const U& indexes) {

	// Allocate a temporary array
  T* tmp = this->allocator.alloc(this->m_size);

  // Fill it with reordered elements
  for (uint i=0; i<m_numberOfElements; ++i) 
    tmp[i] = this->m_block[indexes[i]];

  // Delete original array
  this->allocator.free(this->m_block);
  // Put the new array instead
  this->m_block = tmp;

}


#ifdef __ExtArray_check_bounds


/// @brief Element access operator
/// @version Check bounds
/// @tparam T Encased class
/// @param [in] index Index of the element to access
/// @return Element
template <class T, uint CHUNK, size_t align> 
inline T& ExtArray<T, CHUNK, align>::operator [] (uint index) {

  if (index>=m_numberOfElements) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'ExtArray::operator [] (uint)' : index out of bounds (" 
	     << index << "/" << m_numberOfElements << ") !" << std::endl;
  }

  return this->m_block[index];

}


/// @brief Constant element access operator
/// @version Check bounds
/// @tparam T Encased class
/// @param [in] index Index of the element to access
/// @return Element
template <class T, uint CHUNK, size_t align> 
inline const T& ExtArray<T, CHUNK, align>::operator [] (uint index) const {

  if (index>=m_numberOfElements) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'ExtArray::operator [] (uint) const' : index out of bounds (" 
	     << index << "/" << m_numberOfElements << ") !" << std::endl;
  }

  return this->m_block[index];

}


#endif // __ExtArray_check_bounds

#endif // __EXT_ARRAY_HPP_INCLUDED
