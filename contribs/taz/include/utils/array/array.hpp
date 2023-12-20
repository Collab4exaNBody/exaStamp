/// @file 
/// @brief Definition of the class Array

#ifndef __ARRAY_HPP_INCLUDED
#define __ARRAY_HPP_INCLUDED


#include <immintrin.h>
#include <iostream>

/*
#include <malloc.h>

extern int posix_memalign (void **, size_t, size_t);

static __inline void *
_mm_malloc (size_t size, size_t alignment)
{
  void *ptr;
  if (alignment == 1)
    return malloc (size);
  if (alignment == 2 || (sizeof (void *) == 8 && alignment == 4))
    alignment = sizeof (void *);
  if (posix_memalign (&ptr, alignment, size) == 0)
    return ptr;
  else
    return NULL;
}

static __inline void
_mm_free (void * ptr)
{
  free (ptr);
}
*/
#include "utils/auxMem.hpp"


/// @brief Wrapper class for new and delete of an array in case aligned memory is wanted
///
/// Not the cleanest way to do it ... 
/// @tparam T Encased class
/// @tparam align Alignment
template <class T, size_t align> class __allocator {

public: 

  /// @brief Default constructor
  __allocator() {}
  /// @brief Destructor (nothing to do)
  ~__allocator() {}
  
  /// @brief Allocate an aligned array
  /// @tparam T Encased class
  /// @tparam align Alignment
  /// @param [in] n Size of the array
  /// @return Address to the allocated array
  inline T* alloc(uint n) const {
    return (T*) _mm_malloc(n*sizeof(T), align);
  }

  /// @brief Free the array
  /// @tparam T Encased class
  /// @param [in,out] ptr Pointer to the array
  inline void free(T* ptr) const {
    if (ptr!=nullptr) _mm_free(ptr);
  }

};


/// @brief Wrapper class for new and delete of an array in case aligned memory is not wanted
/// @tparam T Encased class
template <class T> class __allocator<T, 0> {

public: 

  /// @brief Default constructor
  __allocator() {}
  /// @brief Destructor (nothing to do)
  ~__allocator() {}
  
  /// @brief Create an array
  /// @tparam T Encased class
  /// @param [in] n Size of the array
  /// @return Address to the allocated array
  inline T* alloc(uint n) const {
    return new T [n];
  }

  /// @brief Free the array
  /// @tparam T Encased class
  /// @param [in,out] ptr Pointer to the array
  inline void free(T* ptr) const {
    if (ptr!=nullptr) {
      delete [] ptr;
      ptr=nullptr;
    }
    
  }

};


/// @brief A simple one-dimensional array
/// @warning To avoid memory leaks, do not use this class with a template 
/// parameter which has a @c new statement in its default constructor
/// @tparam T Encased class
/// @tparam align Alignment, default 0
template <class T, size_t align=0> class Array {

public:

  Array();
  Array(const uint n);
  Array(const uint n, const T& value);
  Array(const Array<T, align>& array);

  virtual ~Array();

  Array<T, align>& operator = (const Array<T, align>& array);

  /// @brief Accessor to the size
  inline const uint& size() const { return m_size; }

  /// @brief Allocate an Array of specified size (delete pre-existent Array)
  /// @param [in] size Size of the Array
  inline void alloc(const uint size){
    m_size=size; 
    if (m_block!=nullptr) allocator.free(m_block);
    if (m_size>0) m_block = allocator.alloc(m_size);
  }

#ifndef __Array_check_bounds

  /// @brief Element access
  /// @version Default
  /// @tparam T Encased class
  /// @param [in] index Index of the element to access
  /// @return Element
  inline T& operator [] (const uint index) { return m_block[index]; }

  /// @brief Element access (const)
  /// @version Default
  /// @tparam T Encased class
  /// @param [in] index Index of the element to access
  /// @return Element
  inline const T& operator [] (const uint index) const { return m_block[index]; }

#else

  inline T& operator [] (const uint index);
  inline const T& operator [] (const uint index) const;

#endif // __Array_check_bounds

  /// @brief Get a pointer to the data
  /// @tparam T Encased class
  /// @return Pointer to the data
  inline T* data() { return m_block; }
  /// @brief Get a constant pointer to the data
  /// @tparam T Encased class
  /// @return Constant pointer to the data
  inline const T* data() const { return m_block; }

  class iterator;
  iterator begin();
  iterator end();

  void set(const T& val);

  template <class U> void map(const U& indexes);

protected:

  static const __allocator<T, align> allocator; ///< Tool to build and destruct the Array

  uint m_size; ///< Number of elements in the Array
  T* m_block; ///< Pointer to first element

};


template <class T, size_t align> const __allocator<T, align> Array<T, align>::allocator = __allocator<T, align>();


/// @brief Default constructor
///
///
template <class T, size_t align> 
inline Array<T, align>::Array() : m_size(0), m_block(nullptr) {}


/// @brief Constructor from a size
/// @param [in] n Size of the Array
template <class T, size_t align> 
inline Array<T, align>::Array(const uint n) : m_size(n), m_block(nullptr) {
  // Allocate m_size elements
  if (m_size>0) m_block = allocator.alloc(m_size);
}


/// @brief Constructor from a size and value
/// @tparam T Encased class
/// @param [in] n Size of the Array
/// @param [in] value Value for all the elements
template <class T, size_t align> 
inline Array<T, align>::Array(const uint n, const T& value) : m_size(n), m_block(nullptr) {
  if (m_size>0) {
    // Allocate m_size elements
    m_block = allocator.alloc(m_size);
    // Set all elements to value
    this->set(value);
  }
}


/// @brief Copy constructor
/// @param [in] array Another Array
template <class T, size_t align> 
inline Array<T, align>::Array(const Array<T, align>& array) : m_size(array.m_size), m_block(nullptr) {
  if (m_size>0) {
    // Allocate m_size elements
    m_block = allocator.alloc(m_size);
    // Copy the other array into this one
    auxMemCpy(m_block, array.m_block, m_size);
  }
}


/// @brief Destructor
///
///
template <class T, size_t align> 
inline Array<T, align>::~Array() {
  // Set size to 0
  m_size=0;
  // Free the elements
  allocator.free(m_block);
}


/// @brief Assignment operator
/// @param [in] array Array to assign
template <class T, size_t align> 
Array<T, align>& Array<T, align>::operator = (const Array<T, align>& array) {
  // Copy the size
  m_size = array.m_size;

  // Free the pre-existent elements
  allocator.free(m_block);

  if (m_size>0) {
    // Allocate m_size elements
    m_block = allocator.alloc(m_size);
    // Copy the other array into this one
    auxMemCpy(m_block, array.m_block, m_size);
  }
  else {
    m_block = nullptr;
  }

  return *this;
}


#ifdef __Array_check_bounds


/// @brief Element access
/// @version Check bounds
/// @tparam T Encased class
/// @param [in] index Index of the element to access
/// @return Element
template <class T, size_t align> 
inline T& Array<T, align>::operator [] (const uint index) {

  if (index>=m_size) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'Array::operator [] (uint)' : index out of bounds (" 
	     << index << "/" << m_size << ") !" << std::endl;
  }

  return m_block[index];

}


/// @brief Const element access
/// @version Check bounds
/// @tparam T Encased class
/// @param [in] index Index of the element to access
/// @return Element
template <class T, size_t align> 
inline const T& Array<T, align>::operator [] (const uint index) const {

  if (index>=m_size) {
    std::cerr<< "[ERROR] " << __FILE__ << ":" << __LINE__ << ": " 
	     << "in function 'Array::operator [] (uint) const' : index out of bounds (" 
	     << index << "/" << m_size << ") !" << std::endl;
  }

  return m_block[index];

}


#endif // __Array_check_bounds


/// @brief Set all the elements to a value
/// @tparam T Encased class
/// @param [in] val Value for all the elements
template <class T, size_t align>
inline void Array<T, align>::set(const T& val) {
  for (uint i=0; i<m_size; ++i) {
    m_block[i] = val;
  }
}


/// @brief Rearrange elements according to indexes given in an array
/// @tparam T Encased class
/// @tparam U Class of the array containing indexes (must have a constant access operator by [], no check)
/// @param [in] indexes Indexes (must have a minimal size of m_size, no check)
template <class T, size_t align> template <class U>
void Array<T, align>::map(const U& indexes) {

  // Allocate a temporary array
  T* tmp = allocator.alloc(m_size);

  // Fill it with reordered elements
  for (uint i=0; i<m_size; ++i) 
    tmp[i] = m_block[indexes[i]];

  // Delete original array
  allocator.free(m_block);
  // Put the new array instead
  m_block = tmp;

}


// CONTAINS : iterator class and related functions
#include "utils/array/array_iterator.hxx"

#endif // __ARRAY_HPP_INCLUDED
