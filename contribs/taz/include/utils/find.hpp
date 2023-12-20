/// @file 
/// @brief A function to find a value in an array

#ifndef __FIND_HPP_INCLUDED
#define __FIND_HPP_INCLUDED


/// @brief Search a value in an array
/// @tparam T Type of the value
/// @param [in] val Value to find
/// @param [in] ptr Array to search
/// @param [in] size Size of the array
/// @param [out] found Indicates if value was found
/// @param [out] index Indicates where value was found
/// @warning ptr must be sorted (<) for this to work !
template <class T>
inline void find(const T& val, const T* ptr, const size_t size, bool& found, size_t& index) {

  found = false;

  size_t idxMin  = 0;
  size_t idxMax  = size;
  size_t idxSize = size;

  if (val<ptr[0] || val>ptr[size-1]) 
    idxSize = 0;

  while (idxSize>0) {

    index = idxMin + idxSize/2;
    
    if (ptr[index]<val) {
      idxMin = index+1;
    }
    else if (ptr[index]>val) {
      idxMax = index;
    }
    else {
      found=true;
      break;
    };

    idxSize = idxMax-idxMin;

  }

}

#endif // __FIND_HPP_INCLUDED
