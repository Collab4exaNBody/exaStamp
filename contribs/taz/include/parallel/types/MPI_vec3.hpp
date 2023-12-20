/// @file 
/// @brief Templated function for the creation of MPI types for vec3

#ifndef __MPI_VEC3_HPP_INCLUDED
#define __MPI_VEC3_HPP_INCLUDED


#include "parallel/mympi.hpp"
#include "parallel/types/MPI_typesUtils.hpp"

#include "utils/vec3/vec3.hpp"


/// @brief Create the MPI_Datatype for a vec3
/// @tparam T Class inside the vec3
template <class T> void create_MPI_type_vec3() {
  MPI_Datatype newtype;
  MPI_Type_contiguous(VEC3_NDIMS, MPI__Type_get<T>(), &newtype);
  MPI_Type_commit(&newtype);
  MPI__Type_add< vec3<T> >(newtype);
}

#endif // __MPI_VEC3_HPP_INCLUDED
