/// @file 
/// @brief Templated function for the creation of MPI types for mat3

#ifndef __MPI_MAT3_HPP_INCLUDED
#define __MPI_MAT3_HPP_INCLUDED


#include "parallel/mympi.hpp"
#include "parallel/types/MPI_typesUtils.hpp"

#include "utils/mat3/mat3.hpp"


/// @brief Create the MPI_Datatype for a mat3
/// @tparam T Class inside the mat3
template <class T> void create_MPI_type_mat3() {
  MPI_Datatype newtype;
  MPI_Type_contiguous(MAT3_NDIMS*MAT3_NDIMS, MPI__Type_get<T>(), &newtype);
  MPI_Type_commit(&newtype);
  MPI__Type_add< mat3<T> >(newtype);
}

#endif // __MPI_MAT3_HPP_INCLUDED
