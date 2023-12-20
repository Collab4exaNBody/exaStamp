/// @file 
/// @brief Header for the general functions concerning all the MPI_Datatypes

#ifndef __MPI_TYPES_UTILS_HPP_INCLUDED
#define __MPI_TYPES_UTILS_HPP_INCLUDED


#include <map>
#include <typeindex>

#include "parallel/mympi.hpp"
#include "parallel/types/MPI_ghost.hpp"
#include "parallel/types/MPI_inMol.hpp"


/// @brief Get the MPI_Datatype matching type
/// @tparam T Type to match
///
/// The MPI_Datatype must have been created before ...
template <class T> MPI_Datatype MPI__Type_get() {
  extern std::map<std::type_index, MPI_Datatype> MPI__TYPES_MAP;
  return MPI__TYPES_MAP[typeid(T)];
}


/// @brief Add a MPI_Datatype to the MPI types map
/// @tparam T Corresponding non MPI type
template <class T> void MPI__Type_add(MPI_Datatype type) {
  extern std::map<std::type_index, MPI_Datatype> MPI__TYPES_MAP;
  MPI__TYPES_MAP[typeid(T)] = type;
}


/// @brief Create a MPI_Datatype
/// @tparam C Class
template <class C> void MPI__create_type() {
	MPI_Datatype newtype;
  MPI_Type_contiguous(sizeof(C), MPI_CHARACTER, &newtype);
  MPI_Type_commit(&newtype);
  MPI__Type_add<C>(newtype);
}


/// @brief Create alle the MPI_Datatypes for a particle
template <class P> void MPI__create_particle() {
	MPI__create_type<P>();
	MPI__create_type<MPI__Ghost<typename P::Base> >();
	MPI__create_type<MPI__PInMol<P> >();
	MPI__create_type<MPI__Ghost<typename MPI__PInMol<P>::Base > >();
}


void MPI__Init_types();


void MPI__Free_types();

#endif // __MPI_TYPES_UTILS_HPP_INCLUDED
