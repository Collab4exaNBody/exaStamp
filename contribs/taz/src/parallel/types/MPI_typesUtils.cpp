/// @file 
/// @brief Implementation for the general functions concerning all the MPI_Datatypes


#include <map>

#include "io/StampV3LegacyIOStructures.hpp"

#include "parallel/types/MPI_cell.hpp"
#include "parallel/types/MPI_exchangeEAM.hpp"
#include "parallel/types/MPI_exchangeGhostV3.hpp"
#include "parallel/types/MPI_particle.hpp"
#include "parallel/types/MPI_mat3.hpp"
#include "parallel/types/MPI_mesoparticle.hpp"
#include "parallel/types/MPI_smoothparticle.hpp"
#include "parallel/types/MPI_typesUtils.hpp"
#include "parallel/types/MPI_vec3.hpp"


// Definition of some MPI types for old MPI versions
#if MPI_VERSION < 2 || MPI_SUBVERSION < 2

#ifndef MPI_INT8_T
/// @brief MPI type for an integer on 8 bytes
#define MPI_INT8_T MPI_CHAR
#endif

#ifndef MPI_INT16_T
/// @brief MPI type for an integer on 16 bytes
#define MPI_INT16_T MPI_SHORT
#endif

#ifndef MPI_INT32_T
/// @brief MPI type for an integer on 32 bytes
#define MPI_INT32_T MPI_INT
#endif

#ifndef MPI_INT64_T
/// @brief MPI type for an integer on 64 bytes
#define MPI_INT64_T MPI_LONG
#endif

#ifndef MPI_UINT8_T
/// @brief MPI type for an unsigned integer on 8 bytes
#define MPI_UINT8_T MPI_UNSIGNED_CHAR
#endif

#ifndef MPI_UINT16_T
/// @brief MPI type for an unsigned integer on 16 bytes
#define MPI_UINT16_T MPI_SHORT
#endif

#ifndef MPI_UINT32_T
/// @brief MPI type for an unsigned integer on 32 bytes
#define MPI_UINT32_T MPI_UNSIGNED
#endif

#ifndef MPI_UINT64_T
/// @brief MPI type for an unsigned integer on 64 bytes
#define MPI_UINT64_T MPI_UNSIGNED_LONG_LONG
#endif

#endif


extern std::map<std::type_index, MPI_Datatype> MPI__TYPES_MAP;


/// @brief 'Initialization of the MPI_Datatype database
void MPI__Init_types() {

  MPI__Type_add<char>(MPI_CHAR);
  MPI__Type_add<int>(MPI_INT);
  MPI__Type_add<long>(MPI_LONG);

  MPI__Type_add<unsigned char>(MPI_UNSIGNED_CHAR);
  MPI__Type_add<uint>(MPI_UNSIGNED);
  MPI__Type_add<ulong>(MPI_UNSIGNED_LONG);

  MPI__Type_add<int8_t >(MPI_INT8_T );
  MPI__Type_add<int16_t>(MPI_INT16_T);
  MPI__Type_add<int32_t>(MPI_INT32_T);
  MPI__Type_add<int64_t>(MPI_INT64_T);

  MPI__Type_add<uint8_t >(MPI_UINT8_T );
  MPI__Type_add<uint16_t>(MPI_UINT16_T);
  MPI__Type_add<uint32_t>(MPI_UINT32_T);
  MPI__Type_add<uint64_t>(MPI_UINT64_T);

  MPI__Type_add<float> (MPI_FLOAT);
  MPI__Type_add<double>(MPI_DOUBLE);

  create_MPI_type_vec3<int>();
  create_MPI_type_vec3<float>();
  create_MPI_type_vec3<double>();

  create_MPI_type_mat3<int>();
  create_MPI_type_mat3<float>();
  create_MPI_type_mat3<double>();

  MPI__create_particle<MPI__Particle>();
  MPI__create_particle<MPI__Mesoparticle>();
  MPI__create_particle<MPI__SmoothParticle>();

  MPI__create_type<ExchangeEAM>();
  MPI__create_type<ExchangeGhostV3>();

  MPI__create_type<MPI__Ghost<int>>();
  MPI__create_type<MPI__CellEnv>();

  MPI__create_type<LegacyHeaderIOStruct>();
  MPI__create_type<LegacyParticleIOStruct>();
  MPI__create_type<LegacyDPDEParticleIOStruct>();

}


//
/// @brief Delete all MPI types
void MPI__Free_types() {

  // Erase predefined types because it's not possible to free them

  MPI__TYPES_MAP.erase(typeid(char));
  MPI__TYPES_MAP.erase(typeid(int));
  MPI__TYPES_MAP.erase(typeid(long));

  MPI__TYPES_MAP.erase(typeid(unsigned char));
  MPI__TYPES_MAP.erase(typeid(uint));
  MPI__TYPES_MAP.erase(typeid(ulong));

  MPI__TYPES_MAP.erase(typeid(int8_t));
  MPI__TYPES_MAP.erase(typeid(int16_t));
  MPI__TYPES_MAP.erase(typeid(int32_t));
  MPI__TYPES_MAP.erase(typeid(int64_t));

  MPI__TYPES_MAP.erase(typeid(uint8_t));
  MPI__TYPES_MAP.erase(typeid(uint16_t));
  MPI__TYPES_MAP.erase(typeid(uint32_t));
  MPI__TYPES_MAP.erase(typeid(uint64_t));

  MPI__TYPES_MAP.erase(typeid(float));
  MPI__TYPES_MAP.erase(typeid(double));

  for (auto& elem : MPI__TYPES_MAP) 
    MPI_Type_free(&elem.second);

  MPI__TYPES_MAP.clear();

}
