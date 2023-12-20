/// @file 
/// @brief Various utilities to handle neighborhood of a rectangular box (faces, edges and vertexes)

#ifndef __NEIGHBOR_HPP_INCLUDED
#define __NEIGHBOR_HPP_INCLUDED


#include "utils/vec3/vec3.hpp"


/// @brief Namespace containing constants and functions to handle neighborhoods of a rectangular box
/// @note Possible synonyms for neighborhood : side, neighbor, direction
namespace Neighbor {

  constexpr uint8_t num_neighbors = 26;  ///< Total number of neighborhoods
  constexpr uint8_t num_faces     =  6;  ///< Total number of faces
  constexpr uint8_t num_edges     = 12;  ///< Total number of edges
  constexpr uint8_t num_vertices  =  8;  ///< Total number of vertexes

  constexpr uint8_t begin = 0;              ///< Start index 
  constexpr uint8_t end   = num_neighbors;  ///< End index

  constexpr uint8_t begin_face = 0;                       ///< Start index for faces
  constexpr uint8_t end_face   = begin_face + num_faces;  ///< End index for faces

  constexpr uint8_t begin_edge = end_face;                ///< Start index for edges
  constexpr uint8_t end_edge   = begin_edge + num_edges;  ///< End  index for edges

  constexpr uint8_t begin_vertex = end_edge;                     ///< Start index for vertexes
  constexpr uint8_t end_vertex   = begin_vertex + num_vertices;  ///< End index for vertexes

  constexpr int8_t null = -1;  ///< Value for no neighborhood
  constexpr int8_t init = -2;  ///< Init value for neighborhood index (debug purpose)

  inline vec3<int> getCoords(int idx);
  
  inline int getIndex(const vec3<int>& c);

  inline int getIndex(int cx, int cy, int cz);

  inline int computeNumberOfGhostCells(int idx, vec3<int>& numberOfCells);

  constexpr int8_t xCoords[num_neighbors] = { -1, +1,  0,  0,  0,  0, -1, -1, +1, +1, -1, -1, +1, +1,  0,  0,  0,  0, -1, -1, -1, -1, +1, +1, +1, +1 }; ///< X coordinate for each index
  constexpr int8_t yCoords[num_neighbors] = {  0,  0, -1, +1,  0,  0, -1, +1, -1, +1,  0,  0,  0,  0, -1, -1, +1, +1, -1, -1, +1, +1, -1, -1, +1, +1 }; ///< Y coordinate for each index
  constexpr int8_t zCoords[num_neighbors] = {  0,  0,  0,  0, -1, +1,  0,  0,  0,  0, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1, -1, +1 }; ///< Z coordinate for each index

  constexpr int8_t indexes[num_neighbors+1] = { 18, 6, 19, 10, 0, 11, 20, 7, 21, 14, 2, 15, 4, null, 5, 16, 3, 17, 22, 8, 23, 12, 1, 13, 24, 9, 25 }; ///< Index for each set of coordinates

};


/// @brief Get neighborhood coordinates from index
/// @param [in] idx Index (in range [@c begin , @c end])
/// @return A @c vec3<int> (which values are -1, 0 or 1) giving relative
/// direction of the given neighborhood
inline vec3<int> Neighbor::getCoords(int idx) {
  if (idx<0 || idx>=num_neighbors) return vec3<int>(0);
  return vec3<int>(xCoords[idx], yCoords[idx], zCoords[idx]);
}


/// @brief Get index from 3D coordinates
///
/// Inverse function of @c Neighbor::getCoords()
/// @param [in] c 3D coordinates
/// @return Neighborhood index
inline int Neighbor::getIndex(const vec3<int>& c) {
  int base3 = c.x*9 + c.y*3 + c.z + 13;
  if (base3<0 || base3>num_neighbors) return null;
  return indexes[base3];
}


/// @brief Get index from 3 1D coordinates (not used)
///
/// Inverse function of @c Neighbor::getCoords()
/// @param [in] cx X coordinate
/// @param [in] cy Y coordinate
/// @param [in] cz Z coordinate
/// @return Neighborhood index
inline int Neighbor::getIndex(int cx, int cy, int cz) {
  return getIndex(std::move(vec3<int>(cx, cy, cz)));
}


/// @brief Compute number of ghost cells ? (not used)
/// @param [in] numberOfCells Number of cells
/// @param [in] idx Index of the neighborhood
/// @return Calculated 0
inline int Neighbor::computeNumberOfGhostCells(int idx, vec3<int>& numberOfCells) {
  return product((1-auxAbs(getCoords(idx)))*numberOfCells);
}

#endif // __NEIGHBOR_HPP_INCLUDED
