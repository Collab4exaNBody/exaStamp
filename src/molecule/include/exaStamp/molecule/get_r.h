#pragma once

#include <cmath>

#include <exanb/core/basic_types.h>
//#include "exanb/container_utils.h"
#include <exaStamp/molecule/id_map.h>

namespace exaStamp
{
  template<typename CellType>
  static inline std::array<double,3> get_r(const CellType& cells,
                                       const Vec3d& size_box,
                                       const double half_min_size_box,
                                       const IdMap& id_map,
                                       const IdMapGhosts& id_map_ghosts,
                                       const uint64_t id1, const uint64_t id2,
                                       size_t& cell1, size_t& cell2,
                                       size_t& pos1, size_t& pos2,
                                       unsigned int& type1, unsigned int& type2)
  {
    assert(id_map.find(id1)!=id_map.end() || id_map_ghosts.find(id1) != id_map_ghosts.end());
    assert(id_map.find(id2)!=id_map.end() || id_map_ghosts.find(id2) != id_map_ghosts.end());

    uint64_t atom_to_decode_a = atom_from_idmap(id1,id_map,id_map_ghosts);
    uint64_t atom_to_decode_b = atom_from_idmap(id2,id_map,id_map_ghosts);

    // Decode informations (cells, pos and types)
    decode_cell_particle(atom_to_decode_a, cell1, pos1, type1);
    decode_cell_particle(atom_to_decode_b, cell2, pos2, type2);

    Vec3d rij = { cells[cell2][field::rx][pos2] - cells[cell1][field::rx][pos1]
                , cells[cell2][field::ry][pos2] - cells[cell1][field::ry][pos1]
                , cells[cell2][field::rz][pos2] - cells[cell1][field::rz][pos1] };

    if(rij.x>  half_min_size_box) rij.x -= size_box.x;
    if(rij.x< -half_min_size_box) rij.x += size_box.x;
    if(rij.y>  half_min_size_box) rij.y -= size_box.y;
    if(rij.y< -half_min_size_box) rij.y += size_box.y;
    if(rij.z>  half_min_size_box) rij.z -= size_box.z;
    if(rij.z< -half_min_size_box) rij.z += size_box.z;
    assert(norm(rij)<half_min_size_box);

    assert( (rij != Vec3d{0,0,0}) );
    assert(atom_to_decode_a != std::numeric_limits<uint64_t>::max());
    assert(atom_to_decode_b != std::numeric_limits<uint64_t>::max());

    return std::array<double,3>{ rij.x , rij.y , rij.z };
  }


}
