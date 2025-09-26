/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#pragma once

#include <unordered_map>
#include <cstdint>
#include <exanb/core/mt_concurrent_map.h>
#include <onika/thread.h>

namespace exaStamp
{
  using IdMap       = exanb::MultiThreadedConcurrentMap< std::unordered_map     <uint64_t, uint64_t> , onika::max_threads_hint *2 >;
  using IdMapGhosts = exanb::MultiThreadedConcurrentMap< std::unordered_multimap<uint64_t, uint64_t> , onika::max_threads_hint *2 >;

  static inline uint64_t atom_from_idmap(const uint64_t id, const IdMap& id_map, const IdMapGhosts& id_map_ghosts)
  {
    auto it = id_map.find(id);
    if(it!=id_map.end()) { return it->second; }
    else
    {
      auto itg = id_map_ghosts.find(id);
      assert( itg != id_map_ghosts.end() );
      return itg->second;
    }
  }

  static inline uint64_t atom_from_idmap_if_found(const uint64_t id, const IdMap& id_map, const IdMapGhosts& id_map_ghosts, const uint64_t value_if_not_found=std::numeric_limits<uint64_t>::max() )
  {
    auto it = id_map.find(id);
    if(it!=id_map.end()) { return it->second; }
    else
    {
      auto itg = id_map_ghosts.find(id);
      if( itg != id_map_ghosts.end() ) return itg->second;
      else return value_if_not_found;
    }
  }

  static inline size_t all_atoms_from_idmap(const uint64_t id, const IdMap& id_map, const IdMapGhosts& id_map_ghosts, uint64_t* locations, size_t max_locations )
  {
    size_t n_locations = 0;
    auto it = id_map.find(id);
    if(it!=id_map.end())
    {
      assert( n_locations < max_locations );
      locations[ n_locations++ ] = it->second;
    }
    const auto range = id_map_ghosts.meta_bucket_map( id ).equal_range( id );
    for(auto it = range.first; it != range.second; ++it)
    {
      assert( n_locations < max_locations );
      locations[ n_locations++ ] = it->second;
    }
    return n_locations;
  }
  
}
