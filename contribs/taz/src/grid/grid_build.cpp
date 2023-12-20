/// @file
/// @brief Implementation of functions used to build the grid


#include "utils/array/array.hpp"
#include "utils/vec3/vec3.hpp"


/// @brief For a rectangular box neighborhood (edge, face or vertex) get the direction where neighbors can be found
/// @param [in] mov Neighborhood
/// @return Directions of the neighbors
Array< vec3<int> > __getMovList(const vec3<int>& mov) {

  int base3 = mov.x*100 + mov.y*10 + mov.z + 111;

  Array< vec3<int> > result;

  auto init = [&] (Array< vec3<int> >& array, int n) -> void {
    array = Array< vec3<int> >(n);
  };

  switch (base3) {

    // FACES
    
  case  11: 
  case 211: 
  case 101: 
  case 121: 
  case 110: 
  case 112: 
    init(result, 1); 
    result[0] = mov;
    break;
    
    // EDGES

  case   1: 
    init(result, 3);
    result[0] = vec3<int>(-1,0,0);
    result[1] = vec3<int>(0,-1,0);
    result[2] = mov;
    break;

  case  21: 
    init(result, 3);
    result[0] = vec3<int>(-1,0,0);
    result[1] = vec3<int>(0,1,0);
    result[2] = mov;
    break;

  case 201: 
    init(result, 3);
    result[0] = vec3<int>(1,0,0);
    result[1] = vec3<int>(0,-1,0);
    result[2] = mov;
    break;

  case 221: 
    init(result, 3);
    result[0] = vec3<int>(1,0,0);
    result[1] = vec3<int>(0,1,0);
    result[2] = mov;
    break;

  case  10: 
    init(result, 3);
    result[0] = vec3<int>(-1,0,0);
    result[1] = vec3<int>(0,0,-1);
    result[2] = mov;
    break;

  case  12: 
    init(result, 3);
    result[0] = vec3<int>(-1,0,0);
    result[1] = vec3<int>(0,0,1);
    result[2] = mov;
    break;

  case 210: 
    init(result, 3);
    result[0] = vec3<int>(1,0,0);
    result[1] = vec3<int>(0,0,-1);
    result[2] = mov;
    break;

  case 212: 
    init(result, 3);
    result[0] = vec3<int>(1,0,0);
    result[1] = vec3<int>(0,0,1);
    result[2] = mov;
    break;

  case 100: 
    init(result, 3);
    result[0] = vec3<int>(0,-1,0);
    result[1] = vec3<int>(0,0,-1);
    result[2] = mov;
    break;

  case 102: 
    init(result, 3);
    result[0] = vec3<int>(0,-1,0);
    result[1] = vec3<int>(0,0,1);
    result[2] = mov;
    break;

  case 120: 
    init(result, 3);
    result[0] = vec3<int>(0,1,0);
    result[1] = vec3<int>(0,0,-1);
    result[2] = mov;
    break;

  case 122: 
    init(result, 3);
    result[0] = vec3<int>(0,1,0);
    result[1] = vec3<int>(0,0,1);
    result[2] = mov;
    break;

    // VERTEXES

  case   0: 
    init(result, 7);
    result[0] = vec3<int>(-1,0,0);
    result[1] = vec3<int>(0,-1,0);
    result[2] = vec3<int>(0,0,-1);
    result[3] = vec3<int>(-1,-1,0);
    result[4] = vec3<int>(-1,0,-1);
    result[5] = vec3<int>(0,-1,-1);
    result[6] = mov;
    break;

  case   2: 
    init(result, 7);
    result[0] = vec3<int>(-1,0,0);
    result[1] = vec3<int>(0,-1,0);
    result[2] = vec3<int>(0,0,1);
    result[3] = vec3<int>(-1,-1,0);
    result[4] = vec3<int>(-1,0,1);
    result[5] = vec3<int>(0,-1,1);
    result[6] = mov;
    break;

  case  20: 
    init(result, 7);
    result[0] = vec3<int>(-1,0,0);
    result[1] = vec3<int>(0,1,0);
    result[2] = vec3<int>(0,0,-1);
    result[3] = vec3<int>(-1,1,0);
    result[4] = vec3<int>(-1,0,-1);
    result[5] = vec3<int>(0,1,-1);
    result[6] = mov;
    break;

  case  22: 
    init(result, 7);
    result[0] = vec3<int>(-1,0,0);
    result[1] = vec3<int>(0,1,0);
    result[2] = vec3<int>(0,0,1);
    result[3] = vec3<int>(-1,1,0);
    result[4] = vec3<int>(-1,0,1);
    result[5] = vec3<int>(0,1,1);
    result[6] = mov;
    break;

  case 200: 
    init(result, 7);
    result[0] = vec3<int>(1,0,0);
    result[1] = vec3<int>(0,-1,0);
    result[2] = vec3<int>(0,0,-1);
    result[3] = vec3<int>(1,-1,0);
    result[4] = vec3<int>(1,0,-1);
    result[5] = vec3<int>(0,-1,-1);
    result[6] = mov;
    break;

  case 202: 
    init(result, 7);
    result[0] = vec3<int>(1,0,0);
    result[1] = vec3<int>(0,-1,0);
    result[2] = vec3<int>(0,0,1);
    result[3] = vec3<int>(1,-1,0);
    result[4] = vec3<int>(1,0,1);
    result[5] = vec3<int>(0,-1,1);
    result[6] = mov;
    break;

  case 220: 
    init(result, 7);
    result[0] = vec3<int>(1,0,0);
    result[1] = vec3<int>(0,1,0);
    result[2] = vec3<int>(0,0,-1);
    result[3] = vec3<int>(1,1,0);
    result[4] = vec3<int>(1,0,-1);
    result[5] = vec3<int>(0,1,-1);
    result[6] = mov;
    break;

  case 222: 
    init(result, 7);
    result[0] = vec3<int>(1,0,0);
    result[1] = vec3<int>(0,1,0);
    result[2] = vec3<int>(0,0,1);
    result[3] = vec3<int>(1,1,0);
    result[4] = vec3<int>(1,0,1);
    result[5] = vec3<int>(0,1,1);
    result[6] = mov;
    break;

  default:
    break;

  }

  return result;

}
