#pragma once

#include <vector>
#include <map>

namespace exaStamp
{


  struct MolGraphNode
  {
    int m_atom_type;
    int m_label;

    std::vector<int> m_bonds;
  };

}

