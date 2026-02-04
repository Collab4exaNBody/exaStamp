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

#include <exanb/core/grid_fields.h>
#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types.h>
#include <onika/math/quaternion.h>
//#include "exanb/particle_specie.h"
#include <exaStamp/potential_factory/pair_potential.h>

#include <onika/math/quaternion_operators.h>

namespace exaStamp
{
  using namespace exanb;

  namespace details
  {
    template<class FS> struct FieldTupleFromFieldSetT;
    template<typename... ids> struct FieldTupleFromFieldSetT< FieldSet<ids...> > { using field_tuple = onika::soatl::FieldTuple<ids...>; };
    template<class FS> using field_tuple_from_field_set = typename FieldTupleFromFieldSetT<FS>::field_tuple;
  }

  struct RigidMoleculeCompute
  {
    using FieldTupleVirialType = details::field_tuple_from_field_set<PairPotentialVirialFieldSet>;
    using FieldTupleType = details::field_tuple_from_field_set<PairPotentialFieldSet>;
    inline void operator () (
      double& ep,
      double& fx,
      double& fy,
      double& fz,
      const Quaternion& orient,
      Vec3d& couple,
      unsigned int nb_atoms,
      const RigidMoleculeAtom* __restrict__ atoms ,
      const FieldTupleType* __restrict__ molecule_atom_values )
    {
      Vec3d force_lab = {0., 0., 0.};
      //matrice de passage base fixe molecule a base labo
      Mat3d mat_bf_lab;
      mat_bf_lab.m11 = orient.w*orient.w + orient.x*orient.x - orient.y*orient.y - orient.z*orient.z;
      mat_bf_lab.m22 = orient.w*orient.w - orient.x*orient.x + orient.y*orient.y - orient.z*orient.z;
      mat_bf_lab.m33 = orient.w*orient.w - orient.x*orient.x - orient.y*orient.y + orient.z*orient.z;
      mat_bf_lab.m21 = 2.0 * (orient.x*orient.y + orient.w*orient.z ); 
      mat_bf_lab.m12 = 2.0 * (orient.x*orient.y - orient.w*orient.z );
      mat_bf_lab.m31 = 2.0 * (orient.x*orient.z - orient.w*orient.y );
      mat_bf_lab.m13 = 2.0 * (orient.x*orient.z + orient.w*orient.y );
      mat_bf_lab.m32 = 2.0 * (orient.y*orient.z + orient.w*orient.x );
      mat_bf_lab.m23 = 2.0 * (orient.y*orient.z - orient.w*orient.x );

      for(unsigned int i=0;i<nb_atoms;i++)
      {
        //passage position sites de force dans repere labo
        Vec3d pos_lab = mat_bf_lab*atoms[i].m_pos;
        force_lab.x = molecule_atom_values[i][field::fx];
        force_lab.y = molecule_atom_values[i][field::fy];
        force_lab.z = molecule_atom_values[i][field::fz];
        ep += molecule_atom_values[i][field::ep];
        fx += molecule_atom_values[i][field::fx];
        fy += molecule_atom_values[i][field::fy];
        fz += molecule_atom_values[i][field::fz];
        couple += cross (pos_lab , force_lab);
      }
    }
  };

}
