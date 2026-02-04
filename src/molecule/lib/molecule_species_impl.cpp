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

#include <exaStamp/molecule/molecule_species.h>
#include <onika/log.h>
#include <onika/oarray.h>
#include <set>
#include <algorithm>

namespace exaStamp
{
  using namespace exanb;

  void MoleculeSpecies::update_connectivity()
  {
    std::set< onika::oarray_t<int,2> > bonds;
    std::set< onika::oarray_t<int,3> > angles;
    std::set< onika::oarray_t<int,4> > torsions;
    std::set< onika::oarray_t<int,4> > impropers;

    for(int a=0;a<int(m_nb_atoms);a++)
    {
      int bj=0;
      for(int bi=0;bi<4;bi++)
      {
        int b = m_atom_connectivity[a][bi];
        if( b != -1 && b != a )
        {
          if( bi != bj ) { m_atom_connectivity[a][bj] = b; }
          ++bj;
        }
      }
      std::stable_sort( m_atom_connectivity[a] , m_atom_connectivity[a]+bj );
      for(;bj<4;bj++) m_atom_connectivity[a][bj] = -1;
    }
  
    for(int a=0;a<int(m_nb_atoms);a++)
    {
      if( m_atom_connectivity[a][3] != -1 )
      {
        fatal_error() << "Atoms with more than 3 bonds are not supported yet" << std::endl;
      }
      
      if( m_atom_connectivity[a][2] != -1 )
      {
        impropers.insert( { a , m_atom_connectivity[a][0] , m_atom_connectivity[a][1] , m_atom_connectivity[a][2] } );
      }
      
      for(int bi=0;bi<4;bi++)
      {
        int b = m_atom_connectivity[a][bi];
        if( b != -1 && b != a )
        {          
          // 1st order neighbor, we have a bond
          if(a<b) bonds.insert( {a,b} );
          else bonds.insert( {b,a} );
          
          for(int ci=0;ci<4;ci++)
          {
            int c = m_atom_connectivity[b][ci];
            if( c != -1 && c != a && c != b )
            {
              // 2nd order neighbor
              if(a<c) angles.insert( {a,b,c} );
              else angles.insert( {c,b,a} );
              
              for(int di=0;di<4;di++)
              {
                int d = m_atom_connectivity[c][di];
                if( d != -1 && d != a && d != b && d != c )
                {
                  // 3rd order neighbor
                  if(a<d) torsions.insert( {a,b,c,d} );
                  else torsions.insert( {d,c,b,a} );
                }
              }
            }
          }     
        }
      }
    }

    bool fatal_issue = false;
    if( bonds.size()     > MAX_MOLECULE_BONDS     ) { lerr << "Number of bonds ("    <<bonds.size()     <<") exceeds max number of intramolecular bonds ("    <<MAX_MOLECULE_BONDS    <<")"<<std::endl; fatal_issue = true; }
    if( angles.size()    > MAX_MOLECULE_BENDS     ) { lerr << "Number of angles ("   <<angles.size()    <<") exceeds max number of intramolecular angles ("   <<MAX_MOLECULE_BENDS    <<")"<<std::endl; fatal_issue = true; }
    if( torsions.size()  > MAX_MOLECULE_TORSIONS  ) { lerr << "Number of torsions (" <<torsions.size()  <<") exceeds max number of intramolecular torsions (" <<MAX_MOLECULE_TORSIONS <<")"<<std::endl; fatal_issue = true; }
    if( impropers.size() > MAX_MOLECULE_IMPROPERS ) { lerr << "Number of impropers ("<<impropers.size() <<") exceeds max number of intramolecular impropers ("<<MAX_MOLECULE_IMPROPERS<<")"<<std::endl; fatal_issue = true; }
    if( fatal_issue ) { fatal_error()<<"Fatal issue while reading intramolecular data"<<std::endl; }

    m_nb_bonds = 0;
    m_nb_bends = 0;
    m_nb_torsions = 0;
    m_nb_impropers = 0;

    for(const auto& x:bonds)     { assert(m_nb_bonds    <MAX_MOLECULE_BONDS    ); for(int i=0;i<2;i++) { m_bonds    [m_nb_bonds    ][i]=x[i]; } ++m_nb_bonds; }
    for(const auto& x:angles)    { assert(m_nb_bends    <MAX_MOLECULE_BENDS    ); for(int i=0;i<3;i++) { m_bends    [m_nb_bends    ][i]=x[i]; } ++m_nb_bends; }
    for(const auto& x:torsions)  { assert(m_nb_torsions <MAX_MOLECULE_TORSIONS ); for(int i=0;i<4;i++) { m_torsions [m_nb_torsions ][i]=x[i]; } ++m_nb_torsions; }
    for(const auto& x:impropers) { assert(m_nb_impropers<MAX_MOLECULE_IMPROPERS); for(int i=0;i<4;i++) { m_impropers[m_nb_impropers][i]=x[i]; } ++m_nb_impropers; }
  }

}

