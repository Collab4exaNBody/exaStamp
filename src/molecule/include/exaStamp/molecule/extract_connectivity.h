#pragma once

//#include <chrono>
#include <unordered_set>

#include <exanb/core/grid.h>
#include <onika/log.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/grid_fields.h>

#include <exanb/core/std_array_hash.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exaStamp/molecule/id_map.h>

#include <exanb/core/parallel_grid_algorithm.h>

//#define XSTAMP_CHECK_MOL_CONNECTIVITY_DEBUG 1

#ifndef NDEBUG
#ifdef XSTAMP_CHECK_MOL_CONNECTIVITY_DEBUG
#define XSTAMP_CHECK_MOL_CONNECTIVITY 1
#endif
#endif

namespace exaStamp
{

  template<typename GridT> inline  void extract_connectivity_legacy(
    GridT& grid,
   ChemicalBonds& bonds,
   ChemicalAngles& angles,
   ChemicalTorsions& torsions,
   ChemicalImpropers& impropers,
   IdMap& id_map,
   IdMapGhosts& id_map_ghosts)
  {
    using namespace exanb;
    assert( ! omp_in_parallel() );

    size_t n_cells = grid.number_of_cells();
    auto cells = grid.cells();

    bonds    .clear();
    angles   .clear();
    torsions .clear();
    impropers.clear();

//#pragma omp for
    for(size_t cell_i=0;cell_i<n_cells;cell_i++)
    {
      if(grid.is_ghost_cell(cell_i))
        continue;

      for(size_t i =0; i<cells[cell_i].size();++i)
      {
        const uint64_t p_id = cells[cell_i][field::id][i];
        std::array<uint64_t, 4> neigh = cells[cell_i][field::cmol][i];

        // IMPROPERS CASE ---------------------------------------------------------------
        if(neigh[2]!=std::numeric_limits<uint64_t>::max())
        {
          impropers.emplace_back(std::array<uint64_t,4> {p_id, neigh[0], neigh[1], neigh[2]});
        }
        
        if(neigh[3]!=std::numeric_limits<uint64_t>::max())
        {
          lout << "Warning: 4 neighbours not taken into account for impropers: id="<<p_id<<", cmol="<<neigh[0]<<","<<neigh[1]<<","<<neigh[2]<<","<<neigh[3] << std::endl;
        }
        // ------------------------------------------------------------------------------

        // fill bonds, angles, torsions
        // example : entry -> 1 : [2,3,4]
        // bonds  : (1,2) (1,3) (1,4)
        // angles : (2,1,3) (2,1,4) (3,1,4)
        // torsion : to find
        for(size_t a=0; a<neigh.size(); ++a)
        {
          uint64_t id_a = neigh[a];

          // no more neighbour for this atom
          if(id_a==std::numeric_limits<uint64_t>::max()) break;

          size_t cell_a;
          size_t pos_a;

          assert(id_a!=p_id);

          // BONDS CASE -------------------------------------------------------------------
          // here, we sort id to suppress duplicate in bonds after
          // why ? because bond (1,2) involved bond (2,1)
          // so we keep only (1,2)
          if(p_id<id_a) bonds.emplace_back(std::array<uint64_t,2>{p_id, id_a});
          else          bonds.emplace_back(std::array<uint64_t,2>{id_a, p_id});
          // ------------------------------------------------------------------------------


          //-----------------------------PARTICULAR CASE-----------------------------
          // Case for angles
          // one atom with two ghosts neighbour
          // this case is not taken into account
          // (the looking atom can't be a central atom in the following)
          for(size_t z=0; z<cells[cell_i][field::cmol][i].size(); ++z)
          {
            uint64_t id_z = cells[cell_i][field::cmol][i][z];
            if(id_z==std::numeric_limits<uint64_t>::max()) break;
            if( (id_z != id_a) && (id_map_ghosts.find(id_z) != id_map_ghosts.end()) && (id_map_ghosts.find(id_a) != id_map_ghosts.end()) )
            {
              if(id_z<id_a) angles.emplace_back(std::array<uint64_t,3>{id_z,p_id,id_a});
              else          angles.emplace_back(std::array<uint64_t,3>{id_a,p_id,id_z});
            }
          }
          //-------------------------------------------------------------------------


          // get information for second order neighbours
          if(id_map.find(id_a) != id_map.end())
          {
            decode_cell_particle(id_map.at(id_a), cell_a, pos_a);
          }
          // can have several ghosts, but they all have the same id and cmol
          else if(id_map_ghosts.find(id_a) != id_map_ghosts.end())
          {
            decode_cell_particle(id_map_ghosts.find(id_a)->second, cell_a, pos_a);

            // IMPROPERS CASE ---------------------------------------------------------------
            // possible case : one atoms (not the central one) belong to an improper interaction
            // if the central atoms is a ghost in a parallel simulation, we need to take into account the atoms inside the box of simulation
            auto neigh_a = cells[cell_a][field::cmol][pos_a];
            if(neigh_a[2]!=std::numeric_limits<uint64_t>::max())
            {
              impropers.emplace_back(std::array<uint64_t,4> {id_a, neigh_a[0], neigh_a[1], neigh_a[2]});
            }
            if(neigh_a[3]!=std::numeric_limits<uint64_t>::max())
            {
              lerr << "Gestion of the connectivity of atoms : case with 4 neighbours not taken into account for the moment. You must change add impropers in set_connectivity.h." << std::endl;
              std::abort();
            }
            // ------------------------------------------------------------------------------
          }
          else
          {
            lerr << "Set connectivity : " << id_a << " doesn't belong to the box of the MPI process." << std::endl;
            std::abort();
          }

          // second order neighbours
          for(size_t b=0; b<cells[cell_a][field::cmol][pos_a].size(); ++b)
          {
            if(cells[cell_a][field::cmol][pos_a][b] == std::numeric_limits<uint64_t>::max()) break;

            size_t cell_b;
            size_t pos_b;

            uint64_t id_b = cells[cell_a][field::cmol][pos_a][b];

            assert(id_a!=id_b);

            // by definition, atom a is a neighbour of atom b
            if(p_id==id_b) continue;
            assert(p_id!=id_b);
            
            // we order the atoms in order to suppress duplicate after
            if(p_id<id_b)      angles.emplace_back(std::array<uint64_t,3>{p_id,id_a,id_b});
            else if(p_id>id_b) angles.emplace_back(std::array<uint64_t,3>{id_b,id_a,p_id});

            // get information for third order neighbours
            if(id_map.find(id_b) != id_map.end())                    decode_cell_particle(id_map.at(id_b), cell_b, pos_b);
            // can have several ghosts, but they all have the same id and cmol
            else if(id_map_ghosts.find(id_b) != id_map_ghosts.end()) decode_cell_particle(id_map_ghosts.find(id_b)->second, cell_b, pos_b);
            else
            {
              lerr << "Set connectivity : " << id_b << " doesn't belong to the box of the MPI process." << std::endl;
              std::abort();
            }

            //-----------------------------PARTICULAR CASE-----------------------------
            // Second case torsions
            // if the two extrem atoms of the torsion potential are ghosts
            // we need to get the interaction here because
            // they will not be seen in the the next
            for(size_t z=0; z<cells[cell_i][field::cmol][i].size(); ++z)
            {
              uint64_t id_z = cells[cell_i][field::cmol][i][z];
              if(id_z==std::numeric_limits<uint64_t>::max()) break;
              assert(id_z!=id_b);
              if( (id_z != id_a) && (id_map_ghosts.find(id_z) != id_map_ghosts.end()) && (id_map_ghosts.find(id_b) != id_map_ghosts.end()) )
              {
                if(id_z<id_b) torsions.emplace_back(std::array<uint64_t,4> {id_z, p_id, id_a, id_b});
                if(id_z>id_b) torsions.emplace_back(std::array<uint64_t,4> {id_b, id_a, p_id, id_z});
              }
            }
            //-------------------------------------------------------------------------


            // third order neighbours (for torsion)
            for(size_t c=0; c<cells[cell_b][field::cmol][pos_b].size(); ++c)
            {
              uint64_t id_c = cells[cell_b][field::cmol][pos_b][c];
              if(id_c == std::numeric_limits<uint64_t>::max()) break;
              assert(id_c!=id_b);

              if(id_c==p_id) continue; // case of a loop in the molecule. Maybe we should raise an error
              if(id_c==id_a) continue; // by definition, atom b is a neighbour of atom c
              assert(id_c!=id_a);
              
              // we order the atoms to suppress duplicate after
              if(p_id<id_c) torsions.emplace_back(std::array<uint64_t,4> {p_id, id_a, id_b, id_c});
              if(p_id>id_c) torsions.emplace_back(std::array<uint64_t,4> {id_c, id_b, id_a, p_id});
            }

          }
        }

      }
    }

    {
      std::unordered_set< std::array<uint64_t,2> > s;
      for( std::array<uint64_t,2> i : bonds) { s.insert(i); }
      lout << "bonds(L): "<<s.size()<<std::endl;
      bonds.assign( s.begin(), s.end() );
      std::sort( bonds.begin(), bonds.end() );
    }
    
    {
      std::unordered_set< std::array<uint64_t,3> > s;
      for( std::array<uint64_t,3> i : angles) { s.insert(i); }
      angles.assign( s.begin(), s.end() );
      std::sort( angles.begin(), angles.end() );
      lout << "angles(L): "<<s.size()<<std::endl;
    }
    
    {
      std::unordered_set< std::array<uint64_t,4> > s;
      for( std::array<uint64_t,4> i : torsions) { s.insert(i); }
      torsions.assign( s.begin(), s.end() );
      std::sort( torsions.begin(), torsions.end() );
      lout << "torsions(L): "<<s.size()<<std::endl;
    }
    
    {
      std::unordered_set< std::array<uint64_t,4> > s;
      for( std::array<uint64_t,4> i : impropers) { s.insert(i); }
      impropers.assign( s.begin(), s.end() );
      std::sort( impropers.begin(), impropers.end() );
      lout << "impropers(L): "<<s.size()<<std::endl;
    }

  } // end extract_connectivity_legacy(...)


  namespace extract_connectivity_details
  {
    struct PerThreadCount
    {
      size_t bond_count=0;
      size_t angle_count=0;
      size_t torsion_count=0;
      size_t improper_count=0;
    };
    
    template<typename GridT>
    static inline bool is_leader_cell(bool cell_is_inner, const GridT& grid, const IJK& cell_loc, const IJK& domain_cells)
    {
      if( cell_is_inner )
      {
        return true;
      }
      else
      {
        bool is_leader = true;
        for(int di=-1;di<=1;di++) for(int dj=-1;dj<=1;dj++) for(int dk=-1;dk<=1;dk++) if( di!=0 || dj!=0 || dk!=0 )
        {
          IJK loc = cell_loc + IJK{di*domain_cells.i,dj*domain_cells.j,dk*domain_cells.k};
          if( grid.contains(loc) )
          {
            bool loc_is_here = !grid.is_ghost_cell(loc);
            bool lex_order = lexicographic_order(loc,cell_loc);
            if( loc_is_here || lex_order )
            {
              is_leader = false;
            }
          }
        }
        return is_leader;
      }
    }
    
  }
  /*
    Set the connectivity of the molecules
    Add an id to the molecules
  */
  template<typename GridT> inline  void extract_connectivity(
    IJK domain_cells,
    GridT& grid,
    ChemicalBonds& result_bonds,
    ChemicalAngles& result_angles,
    ChemicalTorsions& result_torsions,
    ChemicalImpropers& result_impropers,
    IdMap& id_map,
    IdMapGhosts& id_map_ghosts)
  {
    using namespace extract_connectivity_details;
    static constexpr uint64_t atom_not_found = std::numeric_limits<uint64_t>::max();
//    auto T0 = std::chrono::high_resolution_clock::now();

    assert( ! omp_in_parallel() );

    // size_t n_cells = grid.number_of_cells();
    auto cells = grid.cells_accessor();

    IJK dims = grid.dimension();
    // int gl = grid.ghost_layers();
    
    int max_nt = omp_get_max_threads();
    PerThreadCount per_thread_counts [max_nt];
    ChemicalBonds per_thread_bonds [max_nt];
    ChemicalAngles per_thread_angles [max_nt];
    ChemicalTorsions per_thread_torsions [max_nt];
    ChemicalImpropers per_thread_impropers [max_nt];

    // reuse previously allocated storage if any
    per_thread_bonds[0] = std::move( result_bonds );
    per_thread_bonds[0].clear();
    
    per_thread_angles[0] = std::move( result_angles );
    per_thread_angles[0].clear();
    
    per_thread_torsions[0] = std::move( result_torsions );
    per_thread_torsions[0].clear();

    per_thread_impropers[0] = std::move( result_impropers );
    per_thread_impropers[0].clear();

//    lout << "domain dimensions = "<<domain_cells<<std::endl;
//    lout << "grid dimensions = "<<dims<<std::endl;
//    lout << "grid offset = "<<grid.offset()<<std::endl;
//    lout << "max_nt = "<<max_nt<<std::endl;

    std::set<int64_t> atomIds;
    std::set<int64_t> atomGhostIds;
    const size_t n_cells = grid.number_of_cells();
    for(size_t cell_i=0;cell_i<n_cells;cell_i++)
    {
      size_t n_particles = cells[cell_i].size();
      const auto * __restrict__ ids = cells[cell_i][field::id];
      if( grid.is_ghost_cell(cell_i) ) for(size_t p=0; p<n_particles; ++p) atomGhostIds.insert( ids[p] );
      else for(size_t p=0; p<n_particles; ++p) atomIds.insert( ids[p] );
    }

    auto field_cmol = grid.field_accessor( field::cmol );

#   pragma omp parallel
    {
      int nt = omp_get_num_threads();
      int tid = omp_get_thread_num();
      if(tid>=max_nt || nt>max_nt)
      {
        fatal_error() << "Internal error: bad number of threads" << std::endl;
      }
      auto& bonds = per_thread_bonds[tid];
      auto& angles = per_thread_angles[tid];
      auto& torsions = per_thread_torsions[tid];
      auto& impropers = per_thread_impropers[tid];

      GRID_OMP_FOR_BEGIN(dims,cell_i,cell_loc, schedule(dynamic) nowait )
      {
        bool cell_is_here = ! grid.is_ghost_cell(cell_loc); // 'here' meaning in the central area (not ghost) of the local processor
        // bool cell_is_leader = is_leader_cell(cell_is_here,grid,cell_loc,domain_cells);
        size_t n_particles = cells[cell_i].size();
        const auto * __restrict__ ids = cells[cell_i][field::id];
        const auto * __restrict__ cmols = cells[cell_i].field_pointer_or_null(field_cmol);
        const auto * __restrict__ types = cells[cell_i][field::type];
        
        for(size_t p=0; p<n_particles; ++p)
        {
          uint64_t p_id = ids[p];
          const std::array<uint64_t, 4>& neigh = cmols[p];
          
          if( cell_is_here )
          {
            // first check connectivity consistency (all connex atoms must be found, at least in the ghost area)
            for(int nbh=0;nbh<4;nbh++)
            {
              auto nbh_id = neigh[nbh];
              if(nbh_id != atom_not_found )
              {
                if( atomIds.find(nbh_id)==atomIds.end() && atomGhostIds.find(nbh_id)==atomGhostIds.end() )
                {
#                 pragma omp critical(dbg_mesg)
                  {
                    fatal_error() << "Bad atom connectivity: cell @"<<cell_loc<<" atom #"<<p<<" (id="<<ids[p]<<") : cmol["<<nbh<<"]="<<nbh_id<<" not found in grid"<<std::endl;
                  }
                }
              }
            }
          }
          
          // assert( neigh[3]==std::numeric_limits<uint64_t>::max() );
          const bool p_is_leader = ( encode_cell_particle(cell_i,p,types[p]) == atom_from_idmap(p_id,id_map,id_map_ghosts) );
          if( p_is_leader )
          {        
            // IMPROPERS CASE ---------------------------------------------------------------
            // There are exactly 3 neighbors. In the (common) case where an atom has four neighbors
            // we do not consider (for now) that an improper angle potential can be applied.
            if(neigh[2]!=std::numeric_limits<uint64_t>::max() && neigh[3]==std::numeric_limits<uint64_t>::max())
            {
              // if ghost and a no ghost version exist, ignore this
              // otherwise, detect if we are the preffered ghost cell
              bool n0_is_here = ( id_map.find(neigh[0]) != id_map.end() );
              bool n1_is_here = ( id_map.find(neigh[1]) != id_map.end() );
              bool n2_is_here = ( id_map.find(neigh[2]) != id_map.end() );
              if( cell_is_here || n0_is_here || n1_is_here || n2_is_here )
              {
                impropers.push_back( ChemicalImproper{p_id, neigh[0], neigh[1], neigh[2]} );
              }
            }

            // first neighbor step
            for(auto id_a:neigh)
            {
              assert( id_a != p_id );
              if(id_a==std::numeric_limits<uint64_t>::max()) break;
              bool id_a_is_here = ( id_map.find(id_a) != id_map.end() );
              if( cell_is_here )
              {
                if( p_id < id_a )         bonds.push_back( ChemicalBond{p_id,id_a} );
                else if( ! id_a_is_here ) bonds.push_back( ChemicalBond{id_a,p_id} );
              }
              
              // second neighbor step
              uint64_t enc_atom_a = atom_from_idmap_if_found(id_a,id_map,id_map_ghosts,atom_not_found);
              // can be outside of processor's box if and only if none of preceding atoms in the chain are in the inner box (all are ghosts)
              assert( enc_atom_a!=atom_not_found || !cell_is_here );
              if( enc_atom_a!=atom_not_found )
              {
                size_t cell_a=-1, pos_a=-1;
                decode_cell_particle( enc_atom_a , cell_a, pos_a);
                const auto& neigh_a = cells[cell_a][field_cmol][pos_a];
                for(auto id_b:neigh_a) if(id_b!=p_id)
                {
                  assert( id_b != id_a );
                  if(id_b==std::numeric_limits<uint64_t>::max()) break;                
                  bool id_b_is_here = ( id_map.find(id_b) != id_map.end() );                
                  bool angle_must_be_added = (cell_is_here || id_a_is_here || id_b_is_here) && (p_id!=id_b);
                  bool angle_owner = (cell_is_here && !id_b_is_here) || (cell_is_here==id_b_is_here && p_id<id_b);
                  if( angle_must_be_added && angle_owner )
                  {
                    if( p_id < id_b ) angles.push_back( ChemicalAngle{ p_id , id_a , id_b } );
                    else              angles.push_back( ChemicalAngle{ id_b , id_a , p_id } );
                  }
                  
                  // third neighbor step
                  uint64_t enc_atom_b = atom_from_idmap_if_found(id_b,id_map,id_map_ghosts,atom_not_found);
                  assert( enc_atom_b!=atom_not_found || ( !id_b_is_here && !cell_is_here ) );
                  if( enc_atom_b!=atom_not_found )
                  {
                    size_t cell_b=-1, pos_b=-1;
                    decode_cell_particle( enc_atom_b, cell_b, pos_b);
                    const auto& neigh_b = cells[cell_b][field_cmol][pos_b];
                    for(auto id_c:neigh_b) if(id_c!=id_a && id_c!=p_id)
                    {
                      assert( id_c != id_b );
                      if(id_c==std::numeric_limits<uint64_t>::max()) break;
                      bool id_c_is_here = ( id_map.find(id_c) != id_map.end() );
                      bool add_torsion = (cell_is_here || id_a_is_here || id_b_is_here || id_c_is_here) && (p_id!=id_b) && (id_c!=id_a) && (p_id!=id_c) ;
                      bool torsion_owner = (cell_is_here && !id_c_is_here) || (cell_is_here==id_c_is_here && p_id<id_c);
                      if( add_torsion && torsion_owner )
                      {
                        if( p_id < id_c ) torsions.push_back( ChemicalTorsion{ p_id , id_a , id_b , id_c } );
                        else              torsions.push_back( ChemicalTorsion{ id_c , id_b , id_a , p_id } );

                      }
                    }// for(auto id_c:neigh_b)
                  }// if enc_atom_b!=atom_not_found
                  
                }// for(auto id_b:neigh_a)
              }// if enc_atom_pos_a!=atom_not_found
              
            }

          }
        
        } // if p_is_leader         
      }
      GRID_OMP_FOR_END

      size_t bond_count = bonds.size();
      size_t angle_count = angles.size();
      size_t torsion_count = torsions.size();
      size_t improper_count = impropers.size();
      per_thread_counts[tid] = { bond_count , angle_count , torsion_count , improper_count };
      
#     pragma omp barrier
      
#     pragma omp single nowait
      {
        for(int i=1;i<nt;i++)
        {
          per_thread_counts[i].bond_count     += per_thread_counts[i-1].bond_count;
          per_thread_counts[i].angle_count    += per_thread_counts[i-1].angle_count;
          per_thread_counts[i].torsion_count  += per_thread_counts[i-1].torsion_count;
          per_thread_counts[i].improper_count += per_thread_counts[i-1].improper_count;
        }
        per_thread_bonds    [0].resize( per_thread_counts[nt-1].bond_count );
        per_thread_angles   [0].resize( per_thread_counts[nt-1].angle_count );
        per_thread_torsions [0].resize( per_thread_counts[nt-1].torsion_count );
        per_thread_impropers[0].resize( per_thread_counts[nt-1].improper_count );
      }

#     pragma omp barrier
      
      int bond_offset = 0;
      int angle_offset = 0;
      int torsion_offset = 0;
      int improper_offset = 0;
      if( tid > 0 )
      {
        bond_offset     = per_thread_counts[tid-1].bond_count;
        angle_offset    = per_thread_counts[tid-1].angle_count;
        torsion_offset  = per_thread_counts[tid-1].torsion_count;
        improper_offset = per_thread_counts[tid-1].improper_count;       
      }
      
      for(size_t i=0;i<bond_count    ;i++) per_thread_bonds    [0][bond_offset    +i] = per_thread_bonds    [tid][i];
      for(size_t i=0;i<angle_count   ;i++) per_thread_angles   [0][angle_offset   +i] = per_thread_angles   [tid][i];
      for(size_t i=0;i<torsion_count ;i++) per_thread_torsions [0][torsion_offset +i] = per_thread_torsions [tid][i];
      for(size_t i=0;i<improper_count;i++) per_thread_impropers[0][improper_offset+i] = per_thread_impropers[tid][i];
    }

    result_bonds     = std::move( per_thread_bonds    [0] );
    result_angles    = std::move( per_thread_angles   [0] );
    result_torsions  = std::move( per_thread_torsions [0] );
    result_impropers = std::move( per_thread_impropers[0] );

#   ifdef XSTAMP_CHECK_MOL_CONNECTIVITY
    ChemicalBonds bonds2;
    ChemicalAngles angles2;
    ChemicalTorsions torsions2;
    ChemicalImpropers impropers2;
    extract_connectivity_legacy(grid,bonds2,angles2,torsions2,impropers2,id_map,id_map_ghosts);

    auto& bonds = result_bonds;
    auto& angles = result_angles;
    auto& torsions = result_torsions;
    auto& impropers = result_impropers;

    {
      std::unordered_set< std::array<uint64_t,2> > s;
      for (auto i : bonds)
      {
        if( s.find(i) != s.end() ) { lerr<<"bond not unique as expected : "<<i[0]<<","<<i[1]<<std::endl; std::abort(); }
        s.insert(i);
      }
      lout << "bonds: "<<s.size()<<std::endl;
      bonds.assign( s.begin(), s.end() );
      std::sort( bonds.begin(), bonds.end() );
      if( bonds == bonds2 ) { lout << "bonds equal" << std::endl; }
      else { fatal_error()<<"bonds differ"<<std::endl; }
    }
    
    {
      std::unordered_set< std::array<uint64_t,3> > s;
      for (std::array<uint64_t,3> i : angles)
      {
        if( s.find(i) != s.end() ) { lerr<<"angle not unique as expected : "<<i[0]<<","<<i[1]<<","<<i[2]<<std::endl; std::abort(); }
        s.insert(i);
      }
      lout << "angles: "<<s.size()<<std::endl;
      assert( angles.size() == angles2.size() );
      angles.assign( s.begin(), s.end() );
      std::sort( angles.begin(), angles.end() );
      if( angles == angles2 ) { lout << "angles equal" << std::endl; }
      else { fatal_error()<<"angles differ"<<std::endl; }
    }
    
    {
      std::unordered_set< std::array<uint64_t,4> > s;
      for (std::array<uint64_t,4> i : torsions)
      {
        if( s.find(i) != s.end() ) { lerr<<"torsion not unique as expected : "<<i[0]<<","<<i[1]<<","<<i[2]<<","<<i[3]<<std::endl; std::abort(); }
        s.insert(i);
      }
      lout << "torsions: "<<s.size()<<std::endl;
      assert( torsions.size() == torsions2.size() );
      torsions.assign( s.begin(), s.end() );
      std::sort( torsions.begin(), torsions.end() );
      if( torsions == torsions2 ) { lout << "torsions equal" << std::endl; }
      else { fatal_error()<<"torsions differ"<<std::endl;  }
    }
    
    {
      std::unordered_set< std::array<uint64_t,4> > s;
      for (std::array<uint64_t,4> i : impropers)
      {
        if( s.find(i) != s.end() ) { lerr<<"improper not unique as expected : "<<i[0]<<","<<i[1]<<","<<i[2]<<","<<i[3]<<std::endl; std::abort(); }
        s.insert(i);
      }
      lout << "impropers: "<<s.size()<<std::endl;
      //assert( impropers.size() == impropers2.size() );
      impropers.assign( s.begin(), s.end() );
      std::sort( impropers.begin(), impropers.end() );
      if( impropers == impropers2 ) { lout << "impropers equal" << std::endl; }
      else { lerr<<"impropers differ"<<std::endl; }
    }
#   endif

  } // end extract_connectivity(...)
  
  
} // namespace exaStamp

