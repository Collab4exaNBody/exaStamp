#pragma once

#include <exanb/core/grid.h>
#include <exanb/compute/compute_pair_buffer.h>
#include <exanb/compute/compute_pair_optional_args.h>
#include <exanb/compute/compute_pair_function.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/field_sets.h>
#include <exanb/core/thread.h>
#include <exanb/core/particle_type_pair.h>
#include <onika/soatl/field_id.h>
#include <onika/soatl/field_pointer_tuple.h>
#include <onika/soatl/field_tuple.h>

#include <exanb/particle_neighbors/chunk_neighbors_iterator.h>

#include <vector>
#include <functional>
#include <memory>

#include <exanb/core/quaternion_operators.h>
#include <exanb/core/quaternion_to_matrix.h>

namespace exaStamp
{

  struct IdentityPairIdMap
  {
    inline constexpr int operator [] (int i) { return i; }
  };

  template< typename GridT,
            typename ComputePairBufferT,
            typename MoleculeFuncT,
            typename OptionalArgsT,
            typename... field_ids>
  static inline void compute_pair_rigidmol_cell(
    GridT& grid,
    ComputePairBufferT tab[],
    size_t cell_a,
    const std::vector< double >& rcuts,
    const std::vector< std::shared_ptr< ComputePairOperator< FieldSet<field_ids...> > > >& funcs,
    const unsigned int NTYPES,
    const int * __restrict__ pair_id_map,
    const ParticleSpecies& species,
    MoleculeFuncT molecule_func,
    const OptionalArgsT& optional, // locks are needed if symmetric computation is enabled
    FieldSet< field_ids... >
    )
  {
    static constexpr auto const_1 = std::integral_constant<unsigned int,1>();
    static constexpr auto const_4 = std::integral_constant<unsigned int,4>();
    static constexpr auto const_8 = std::integral_constant<unsigned int,8>();

    const unsigned int chunk_size = optional.nbh.m_chunk_size;
    switch( chunk_size )
    {
      case 1 : compute_pair_rigidmol_cell(grid,tab,cell_a,rcuts,funcs,NTYPES,pair_id_map,species,molecule_func,optional, const_1, FieldSet<field_ids...>{}); break;
      case 4 : compute_pair_rigidmol_cell(grid,tab,cell_a,rcuts,funcs,NTYPES,pair_id_map,species,molecule_func,optional, const_4, FieldSet<field_ids...>{}); break;
      case 8 : compute_pair_rigidmol_cell(grid,tab,cell_a,rcuts,funcs,NTYPES,pair_id_map,species,molecule_func,optional, const_8, FieldSet<field_ids...>{}); break;
      default:
        lerr << "compute_cell_particle_pairs: chunk size "<<chunk_size<<" not supported. Accepted values are 1, 4, 8." << std::endl;
        std::abort();
        break;
    }

  }

  template< typename GridT,
            typename CST,
            typename ComputePairBufferT,
            typename MoleculeFuncT,
            typename OptionalArgsT,
            typename... field_ids>
  static inline void compute_pair_rigidmol_cell(
    GridT& grid,
    ComputePairBufferT tab[],
    size_t cell_a,
    const std::vector< double >& rcuts,
    const std::vector< std::shared_ptr< ComputePairOperator< FieldSet<field_ids...> > > >& funcs,
    const unsigned int NTYPES,
    const int * __restrict__  pair_id_map,
    const ParticleSpecies& species,
    MoleculeFuncT molecule_func,
    const OptionalArgsT& optional,
    CST CS,
    FieldSet< field_ids... >
    )
  {
    using exanb::chunknbh_stream_to_next_particle;
    using exanb::chunknbh_stream_info;
    using exanb::decode_cell_index;

    using namespace exanb;
    using BaseFields = FieldSet<field::_rx,field::_ry,field::_rz,field::_type,field::_orient>;
    using BaseFieldsCouple = FieldSet<field::_rx,field::_ry,field::_rz,field::_type,field::_orient,field::_couple>;
    using CellAFields = MergeFieldSet< BaseFieldsCouple , FieldSet<field_ids...> >;
    using CellAPointerTuple = GridFieldSetPointerTuple<GridT,CellAFields>;
    using CellBPointerTuple = GridFieldSetPointerTuple<GridT,BaseFields>;
    using FieldTupleType = onika::soatl::FieldTuple< field_ids... >;

    static constexpr bool Symetric = false;
    // const unsigned int N_PAIRS = unique_pair_count(NTYPES);
    const unsigned int n_reduced_pairs = rcuts.size();

    assert( rcuts.size() == funcs.size() );
    assert( rcuts.size() <= unique_pair_count(NTYPES) );
//    assert( funcs.size() <= N_PAIRS );

    FieldTupleType molecule_atom_values[MAX_RIGID_MOLECULE_ATOMS];

    auto cells = grid.cells();
    IJK dims = grid.dimension();
    IJK loc_a = grid_index_to_ijk( dims, cell_a );

    const size_t cell_a_particles = cells[cell_a].size();
    CellAPointerTuple cell_a_pointers;
    cells[cell_a].capture_pointers(cell_a_pointers);
    // auto& cell_a_locks = optional.locks[cell_a];

    auto nbh_data_ctx = optional.nbh_data.make_ctx();

    const double*     __restrict__ rx_a       = cell_a_pointers[field::rx];     ONIKA_ASSUME_ALIGNED(rx_a);
    const double*     __restrict__ ry_a       = cell_a_pointers[field::ry];     ONIKA_ASSUME_ALIGNED(ry_a);
    const double*     __restrict__ rz_a       = cell_a_pointers[field::rz];     ONIKA_ASSUME_ALIGNED(rz_a);
    const uint8_t*    __restrict__ type_a_ptr = cell_a_pointers[field::type];   ONIKA_ASSUME_ALIGNED(type_a_ptr);
    const Quaternion* __restrict__ orient_a   = cell_a_pointers[field::orient]; ONIKA_ASSUME_ALIGNED(orient_a);
          Vec3d*      __restrict__ couple_a   = cell_a_pointers[field::couple]; ONIKA_ASSUME_ALIGNED(couple_a);

    for(unsigned int i=0;i<n_reduced_pairs;i++)
    {
      tab[i].cell = cell_a;
      tab[i].count = 0;
    }

    const auto stream_info = chunknbh_stream_info( optional.nbh.m_nbh_streams[cell_a] , cell_a_particles );
    const uint16_t* stream_base = stream_info.stream;
    const uint16_t* __restrict__ stream = stream_base;
    //const uint32_t* __restrict__ particle_offset = stream_info.offset;
    //const int32_t poffshift = stream_info.shift;

    for(size_t p_a=0;p_a<cell_a_particles;p_a++)
    {
      const double* __restrict__ rx_b = nullptr; ONIKA_ASSUME_ALIGNED(rx_b);
      const double* __restrict__ ry_b = nullptr; ONIKA_ASSUME_ALIGNED(ry_b);
      const double* __restrict__ rz_b = nullptr; ONIKA_ASSUME_ALIGNED(rz_b);
      const Quaternion* __restrict__ orient_b = nullptr; ONIKA_ASSUME_ALIGNED(orient_b);
      const uint8_t* __restrict__ type_b_ptr = nullptr; ONIKA_ASSUME_ALIGNED(type_b_ptr);

      const unsigned int mol_type_a = type_a_ptr[p_a];
      assert( mol_type_a < species.size() );

      for(unsigned int i=0;i<n_reduced_pairs;i++)
      {
        tab[i].part = p_a;
//        tab[i].ta = type_a;
      }

      const Vec3d mola_r = optional.xform.transformCoord( Vec3d{ rx_a[p_a] , ry_a[p_a] , rz_a[p_a] } );
      //dr = optional.xform.transformCoord( dr );
      const unsigned int nb_sub_atom_a = species[mol_type_a].m_rigid_atom_count;
      const auto * rigid_atoms_a = species[mol_type_a].m_rigid_atoms;
      //molecule_atom_values.resize( nb_sub_atom_a );

      const auto stream_backup = stream;
      //const auto p_nbh_index_backup = p_nbh_index;

      // needs to take rotation into account here
      const Mat3d mat_bf_lab_a = quaternion_to_matrix( orient_a[p_a] );

      for(unsigned int cur_sub_atom_a=0; cur_sub_atom_a<nb_sub_atom_a; cur_sub_atom_a++)
      {
        // initialize neighbor list traversal
        size_t p_nbh_index = 0;

        stream = stream_backup;
        //p_nbh_index = p_nbh_index_backup;
        molecule_atom_values[cur_sub_atom_a].zero(); // zero all fields before accumulating forces and energies
        
        const Vec3d ra = mola_r + mat_bf_lab_a*rigid_atoms_a[cur_sub_atom_a].m_pos;
        const unsigned int type_a = rigid_atoms_a[cur_sub_atom_a].m_atom_type;
        assert( type_a < NTYPES );

        unsigned int cell_groups = *(stream++); // number of cell groups for this neighbor list
        size_t cell_b = cell_a;
        unsigned int chunk = 0;
        unsigned int nchunks = 0;
        unsigned int cg = 0; // cell group index.
        bool symcont = true;

        for(cg=0; cg<cell_groups && symcont ;cg++)
        {
          uint16_t cell_b_enc = *(stream++);
          IJK loc_b = loc_a + decode_cell_index(cell_b_enc);
          cell_b = grid_ijk_to_index( dims , loc_b );
          unsigned int nbh_cell_particles = cells[cell_b].size();
          CellBPointerTuple cell_b_pointers;
          cells[cell_b].capture_pointers(cell_b_pointers);
          rx_b       = cell_b_pointers[field::rx]; ONIKA_ASSUME_ALIGNED(rx_b);
          ry_b       = cell_b_pointers[field::ry]; ONIKA_ASSUME_ALIGNED(ry_b);
          rz_b       = cell_b_pointers[field::rz]; ONIKA_ASSUME_ALIGNED(rz_b);
          orient_b   = cell_b_pointers[field::orient]; ONIKA_ASSUME_ALIGNED(orient_b);
          type_b_ptr = cell_b_pointers[field::type]; ONIKA_ASSUME_ALIGNED(type_b_ptr);

          nchunks = *(stream++);
          for(chunk=0;chunk<nchunks && symcont;chunk++)
          {
            unsigned int chunk_start = static_cast<unsigned int>( *(stream++) ) * CS;
            for(unsigned int i=0;i<CS && symcont;i++)
            {
              unsigned int p_b = chunk_start + i;
              if( Symetric && ( cell_b>cell_a || ( cell_b==cell_a && p_b>=p_a ) ) )
              {
                symcont = false;
              }
              else if( p_b<nbh_cell_particles && (cell_b!=cell_a || p_b!=p_a) )
              {
                const unsigned int mol_type_b = type_b_ptr[p_b];
                assert( mol_type_b < species.size() );

                const Vec3d molb_r = optional.xform.transformCoord( Vec3d{ rx_b[p_b] , ry_b[p_b] , rz_b[p_b] } );
                const unsigned int nb_sub_atom_b = species[mol_type_b].m_rigid_atom_count;
                const auto * rigid_atoms_b = species[mol_type_b].m_rigid_atoms;
                const Mat3d mat_bf_lab_b = quaternion_to_matrix( orient_b[p_b] );

                for(unsigned int cur_sub_atom_b=0; cur_sub_atom_b<nb_sub_atom_b; cur_sub_atom_b++)
                {
                  const Vec3d rb = molb_r + mat_bf_lab_b*rigid_atoms_b[cur_sub_atom_b].m_pos;
//                  const Vec3d rb = molb_r + rotate( orient_b[p_b] , rigid_atoms_b[cur_sub_atom_b].m_pos );
                  const unsigned int type_b = rigid_atoms_b[cur_sub_atom_b].m_atom_type;
                  assert( type_b < NTYPES );
                  const Vec3d dr = rb - ra; // optional.xform.transformCoord( rb - ra );

                  unsigned int pair_ab_id = unique_pair_id(type_a,type_b);
                  assert( pair_ab_id < unique_pair_count(NTYPES) );
                  if(pair_id_map!=nullptr) pair_ab_id = pair_id_map[pair_ab_id];
                  assert( pair_ab_id < n_reduced_pairs );

                  const double rcut = rcuts[pair_ab_id];
                  const double rcut2 = rcut*rcut;
                  const double d2 = norm2(dr);

                  const auto w = optional.nbh_data.get(cell_a,p_a,p_nbh_index,nbh_data_ctx);
                  if( d2 <= rcut2 && w )
                  {
                    tab[pair_ab_id].ta = type_a;
                    tab[pair_ab_id].tb = type_b;
                    tab[pair_ab_id].process_neighbor(tab[pair_ab_id], dr, d2, cells, cell_b, p_b, w );
                  }
                  /*
                  const auto w = optional.nbh_data.get(cell_a,p_a,p_nbh_index,nbh_data_ctx);
                  if( d2 <= rcut2 && w>0.0 )
                  {
                    unsigned int bi = tab[pair_ab_id].count;
                    assert( bi < tab[pair_ab_id].MaxNeighbors );
                    assert( bi==0 || (tab[pair_ab_id].ta==type_a && tab[pair_ab_id].tb==type_b) );

                    tab[pair_ab_id].ta = type_a;
                    tab[pair_ab_id].tb = type_b;

                    tab[pair_ab_id].d2[bi] = d2;
                    tab[pair_ab_id].drx[bi] = dr.x;
                    tab[pair_ab_id].dry[bi] = dr.y;
                    tab[pair_ab_id].drz[bi] = dr.z;

                    tab[pair_ab_id].nbh.set( bi , cell_b, p_b );
                    tab[pair_ab_id].nbh_data.set( bi , w );
                    tab[pair_ab_id].process_neighbor( tab[pair_ab_id], bi, cells, cell_b, p_b );
                    ++ tab[pair_ab_id].count;
                  }
                  */
                } // end of loop on sub atoms of B
                ++ p_nbh_index;
              }
            }
          }
        } // end of loop on cell groups (end of traversal for p_a neighbors)

        // call compute function of different type pairs
        for(unsigned int pair_i=0;pair_i<n_reduced_pairs;pair_i++)
        {
          if( tab[pair_i].count > 0 )
          {
            (*(funcs[pair_i])) ( tab[pair_i], molecule_atom_values[cur_sub_atom_a][onika::soatl::FieldId<field_ids>()] ... );
            tab[pair_i].count = 0;
          }
        }

      } // end of loop for sub atoms A

      // end of loop on sub atoms of particle A
      molecule_func( cell_a_pointers[onika::soatl::FieldId<field_ids>()][p_a] ...
                   , orient_a[p_a] , couple_a[p_a]
                   , nb_sub_atom_a, species[mol_type_a].m_rigid_atoms, molecule_atom_values );

    } // end of loop on cell_a particles

  } // end of compute_pair_rigidmol_cell


  template< typename GridT,
            typename MoleculeFuncT,
            typename OptionalArgsT,
            typename ComputePairBufferFactoryT,
            typename... field_ids>
  static inline void compute_pair_rigidmol(
    GridT& grid,
    const std::vector< double >& rcuts,
    const std::vector< std::shared_ptr< ComputePairOperator< FieldSet<field_ids...> > > >& funcs,
    const unsigned int NTYPES,
    const std::vector<int>& pair_id_map,
    const ParticleSpecies& species,
    bool enable_ghosts,
    MoleculeFuncT molecule_func,
    OptionalArgsT optional_in, // locks are needed if symmetric computation is enabled
    ComputePairBufferFactoryT cpbuf_factory,
    FieldSet< field_ids... >
    )
  {
    using ComputeBuffer = typename ComputePairBufferFactoryT::ComputePairBuffer;

    //const unsigned int N_PAIRS = unique_pair_count(NTYPES);
//    const double rcut2 = rcut * rcut;
    const IJK grid_dims = grid.dimension();
    int gl = grid.ghost_layers();
    if( enable_ghosts ) { gl = 0; }
    const IJK block_dims = grid_dims - (2*gl);

    assert( rcuts.size() == funcs.size() );
    assert( !pair_id_map.empty() || rcuts.size() == unique_pair_count(NTYPES) );
    assert( pair_id_map.empty() || pair_id_map.size() == unique_pair_count(NTYPES) );

    const int * __restrict__ idmapptr = pair_id_map.data();
    unsigned int nb_reduced_pairs = rcuts.size();

#   pragma omp parallel
    {
      // create per tread local objects
      //auto nbh_it = optional.nbh;
      OptionalArgsT optional = optional_in;
      ComputeBuffer tab[ nb_reduced_pairs ];

      for(unsigned int i=0;i<nb_reduced_pairs;i++)
      {
        cpbuf_factory.init(tab[i]);
        tab[i].ta = particle_id_codec::MAX_PARTICLE_TYPE;
        tab[i].tb = particle_id_codec::MAX_PARTICLE_TYPE;
      }

      GRID_OMP_FOR_BEGIN(block_dims,_,block_cell_a_loc, schedule(dynamic) )
      {
        IJK cell_a_loc = block_cell_a_loc + gl;
        size_t cell_a = grid_ijk_to_index( grid_dims , cell_a_loc );
        compute_pair_rigidmol_cell(grid,tab,cell_a,rcuts,funcs,NTYPES,idmapptr,species,molecule_func,optional, FieldSet< field_ids... >() );
      }
      GRID_OMP_FOR_END

    }

  }

}
