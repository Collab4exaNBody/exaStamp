#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/geometry.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/log.h>
#include <exanb/core/thread.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/std_array_hash.h>

#include <memory>
#include <vector>
#include <atomic>

#include <iostream>

namespace exaStamp
{
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_id, field::_idmol, field::_cmol>
    >
  class ReplicateDomainForMoleculesNode : public OperatorNode
  {

    // compile time constant indicating if grid has id field
    using has_id_field_t = typename GridT::CellParticles::template HasField < field::_id > ;
    static constexpr bool has_id_field = has_id_field_t::value;

    // compile time constant indicating if grid has id field
    using has_id_mol_field_t = typename GridT::CellParticles::template HasField < field::_idmol > ;
    static constexpr bool has_id_mol_field = has_id_mol_field_t::value;

    ADD_SLOT( GridT             , grid               , INPUT_OUTPUT );
    ADD_SLOT( Domain            , domain             , INPUT_OUTPUT );
    ADD_SLOT( IJK               , repeat             , INPUT );
    ADD_SLOT( ChemicalBonds     , chemical_bonds     , INPUT_OUTPUT );
    ADD_SLOT( ChemicalAngles    , chemical_angles    , INPUT_OUTPUT );
    ADD_SLOT( ChemicalTorsions  , chemical_torsions  , INPUT_OUTPUT );
    ADD_SLOT( ChemicalImpropers , chemical_impropers , INPUT_OUTPUT );

  public:
    inline void execute () override final
    {
      GridT& grid = *(this->grid);
      Domain& domain = *(this->domain);
      IJK repeat = *(this->repeat);
      ChemicalBonds&     bonds     = *chemical_bonds;
      ChemicalAngles&    angles    = *chemical_angles;
      ChemicalTorsions&  torsions  = *chemical_torsions;
      ChemicalImpropers& impropers = *chemical_impropers;

      if(!has_id_field)
        {
          lerr << "Error in replicate domaine for molecules : no field id in grid flavor." << std::endl;
          abort();
        }

      if(!has_id_mol_field)
        {
          lerr << "Error in replicate domaine for molecules : no field id_mol in grid flavor." << std::endl;
          abort();
        }

      Vec3d dom_size = bounds_size(domain.bounds());
      IJK dom_dims = domain.grid_dimension();
      domain.set_bounds( { domain.bounds().bmin , domain.bounds().bmin + dom_size * repeat } );
      domain.set_grid_dimension( domain.grid_dimension() * repeat );

      ssize_t ghost_layers = grid.ghost_layers();
      IJK grid_dims_with_ghosts = grid.dimension();
      IJK grid_dims = grid.dimension() - (2*ghost_layers);

      GridT ngrid;
      ngrid.set_offset( grid.offset() + ghost_layers );
      ngrid.set_origin( grid.origin() );
      ngrid.set_cell_size( grid.cell_size() );
      IJK ngrid_dims = grid_dims + dom_dims * (repeat-1);
      ngrid.set_dimension( ngrid_dims );
      ngrid.set_max_neighbor_distance( 0.0 );

      auto src_cells = grid.cells();
      auto dst_cells = ngrid.cells();

      int bonds_size     = bonds.size();
      int angles_size    = angles.size();
      int torsions_size  = torsions.size();
      int impropers_size = impropers.size();

      //set an id of molecules
      //we follow the connectivity list
      uint64_t max_id=std::numeric_limits<uint64_t>::lowest();
      uint64_t max_molecules_id=std::numeric_limits<uint64_t>::lowest();
      size_t n_cells = grid.number_of_cells();
      for(size_t cell_i=0;cell_i<n_cells;cell_i++)
        {
          for(size_t i=0;i<src_cells[cell_i].size();i++)
            {
              if(src_cells[cell_i][field::id][i]>max_id)
                {
                  max_id=src_cells[cell_i][field::id][i];
                }

              if(src_cells[cell_i][field::idmol][i]>max_molecules_id)
                {
                  max_molecules_id=src_cells[cell_i][field::idmol][i];
                }
            }
        }


#     ifndef NDEBUG
      size_t src_ncells = grid.number_of_cells();
      size_t dst_ncells = ngrid.number_of_cells();
#     endif

      //size_t local_particles = 0;
      for(ssize_t rk=0;rk<repeat.k;rk++)
      {
        for(ssize_t rj=0;rj<repeat.j;rj++)
        {
          for(ssize_t ri=0;ri<repeat.i;ri++)
          {
            ssize_t index = grid_ijk_to_index(repeat, IJK{ri, rj, rk});

            IJK grid_displ = dom_dims * IJK{ri,rj,rk};
            Vec3d r_shift = dom_size * IJK{ri,rj,rk};
            GRID_FOR_BEGIN(grid_dims,_,cell_loc)
            {
              IJK src_loc = cell_loc + ghost_layers;
              assert( grid.contains( src_loc ) );

              IJK dst_loc = cell_loc + grid_displ;
              assert( ngrid.contains( dst_loc ) );

              size_t src_cell_i = grid_ijk_to_index( grid_dims_with_ghosts , src_loc );
              assert( src_cell_i>=0 && src_cell_i<src_ncells );

              size_t dst_cell_i = grid_ijk_to_index( ngrid_dims , dst_loc );
              assert( dst_cell_i>=0 && dst_cell_i<dst_ncells );

              assert( dst_cells[dst_cell_i].size() == 0 );
              size_t n = src_cells[src_cell_i].size();
              dst_cells[dst_cell_i].resize( n );

              for(size_t p_i=0;p_i<n;p_i++)
              {
                auto p = src_cells[src_cell_i][p_i];
                p[ field::rx ] += r_shift.x;
                p[ field::ry ] += r_shift.y;
                p[ field::rz ] += r_shift.z;

                //std::cout << p[ field::vx ] << " " << p[ field::vz ] << " " << p[ field::vz ] << " " << std::endl;

                p[ field::id ] += max_id*index;

                p[ field::idmol ] += max_molecules_id*index;
                for(size_t ci=0; ci<p[ field::cmol ].size();ci++)
                  if(p[ field::cmol ][ci] != std::numeric_limits<uint64_t>::max())
                    p[ field::cmol ][ci] += max_id*index;

                dst_cells[dst_cell_i].write_tuple( p_i, p );
              }
              GRID_FOR_END

              //bonds transfer
              for(int b=0; b<bonds_size;++b)
                {
                  std::array<uint64_t, 2> bond_to_add;
                  for(int pos=0; pos<2;++pos)
                    {
                      bond_to_add[pos] = bonds[b][pos]+max_id*index;
                    }
                  if(index!=0)
                    bonds.push_back(bond_to_add);
                }

              //angles transfer
              for(int a=0; a<angles_size;++a)
                {
                  std::array<uint64_t, 3> angle_to_add;
                  for(int pos=0; pos<3;++pos)
                    {
                      angle_to_add[pos] = angles[a][pos]+max_id*index;
                    }
                  if(index!=0)
                    angles.push_back(angle_to_add);
                }

              //torsions transfer
              for(int t=0; t<torsions_size;++t)
                {
                  std::array<uint64_t, 4> torsion_to_add;
                  for(int pos=0; pos<4;++pos)
                    {
                      torsion_to_add[pos] = torsions[t][pos]+max_id*index;
                    }
                  if(index!=0)
                    torsions.push_back(torsion_to_add);
                }

              //impropers transfer
              for(int i=0; i<impropers_size;++i)
                {
                  std::array<uint64_t, 4> improper_to_add;
                  for(int pos=0; pos<4;++pos)
                    {
                      improper_to_add[pos] = impropers[i][pos]+max_id*index;
                    }
                  if(index!=0)
                    impropers.push_back(improper_to_add);
                }
            }
          }
        }
      }

      ngrid.rebuild_particle_offsets();
      grid = std::move( ngrid );
    }

  };

  template<class GridT> using ReplicateDomainForMoleculesNodeTmpl = ReplicateDomainForMoleculesNode<GridT>;

   // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "replicate_domain_for_molecules", make_grid_variant_operator< ReplicateDomainForMoleculesNodeTmpl > );
  }

}
