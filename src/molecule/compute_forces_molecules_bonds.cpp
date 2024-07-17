#include <yaml-cpp/yaml.h>
#include <memory>
//#include <utility>// std::pair
#include <iomanip>

#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>

#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/basic_types_yaml.h>
#include <exaStamp/molecule/mol_connectivity.h>
#include <exanb/core/particle_id_codec.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/log.h>

#include <exaStamp/molecule/bonds_potentials_parameters.h>
#include <exaStamp/molecule/periodic_r_delta.h>
#include <exanb/core/concurent_add_contributions.h>

#include <string>
#include <iostream>
#include <sstream>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class ComputeForcesBondsNode : public OperatorNode
  {

    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( Domain                  , domain                , INPUT );
    ADD_SLOT( ChemicalBonds           , chemical_bonds        , INPUT, OPTIONAL );
    ADD_SLOT( ParticleSpecies         , species               , INPUT, REQUIRED );
    ADD_SLOT( BondsPotentialParameters, potentials_for_bonds  , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( GridT                   , grid                  , INPUT_OUTPUT );
    ADD_SLOT( GridParticleLocks       , particle_locks        , INPUT_OUTPUT);

  public:
    inline void execute() override final
    {    
      // compile time constant indicating if grid has virial field
      static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;
      using CellLockT = decltype( (*particle_locks)[0][0] );

      if( ! grid.has_value() || grid->number_of_cells()==0 )
      {
        return;
      }

      if( ! chemical_bonds.has_value() )
      {
        lerr << "chemical_bonds input missing" << std::endl;
        std::abort();
      }
      
      ChemicalBonds&            bonds_list    = *chemical_bonds;

      Mat3d xform = domain->xform();
      bool xform_is_identity = domain->xform_is_identity();

      // Re build the map containing the bond potential parameters, in case it was not done before
      // by a call to read_bond_potentials
      if( potentials_for_bonds->m_type_pair_to_potential.empty() )
      {
        // Convert atom type name (string) to type (unsigned int)
        std::unordered_map<std::string, unsigned int> map_species_name_id;
        for(size_t i=0;i<species->size(); ++i)
        {
          map_species_name_id[species->at(i).m_name] = i;
        }

        // pre-build map of bond type pair to bond potential function
        for( const auto& b_type : potentials_for_bonds->m_bond_desc )
        {
          uint64_t a = map_species_name_id.at(b_type.species.at(0));
          uint64_t b = map_species_name_id.at(b_type.species.at(1));
          if( a > b ) std::swap( a , b );
          assert( a < 65536 );
          assert( b < 65536 );
          potentials_for_bonds->m_type_pair_to_potential[ (a<<16) | b ] = b_type.potential;
          const auto gp = b_type.potential->generic_parameters();
          ldbg << "Bond potential for "<<b_type.species.at(0)<<" / "<<b_type.species.at(1)<<" : p0="<<gp.p0<<", p1="<<gp.p1<<", p2="<<gp.p2<<", x0="<<gp.x0<<", coeff="<<gp.coeff<<std::endl;
        }
      }

      const auto& map_bonds_potential = potentials_for_bonds->m_type_pair_to_potential;

      auto cells = grid->cells();

      const Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 

#     ifndef NDEBUG
      if( ! xform_is_identity )
      {
        Vec3d tmp = xform * size_box;
        if( fabs(tmp.x) <= 1.5 || fabs(tmp.y) <= 1.5 || fabs(tmp.z) <= 1.5 )
        {
          fatal_error()<<"xform="<<xform<<", size_box="<<size_box<<", tmp="<<tmp<<std::endl;
        }
      }
#     endif

      size_t n_bonds = bonds_list.size();
      ldbg<<"n_bonds = "<<n_bonds<<std::endl;
      
#     pragma omp parallel //for shared(map_bonds_potential, cells) private(cell,pos,type)
      {
        size_t       cell[2];
        size_t       pos[2];
        unsigned int type[2];

#       pragma omp for
        for(size_t i=0; i<n_bonds; ++i)
        {
          // ---- decode using local ids in bonds -----
          const uint64_t atom_to_decode_a = bonds_list[i][0]; //atom_from_idmap( bonds_list[i][0] ,id_map, id_map_ghosts ); 
          const uint64_t atom_to_decode_b = bonds_list[i][1]; //atom_from_idmap( bonds_list[i][1] ,id_map, id_map_ghosts ); 
          
          decode_cell_particle(atom_to_decode_a, cell[0], pos[0], type[0]);
          decode_cell_particle(atom_to_decode_b, cell[1], pos[1], type[1]);
          const Vec3d ra = Vec3d { cells[cell[0]][field::rx][pos[0]] , cells[cell[0]][field::ry][pos[0]] , cells[cell[0]][field::rz][pos[0]] };
          const Vec3d rb = Vec3d { cells[cell[1]][field::rx][pos[1]] , cells[cell[1]][field::ry][pos[1]] , cells[cell[1]][field::rz][pos[1]] };
          Vec3d r = periodic_r_delta( ra , rb , size_box , half_min_size_box );

          if( ! xform_is_identity ) { r = xform * r; }
          // --------------------------------------------

          //-----------------------------READ POTENTIAL--------------------------------------------
          // Get potential type from species of atoms pair
          uint64_t type_a = type[0];
          uint64_t type_b = type[1];
          if( type_a > type_b ) std::swap( type_a , type_b );          
          auto it = map_bonds_potential.find( ( type_a << 16 ) | type_b );
          assert( it != map_bonds_potential.end() );
          //---------------------------------------------------------------------------------------

          //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
          double norm_r = norm(r);
          assert(norm_r>0);

          // compute the pair dEp/dr , E
          const auto fe = it->second->force_energy( norm_r );

          // Compute energy
          double e = fe.second;

          // Compute forces
          // F = dEp/dr nij
          double dEp_dr = fe.first;
          Vec3d F = ( r * dEp_dr ) / norm_r;

#if 0
          // Debug Th. C.
#         pragma omp critical(dbg_mesg)
          {
            const uint64_t atid1 = cells[cell[0]][field::id][pos[0]];
            const uint64_t atid2 = cells[cell[1]][field::id][pos[1]];
            if( atid1==846 || atid1==3004 || atid1==11323 )
            {
              printf("BOND %05ld - %05ld : x1=(% .5e,% .5e,% .5e) x2=(% .5e,% .5e,% .5e) dr=(% .5e,% .5e,% .5e) f=(% .5e,% .5e,% .5e)\n"
                    , atid1 , atid2
                    , ra.x , ra.y , ra.z
                    , rb.x , rb.y , rb.z
                    , r.x , r.y , r.z
                    , F.x , F.y , F.z );
            }
          }
          /*******************/
#endif
          
/*
          force_buf.m_force_enegery[i].m_force = F;
          force_buf.m_force_enegery[i].m_energy = e * 0.5;
          if( has_virial_field )
          {
            force_buf.m_virial[i] = tensor(F,r) * 0.5;
          }
*/
          if constexpr ( has_virial_field )
          {
            auto virial = tensor(F,r) * 0.5;
            concurent_add_contributions<CellLockT,false,true,double,double,double,double,Mat3d> (
                (*particle_locks)[cell[0]][pos[0]]
              , cells[cell[0]][field::fx][pos[0]], cells[cell[0]][field::fy][pos[0]], cells[cell[0]][field::fz][pos[0]], cells[cell[0]][field::ep][pos[0]], cells[cell[0]][field::virial][pos[0]]
              , F.x, F.y, F.z, e*0.5, virial );

            concurent_add_contributions<CellLockT,false,true,double,double,double,double,Mat3d> (
                (*particle_locks)[cell[1]][pos[1]]
              , cells[cell[1]][field::fx][pos[1]], cells[cell[1]][field::fy][pos[1]], cells[cell[1]][field::fz][pos[1]], cells[cell[1]][field::ep][pos[1]], cells[cell[1]][field::virial][pos[1]]
              , -F.x, -F.y, -F.z, e*0.5, virial );
          }
          else
          {
            concurent_add_contributions<CellLockT,false,true,double,double,double,double> (
                (*particle_locks)[cell[0]][pos[0]]
              , cells[cell[0]][field::fx][pos[0]], cells[cell[0]][field::fy][pos[0]], cells[cell[0]][field::fz][pos[0]], cells[cell[0]][field::ep][pos[0]]
              , F.x, F.y, F.z, e*0.5 );

            concurent_add_contributions<CellLockT,false,true,double,double,double,double> (
                (*particle_locks)[cell[1]][pos[1]]
              , cells[cell[1]][field::fx][pos[1]], cells[cell[1]][field::fy][pos[1]], cells[cell[1]][field::fz][pos[1]], cells[cell[1]][field::ep][pos[1]]
              , -F.x, -F.y, -F.z, e*0.5 );
          }
        }
      }

    }

  };

  template<class GridT> using ComputeForcesBondsNodeTmpl = ComputeForcesBondsNode<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "compute_force_bond", make_grid_variant_operator< ComputeForcesBondsNodeTmpl > );
  }

}

