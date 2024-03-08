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
#include <exaStamp/compute/virial_add_contribution.h>


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
    ADD_SLOT( GridT                   , grid                  , INPUT_OUTPUT );
    ADD_SLOT( Domain                  , domain                , INPUT );
    ADD_SLOT( ChemicalBonds           , chemical_bonds        , INPUT, OPTIONAL );
    ADD_SLOT( GridAtomsToChemicalChain, atoms_to_bonds        , INPUT, OPTIONAL );
    ADD_SLOT( BondComputeBuffer       , chemical_bonds_force_buffer , INPUT_OUTPUT );
    ADD_SLOT( BondsPotentialParameters, potentials_for_bonds  , INPUT_OUTPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies         , species               , INPUT, REQUIRED );

  public:
    inline void execute() override final
    {    
      // compile time constant indicating if grid has virial field
      static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

      if( ! grid.has_value() || grid->number_of_cells()==0 )
      {
        return;
      }

      if( ! chemical_bonds.has_value() )
      {
        lerr << "chemical_bonds input missing" << std::endl;
        std::abort();
      }
      
      if( ! atoms_to_bonds.has_value() )
      {
        lerr << "atoms_to_bonds input missing" << std::endl;
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
          lerr<<"xform="<<xform<<", size_box="<<size_box<<", tmp="<<tmp<<std::endl;
          std::abort();
        }
      }
#     endif

      size_t n_bonds = bonds_list.size();
      
      auto& force_buf = *chemical_bonds_force_buffer;
      force_buf.m_force_enegery.resize( n_bonds );
      if( has_virial_field )
      {
        force_buf.m_virial.resize( n_bonds );
      }
      
      // std::cout<<"n_bonds = "<<n_bonds<<std::endl;
      
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

          force_buf.m_force_enegery[i].m_force = F;
          force_buf.m_force_enegery[i].m_energy = e * 0.5;
          if( has_virial_field )
          {
            force_buf.m_virial[i] = tensor(F,r) * 0.5;
          }

        }
      }

      size_t n_cells = grid->number_of_cells();

#     pragma omp parallel for schedule(dynamic)
      for(size_t i=0;i<n_cells;i++)
      {
        size_t n_particles = cells[i].size();
        auto* __restrict__ cell_ep = cells[i][field::ep];
        auto* __restrict__ cell_ax = cells[i][field::ax];
        auto* __restrict__ cell_ay = cells[i][field::ay];
        auto* __restrict__ cell_az = cells[i][field::az];
        auto& cell_ab = (*atoms_to_bonds)[i];
        size_t k = 0;
        for(size_t j=0;j<n_particles;j++)
        {

#         ifndef NDEBUG
          if( cell_ab[k] != j ) { lerr<<"Bad particle index : found "<<cell_ab[k]<<", expected "<<j<<std::endl; std::abort(); }
          ++k;
#         endif

          double ep = 0.0;
          Vec3d F {0.,0.,0.};
          Mat3d virial; // initialized to all 0
          int count = cell_ab[k++];
          for(int b=0;b<count;b++)
          {
            uint64_t eb = cell_ab[k++];
            size_t bi = eb / 4;
            double fs = 1.0 - ( ( eb % 4 )*2 ) ; // force sign
            assert( fs==-1.0 || fs==1.0 );
            ep += force_buf.m_force_enegery[bi].m_energy ;
            F += force_buf.m_force_enegery[bi].m_force * fs;
            if( has_virial_field )
            {
              virial += force_buf.m_virial[bi];
            }
          }
          cell_ep[j] += ep;
          cell_ax[j] += F.x;
          cell_ay[j] += F.y;
          cell_az[j] += F.z;
          virial_add_contribution(cells[i],j,virial);
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

