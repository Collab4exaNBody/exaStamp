#include <yaml-cpp/yaml.h>
#include <memory>
#include <utility>// std::pair

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
#include <exanb/core/thread.h> // GridParticleLocks

#include <exaStamp/molecule/bends_potentials_parameters.h>
#include <exaStamp/molecule/periodic_r_delta.h>
#include <exaStamp/compute/virial_add_contribution.h>

namespace exaStamp
{

  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_ep, field::_ax, field::_ay, field::_az >
    >
  class ComputeForcesBendsNode : public OperatorNode
  {    
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( GridT                    , grid                  , INPUT_OUTPUT );
    ADD_SLOT( Domain                   , domain                , INPUT );
    ADD_SLOT( ChemicalAngles           , chemical_angles       , INPUT, OPTIONAL );
    ADD_SLOT( BendsPotentialParameters , potentials_for_angles , INPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies          , species               , INPUT, REQUIRED );

#   ifdef XSTAMP_USE_BEND_FORCE_BUFFER
    ADD_SLOT( GridAtomsToChemicalChain , atoms_to_angles       , INPUT, OPTIONAL );
    ADD_SLOT( BendComputeBuffer        , chemical_angles_force_buffer , INPUT_OUTPUT );
#   else
    ADD_SLOT( GridParticleLocks        , particle_locks        , INPUT_OUTPUT);
#   endif

  public:
    inline void execute ()  override final
    {
      // static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;

      if( ! grid.has_value() || grid->number_of_cells()==0 )
      {
        return;
      }

      if( ! chemical_angles.has_value() )
      {
        lerr << "chemical_angles input missing" << std::endl;
        std::abort();
      }
      
#     ifdef XSTAMP_USE_BEND_FORCE_BUFFER
      if( ! atoms_to_angles.has_value() )
      {
        lerr << "atoms_to_angles input missing" << std::endl;
        std::abort();        
      }
#     endif

      ChemicalAngles&           bends_list  = *chemical_angles;

      Mat3d xform = domain->xform();
      bool xform_is_identity = domain->xform_is_identity();

      // build type pair to bend potential map if needed
      if( potentials_for_angles->m_type_to_potential.empty() )
      {
        std::unordered_map<std::string, unsigned int> map_species_name_id;
        for(size_t i=0;i<species->size(); ++i)
        {
          map_species_name_id[species->at(i).m_name] = i;
        }        
        for(const auto& b_type : potentials_for_angles->m_potentials)
        {          
          uint64_t type_a = map_species_name_id.at(b_type.species.at(0));
          uint64_t type_b = map_species_name_id.at(b_type.species.at(1));
          uint64_t type_c = map_species_name_id.at(b_type.species.at(2));
          if( type_a > type_c ) std::swap( type_a, type_c );
          assert( type_a < 65536 );
          assert( type_b < 65536 );
          assert( type_c < 65536 );
          uint64_t key = (type_a<<32) | (type_b<<16) | type_c ;
          assert( potentials_for_angles->m_type_to_potential.find(key) == potentials_for_angles->m_type_to_potential.end() );
          potentials_for_angles->m_type_to_potential[key] = b_type.m_potential_function;
        }
      }
      const auto& map_bends_potential = potentials_for_angles->m_type_to_potential;

      auto cells = grid->cells();

      const Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 

      const size_t n_bends = bends_list.size();
      //std::cout<<"n_bends = "<<n_bends<<std::endl;

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


#     ifdef XSTAMP_USE_BEND_FORCE_BUFFER
      auto& force_buf = *chemical_angles_force_buffer;
      force_buf.m_force_enegery.resize( n_bends );
      if( has_virial_field )
      {
        force_buf.m_virial.resize( n_bends );
      }
#     endif

#     pragma omp parallel
      {
        size_t cell[3];
        size_t pos[3];
        unsigned int type[3];

#       pragma omp for schedule(dynamic)
        for(size_t i=0; i<n_bends; ++i)
        {
          //-----------------------------DECODE----------------------------------------------------
          uint64_t atom_to_decode_a = bends_list[i][0]; // atom_from_idmap(bends_list[i][0],id_map,id_map_ghosts); 
          uint64_t atom_to_decode_b = bends_list[i][1]; // atom_from_idmap(bends_list[i][1],id_map,id_map_ghosts); 
          uint64_t atom_to_decode_c = bends_list[i][2]; // atom_from_idmap(bends_list[i][2],id_map,id_map_ghosts); 
          assert(atom_to_decode_a != std::numeric_limits<uint64_t>::max());
          assert(atom_to_decode_b != std::numeric_limits<uint64_t>::max());
          assert(atom_to_decode_c != std::numeric_limits<uint64_t>::max());

          decode_cell_particle(atom_to_decode_a, cell[0], pos[0], type[0]);
          decode_cell_particle(atom_to_decode_b, cell[1], pos[1], type[1]);
          decode_cell_particle(atom_to_decode_c, cell[2], pos[2], type[2]);

          // bend atoms
          const Vec3d ra = Vec3d { cells[cell[0]][field::rx][pos[0]] , cells[cell[0]][field::ry][pos[0]] , cells[cell[0]][field::rz][pos[0]] };
          const Vec3d rb = Vec3d { cells[cell[1]][field::rx][pos[1]] , cells[cell[1]][field::ry][pos[1]] , cells[cell[1]][field::rz][pos[1]] };
          const Vec3d rc = Vec3d { cells[cell[2]][field::rx][pos[2]] , cells[cell[2]][field::ry][pos[2]] , cells[cell[2]][field::rz][pos[2]] };
          
          // first harm
          Vec3d r1 = periodic_r_delta( rb , ra , size_box , half_min_size_box );

          // second harm
          Vec3d r2 = periodic_r_delta( rb , rc , size_box , half_min_size_box ); // rc - rb;
          
          //ldbg << "r1 avant : " << r1 << std::endl;
          if( ! xform_is_identity )
          {
            r1 = xform * r1;
            r2 = xform * r2;
          }
          //ldbg << "r1 après : " << r1 << std::endl;

          //---------------------------------------------------------------------------------------


          //-----------------------------READ POTENTIAL--------------------------------------------
          // Get potential type from species of atoms pair
          uint64_t type_a = type[0];
          uint64_t type_b = type[1];
          uint64_t type_c = type[2];
          if( type_a > type_c ) std::swap( type_a, type_c );
          assert( type_a < 65536 );
          assert( type_b < 65536 );
          assert( type_c < 65536 );
          uint64_t key = (type_a<<32) | (type_b<<16) | type_c ;
          auto it = map_bends_potential.find( key );
          assert( it != map_bends_potential.end() );
          //---------------------------------------------------------------------------------------


          //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
          const double norm_r1 = norm(r1);
          const double norm_r2 = norm(r2);
          assert(norm_r1>0);
          assert(norm_r2>0);

          //---------------------------------------------------------------------------------------------
          // Compute directions of Fa and Fb
          Vec3d tmp = cross(r1, r2);
          assert( tmp.x!=0. || tmp.y!=0. || tmp.z!=0. );
          const Vec3d p_a = cross(r1 , tmp);
          const Vec3d p_b = cross(tmp, r2 );
          const double norm_pa = norm(p_a);
          const Vec3d pa_r1 = p_a / norm_pa / norm_r1;
          const double norm_pb = norm(p_b);
          const Vec3d pb_r2 = p_b / norm_pb / norm_r2;
          //---------------------------------------------------------------------------------------------

          //---------------------------------------------------------------------------------------------
          // Compute energy
          const double theta = angle(r1,r2);
          const auto fe = it->second->force_energy( theta );
          const double e = fe.second; //b.f_energy(theta);
          //---------------------------------------------------------------------------------------------

          //---------------------------------------------------------------------------------------------
          // Compute forces
          // F = - grad(Ep)
          const double dep_on_dtheta = fe.first; //b.f_forces(theta);

          const Vec3d F1 = pa_r1 * (-dep_on_dtheta);
          const Vec3d F2 = pb_r2 * (-dep_on_dtheta);

          const Mat3d virial1 = tensor (F1,r1) * 0.5;
          const Mat3d virial2 = tensor (F2,r2) * 0.5;

//          compute_bend_interaction(r1, r2, e, F1, F2, * (it->second) );
          // add force, energy and virial contributions to 3 corresponding atoms

#         ifdef XSTAMP_USE_BEND_FORCE_BUFFER
          force_buf.m_force_enegery[i].m_force = { F1 , F2 };
          force_buf.m_force_enegery[i].m_energy = e / 3.;
          if( has_virial_field ) { force_buf.m_virial[i] = { virial1 , virial2 }; }
#         else
          const Vec3d Fo = - ( F1 + F2 );
          const Mat3d virialo = virial1 + virial2;

          (*particle_locks)[cell[0]][pos[0]].lock();
          cells[cell[0]][field::ep][pos[0]] += e/3;
          cells[cell[0]][field::fx][pos[0]] += F1.x; 
          cells[cell[0]][field::fy][pos[0]] += F1.y;
          cells[cell[0]][field::fz][pos[0]] += F1.z; 
          virial_add_contribution( cells[cell[0]] , pos[0] , virial1 );
          (*particle_locks)[cell[0]][pos[0]].unlock();

          (*particle_locks)[cell[1]][pos[1]].lock();
          cells[cell[1]][field::ep][pos[1]] += e/3;
          cells[cell[1]][field::fx][pos[1]] += Fo.x;
          cells[cell[1]][field::fy][pos[1]] += Fo.y;
          cells[cell[1]][field::fz][pos[1]] += Fo.z;
          virial_add_contribution( cells[cell[1]] , pos[1] , virialo );
          (*particle_locks)[cell[1]][pos[1]].unlock();

          (*particle_locks)[cell[2]][pos[2]].lock();
          cells[cell[2]][field::ep][pos[2]] += e/3;
          cells[cell[2]][field::fx][pos[2]] += F2.x;
          cells[cell[2]][field::fy][pos[2]] += F2.y;
          cells[cell[2]][field::fz][pos[2]] += F2.z;
          virial_add_contribution( cells[cell[2]] , pos[2] , virial2 );
          (*particle_locks)[cell[2]][pos[2]].unlock();
#       endif

        }
      }

#     ifdef XSTAMP_USE_BEND_FORCE_BUFFER
      size_t n_cells = grid->number_of_cells();
#     pragma omp parallel for schedule(dynamic)
      for(size_t i=0;i<n_cells;i++)
      {
        size_t n_particles = cells[i].size();
        auto* __restrict__ cell_ep = cells[i][field::ep];
        auto* __restrict__ cell_fx = cells[i][field::fx];
        auto* __restrict__ cell_fy = cells[i][field::fy];
        auto* __restrict__ cell_fz = cells[i][field::fz];
        auto& cell_aa = (*atoms_to_angles)[i];
        size_t k = 0;
        for(size_t j=0;j<n_particles;j++)
        {

#         ifndef NDEBUG
          if( cell_aa[k] != j ) { lerr<<"Bad particle index : found "<<cell_aa[k]<<", expected "<<j<<std::endl; std::abort(); }
          ++k;
#         endif

          double ep = 0.0;
          Vec3d F {0.,0.,0.};
          Mat3d virial; // initialized to all 0
          
          int count = cell_aa[k++];
          for(int b=0;b<count;b++)
          {
            uint64_t eb = cell_aa[k++];
            size_t bi = eb / 4;
            assert( bi < n_bends );
            unsigned int p = eb % 4;
            assert( p < 3 );
            const Vec3d F1 = force_buf.m_force_enegery[bi].m_force[0];
            const Vec3d F2 = force_buf.m_force_enegery[bi].m_force[1];
            const Vec3d force[3] = { F1 , -(F1+F2) , F2 };
            ep += force_buf.m_force_enegery[bi].m_energy ;
            F += force[p];
            if( has_virial_field )
            {
              const Mat3d vir1 = force_buf.m_virial[bi][0];
              const Mat3d vir2 = force_buf.m_virial[bi][1];
              const Mat3d vir[3] = { vir1 , vir1+vir2 , vir2 };
              virial += vir[p];
            }
          }
          
          cell_ep[j] += ep;
          cell_fx[j] += F.x;
          cell_fy[j] += F.y;
          cell_fz[j] += F.z;
          virial_add_contribution(cells[i],j,virial);
        }
        assert( k == cell_aa.size() );
      }
#     endif      
     
    }

  };

  template<class GridT> using ComputeForcesBendsNodeTmpl = ComputeForcesBendsNode<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "compute_force_bend", make_grid_variant_operator< ComputeForcesBendsNodeTmpl > );
  }

}

