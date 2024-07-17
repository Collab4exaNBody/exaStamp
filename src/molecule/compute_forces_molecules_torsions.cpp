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
#include <exanb/core/thread.h>  // GridParticleLocks

#include <exaStamp/molecule/torsions_potentials_parameters.h>
#include <exaStamp/molecule/periodic_r_delta.h>

namespace exaStamp
{

  ///**
  // * Signum function : give the sign of an expression (+1, -1 or 0)
  // */
  //template <typename T>
  //typename std::enable_if<std::is_unsigned<T>::value, int>::type
  //inline constexpr signum(T x) noexcept {
  //  return T(0) < x;
  //}

  //template <typename T>
  //typename std::enable_if<std::is_signed<T>::value, int>::type
  //inline constexpr signum(T x) noexcept {
  //  return (T(0) < x) - (x < T(0));
  //}


  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_fx, field::_fy, field::_fz, field::_ep >
    >
  class ComputeForcesTorsionsNode : public OperatorNode
  {
    // -----------------------------------------------
    // Operator slots
    // -----------------------------------------------
    ADD_SLOT( GridT                       , grid                    , INPUT_OUTPUT );
    ADD_SLOT( Domain                      , domain                  , INPUT );
    ADD_SLOT( ChemicalTorsions            , chemical_torsions       , INPUT, OPTIONAL );
    ADD_SLOT( TorsionsPotentialParameters , potentials_for_torsions , INPUT, REQUIRED );
    ADD_SLOT( ParticleSpecies             , species                 , INPUT, REQUIRED );
    
    ADD_SLOT( GridParticleLocks        , particle_locks        , INPUT_OUTPUT);

    inline void execute ()  override final
    {
      static constexpr bool has_virial_field = GridHasField<GridT,field::_virial>::value;
      // using CellLockT = decltype( (*particle_locks)[0][0] );

      if( ! grid.has_value() || grid->number_of_cells()==0 )
      {
        return;
      }

      if( ! chemical_torsions.has_value() )
      {
        lerr << "chemical_torsions input missing" << std::endl;
        std::abort();
      }
      
      GridT& grid = *(this->grid);

      ChemicalTorsions&     torsions_list         = *chemical_torsions;
      ParticleSpecies& species = *(this->species);

      Mat3d xform = domain->xform();
      bool xform_is_identity = domain->xform_is_identity();

      if( potentials_for_torsions->m_type_to_potential.empty() )
      {
        // Convert name (string) to type (uint_8)
        std::map<std::string, uint8_t> map_species_name_id;
        for(size_t i=0;i<species.size(); ++i)
        {
          map_species_name_id[species.at(i).m_name] = i;
        }

        // Creation map torsions potentials
        //std::map<std::array<unsigned int, 4>, std::shared_ptr<IntraMolecularPotentialFunctional> > map_torsions_potential;
        for(const auto& b_type : potentials_for_torsions->m_potentials )
        {
          uint64_t t0 = map_species_name_id.at(b_type.species.at(0));
          uint64_t t1 = map_species_name_id.at(b_type.species.at(1));
          uint64_t t2 = map_species_name_id.at(b_type.species.at(2));
          uint64_t t3 = map_species_name_id.at(b_type.species.at(3));
          uint64_t key = (t0<<48) | (t1<<32) | (t2<<16) | t3;
          uint64_t alt_key = (t3<<48) | (t2<<32) | (t1<<16) | t0;
          if( alt_key < key ) key = alt_key;
          potentials_for_torsions->m_type_to_potential[key] = b_type.m_potential_function;

          const auto gp = b_type.m_potential_function->generic_parameters();
          ldbg << "Torsion potential for "<<b_type.species.at(0)<<"/"<<b_type.species.at(1)<<"/"<<b_type.species.at(2)<<"/"<<b_type.species.at(3)<< " : p0="<<gp.p0<<", p1="<<gp.p1<<", p2="<<gp.p2<<", x0="<<gp.x0<<", coeff="<<gp.coeff<<std::endl;
        }
      }
      auto cells = grid.cells();

      Vec3d size_box {std::abs(domain->extent().x - domain->origin().x),
                      std::abs(domain->extent().y - domain->origin().y),
                      std::abs(domain->extent().z - domain->origin().z)};
      double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 

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

      size_t n_torsion = torsions_list.size();

#     pragma omp parallel 
      {
        size_t cell[4];
        size_t pos[4];
        unsigned int type[4];
        
#       pragma omp for schedule(static)
        for(size_t i=0; i<n_torsion; ++i)
        {
          //-----------------------------DECODE----------------------------------------------------
          // const int id1 = torsions_list[i][0];
          // const int id2 = torsions_list[i][1];
          // const int id3 = torsions_list[i][2];
          // const int id4 = torsions_list[i][3];

          // At least one atom have to belong to the simulation box
          uint64_t atom_to_decode_a = torsions_list[i][0]; // atom_from_idmap(id1,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_b = torsions_list[i][1]; // atom_from_idmap(id2,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_c = torsions_list[i][2]; // atom_from_idmap(id3,id_map,id_map_ghosts); 
          uint64_t atom_to_decode_d = torsions_list[i][3]; // atom_from_idmap(id4,id_map,id_map_ghosts); 

          decode_cell_particle(atom_to_decode_a, cell[0], pos[0], type[0]);
          decode_cell_particle(atom_to_decode_b, cell[1], pos[1], type[1]);
          decode_cell_particle(atom_to_decode_c, cell[2], pos[2], type[2]);
          decode_cell_particle(atom_to_decode_d, cell[3], pos[3], type[3]);
          
          const Vec3d r0 = Vec3d { cells[cell[0]][field::rx][pos[0]] , cells[cell[0]][field::ry][pos[0]] , cells[cell[0]][field::rz][pos[0]] };
          const Vec3d r1 = Vec3d { cells[cell[1]][field::rx][pos[1]] , cells[cell[1]][field::ry][pos[1]] , cells[cell[1]][field::rz][pos[1]] };
          const Vec3d r2 = Vec3d { cells[cell[2]][field::rx][pos[2]] , cells[cell[2]][field::ry][pos[2]] , cells[cell[2]][field::rz][pos[2]] };
          const Vec3d r3 = Vec3d { cells[cell[3]][field::rx][pos[3]] , cells[cell[3]][field::ry][pos[3]] , cells[cell[3]][field::rz][pos[3]] };
          
          Vec3d r10 = periodic_r_delta( r1 , r0 , size_box , half_min_size_box );
          Vec3d r12 = periodic_r_delta( r1 , r2 , size_box , half_min_size_box );
          Vec3d r23 = periodic_r_delta( r2 , r3 , size_box , half_min_size_box );

          if( ! xform_is_identity )
          {
            r10 = xform * r10;
            r12 = xform * r12;
            r23 = xform * r23;
          }
          //---------------------------------------------------------------------------------------

          //-----------------------------READ POTENTIAL--------------------------------------------
          uint64_t t0 = type[0];
          uint64_t t1 = type[1];
          uint64_t t2 = type[2];
          uint64_t t3 = type[3];
          uint64_t key = (t0<<48) | (t1<<32) | (t2<<16) | t3;
          uint64_t alt_key = (t3<<48) | (t2<<32) | (t1<<16) | t0;
          if( alt_key < key ) key = alt_key;
          auto it = potentials_for_torsions->m_type_to_potential.find(key);
          assert( it != potentials_for_torsions->m_type_to_potential.end() );
          //---------------------------------------------------------------------------------------

          //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
          // Compute energy
          const Vec3d V = cross(r10,r12);
          const Vec3d W = cross(r23,r12);

          assert(V!=(Vec3d{0,0,0}));//case r10 and r12 aligned
          assert(W!=(Vec3d{0,0,0}));//case r12 and r23 aligned

          const double phi = signum(dot(cross(V,W), r12)) * angle(V, W);
          const auto force_energy = it->second->force_energy( phi );
          const double e = force_energy.second / 4.; //t.f_energy(phi);
          //---------------------------------------------------------------------------------------------

          const double norm_r2 = norm(r12);
          assert(norm_r2>0);

          //pa = V/norm(V)
          const double sqr_norm_V = dot(V,V);
          const Vec3d pa_r1sintheta1 = V * norm_r2 / sqr_norm_V;

          //pd = m_W/norm(m_W)
          const double sqr_norm_W = dot(W,W);
          const Vec3d pd_r3sintheta2 = - W * norm_r2 / sqr_norm_W ;

          //---------------------------------------------------------------------------------------------
          // Compute forces
          // F = - grad(Ep)
          const double dep_on_dphi = force_energy.first; //t.f_forces(phi);

          //Force on A
          const Vec3d Fa = -dep_on_dphi * pa_r1sintheta1;
          
          //Force on D
          const Vec3d Fd = -dep_on_dphi * pd_r3sintheta2;

          //Force on C
          //-1/norm(r12)^2 ((r23 + r12) ^ FD + r1 ^ FA) ^ r12    
          const double sqr_norm_r2 = dot(r12,r12);
          const Vec3d  Fc = - cross( cross(r23+r12,Fd) + cross(r10,Fa) , r12 ) / sqr_norm_r2;

          //e /= 4.;
          const Vec3d Fo = - (Fa+Fc+Fd);

          const Mat3d vira = tensor (Fa,r10) * 0.5;
          const Mat3d virc = tensor (Fc,r12) * 0.5;
          const Mat3d vird = tensor (Fd,r23) * 0.5;
          
          const Mat3d virac = vira + virc;
          const Mat3d vircd = virc + vird;

#if 0 
          // Debug Th. C.
#         pragma omp critical(dbgèmesg)
          {
            const uint64_t atid1 = cells[cell[0]][field::id][pos[0]];
            const uint64_t atid2 = cells[cell[1]][field::id][pos[1]];
            const uint64_t atid3 = cells[cell[2]][field::id][pos[2]];
            const uint64_t atid4 = cells[cell[3]][field::id][pos[3]];
            if( atid1==846 || atid1==3004 || atid1==11323 || atid4==846 || atid4==3004 || atid4==11323 )
            {
              printf("TORSION %05ld - %05ld - %05ld - %05ld : x1=(% .5e,% .5e,% .5e) x2=(% .5e,% .5e,% .5e) x3=(% .5e,% .5e,% .5e) x4=(% .5e,% .5e,% .5e) f=(% .5e,% .5e,% .5e)\n"
                    , atid1 , atid2 , atid3 , atid4
                    , r0.x , r0.y , r0.z
                    , r1.x , r1.y , r1.z
                    , r2.x , r2.y , r2.z
                    , r3.x , r3.y , r3.z
                    , Fa.x , Fa.y , Fa.z );
            }
          }
          /*******************/
#endif

//#         pragma omp critical
//          {

          (*particle_locks)[cell[0]][pos[0]].lock();
          cells[cell[0]][field::ep][pos[0]] += e;
          cells[cell[0]][field::fx][pos[0]] += Fa.x; 
          cells[cell[0]][field::fy][pos[0]] += Fa.y;
          cells[cell[0]][field::fz][pos[0]] += Fa.z; 
          if constexpr (has_virial_field) cells[cell[0]][field::virial][pos[0]] += vira;
          (*particle_locks)[cell[0]][pos[0]].unlock();

          (*particle_locks)[cell[1]][pos[1]].lock();
          cells[cell[1]][field::ep][pos[1]] += e;
          cells[cell[1]][field::fx][pos[1]] += Fo.x; 
          cells[cell[1]][field::fy][pos[1]] += Fo.y;
          cells[cell[1]][field::fz][pos[1]] += Fo.z; 
          if constexpr (has_virial_field) cells[cell[1]][field::virial][pos[1]] += virac;
          (*particle_locks)[cell[1]][pos[1]].unlock();


          (*particle_locks)[cell[2]][pos[2]].lock();
          cells[cell[2]][field::ep][pos[2]] += e;
          cells[cell[2]][field::fx][pos[2]] += Fc.x; 
          cells[cell[2]][field::fy][pos[2]] += Fc.y;
          cells[cell[2]][field::fz][pos[2]] += Fc.z; 
          if constexpr (has_virial_field) cells[cell[2]][field::virial][pos[2]] += vircd;
          (*particle_locks)[cell[2]][pos[2]].unlock();

          (*particle_locks)[cell[3]][pos[3]].lock();
          cells[cell[3]][field::ep][pos[3]] += e;
          cells[cell[3]][field::fx][pos[3]] += Fd.x; 
          cells[cell[3]][field::fy][pos[3]] += Fd.y;
          cells[cell[3]][field::fz][pos[3]] += Fd.z; 
          if constexpr (has_virial_field) cells[cell[3]][field::virial][pos[3]] += vird;
          (*particle_locks)[cell[3]][pos[3]].unlock();
//          }
        }
      }
    }

  };


  template<class GridT> using ComputeForcesTorsionsNodeTmpl = ComputeForcesTorsionsNode<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "compute_force_torsion", make_grid_variant_operator< ComputeForcesTorsionsNodeTmpl> );
  }

}


