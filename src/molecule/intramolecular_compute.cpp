#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/particle_id_codec.h>
#include <exanb/core/domain.h>
#include <exaStamp/molecule/molecule_list.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/molecule_compute_buffer.h>
#include <exaStamp/molecule/periodic_r_delta.h>
#include <exaStamp/molecule/potential_functional.h>

namespace exaStamp
{
  template<
    class GridT,
    class = AssertGridHasFields< GridT, field::_idmol, field::_fx, field::_fy, field::_fz, field::_ep >
    >
  class IntramolecularForceCompute : public OperatorNode
  {
    using has_virial_field_t = typename GridT::CellParticles::template HasField < field::_virial > ;
    static constexpr bool has_virial_field = has_virial_field_t::value;
    using ComputeBufferT = std::conditional_t< has_virial_field , MoleculeVirialComputeBuffer, MoleculeComputeBuffer >;
  
    ADD_SLOT( GridT                    , grid                        , INPUT_OUTPUT);
    ADD_SLOT( Domain                   , domain                      , INPUT , REQUIRED );
    ADD_SLOT( MoleculeSpeciesVector    , molecules                   , INPUT , MoleculeSpeciesVector{}    , DocString{"Molecule descriptions"} );
    ADD_SLOT( MoleculeLists            , molecule_list               , INPUT , MoleculeLists{}            , DocString{"List of atoms for each molecule"} );
    ADD_SLOT( MoleculeSetComputeParams , molecule_compute_parameters , INPUT , MoleculeSetComputeParams{} , DocString{"Intramolecular functionals' parameters"} );

  public:
    inline void execute ()  override final
    {
      if( grid->number_of_cells() == 0 || molecules->m_molecules.empty() || molecule_list->number_of_molecules() == 0 ) return;
    
      auto cells = grid->cells();
      const size_t nmol = molecule_list->number_of_molecules();
      const Vec3d size_box { std::abs(domain->extent().x - domain->origin().x)
                           , std::abs(domain->extent().y - domain->origin().y)
                           , std::abs(domain->extent().z - domain->origin().z) };
      const double half_min_size_box = std::min( std::min(size_box.x,size_box.y) , size_box.z) / 2.0; 
      const auto xform = domain->xform();
      
      size_t n_cells = grid->number_of_cells();
      size_t nghosts = 0;
      for(size_t i=0;i<n_cells;i++) if(grid->is_ghost_cell(i)) nghosts += grid->cell(i).size();
      ldbg << nghosts << " ghost particles"<<std::endl;
      
#     pragma omp parallel
      {      
        ComputeBufferT buf;
        
#       pragma omp for schedule(dynamic)        
        for(size_t m=0;m<nmol;m++)
        {
          const unsigned int mtype = molecule_list->type(m);
          const unsigned int natoms = molecules->m_molecules.at(mtype).m_nb_atoms;
          
          // 1. initialize atom buffer
          for(unsigned int a=0;a<natoms;a++)
          {
            size_t c=0 , p=0;
            decode_cell_particle( molecule_list->atom_data(m,a) , c, p );
            buf.rx[a] = cells[c][field::rx][p];
            buf.ry[a] = cells[c][field::ry][p];
            buf.rz[a] = cells[c][field::rz][p];
            buf.fx[a] = 0.0;
            buf.fy[a] = 0.0;
            buf.fz[a] = 0.0;
            buf.ep[a] = 0.0;
            if constexpr ( has_virial_field )
            {
              buf.virial[a] = Mat3d{};
            }
          }
          
          // 2. compute bonds intramolecular forces
          const unsigned int nbonds = molecule_compute_parameters->m_molecules[mtype].m_nb_bonds;
          for(unsigned int i=0;i<nbonds;i++)
          {
            const uint32_t x = molecule_compute_parameters->m_molecules[mtype].m_bonds[i];
            const unsigned int a = ( x >> 16 ) & 0xFF;
            const unsigned int b = ( x >> 8 ) & 0xFF;
            const unsigned int pidx = x & 0xFF;

            const Vec3d ra = { buf.rx[a] , buf.ry[a] , buf.rz[a] };
            const Vec3d rb = { buf.rx[b] , buf.ry[b] , buf.rz[b] };
            const Vec3d r = xform * periodic_r_delta( ra , rb , size_box , half_min_size_box );
            
            const double norm_r = norm(r);
            assert(norm_r>0);

            // compute the pair dEp/dr , E
            std::pair<double,double> fe = intramolecular_quar( molecule_compute_parameters->m_func_params[pidx] , norm_r );

            // Compute energy
            const double e = fe.second;

            // Compute forces
            // F = dEp/dr nij
            const double dEp_dr = fe.first;
            const Vec3d F = ( r * dEp_dr ) / norm_r;

            buf.fx[a] += F.x;
            buf.fy[a] += F.y;
            buf.fz[a] += F.z;
            buf.ep[a] += e * 0.5;

            buf.fx[b] -= F.x;
            buf.fy[b] -= F.y;
            buf.fz[b] -= F.z;
            buf.ep[b] += e * 0.5;

            if constexpr ( has_virial_field )
            {
              const auto vir = tensor(F,r) * 0.5;
              buf.virial[a] += vir;
              buf.virial[b] += vir;
            }
          }
          

          // 3. compute bends(angles) intramolecular forces
          const unsigned int nbends = molecule_compute_parameters->m_molecules[mtype].m_nb_bends;
          for(unsigned int i=0;i<nbends;i++)
          {
            const uint32_t x = molecule_compute_parameters->m_molecules[mtype].m_bends[i];
            const unsigned int a = ( x >> 24 ) & 0xFF;
            const unsigned int b = ( x >> 16 ) & 0xFF;
            const unsigned int c = ( x >> 8 ) & 0xFF;
            const unsigned int pidx = x & 0xFF;

            const Vec3d ra = { buf.rx[a] , buf.ry[a] , buf.rz[a] };
            const Vec3d rb = { buf.rx[b] , buf.ry[b] , buf.rz[b] };
            const Vec3d rc = { buf.rx[c] , buf.ry[c] , buf.rz[c] };
            
            // first harm
            const Vec3d r1 = xform * periodic_r_delta( rb , ra , size_box , half_min_size_box );

            // second harm
            const Vec3d r2 = xform * periodic_r_delta( rb , rc , size_box , half_min_size_box );

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
            const auto fe = intramolecular_quar( molecule_compute_parameters->m_func_params[pidx] , theta );
            const double e = fe.second; //b.f_energy(theta);
            //---------------------------------------------------------------------------------------------

            //---------------------------------------------------------------------------------------------
            // Compute forces
            // F = - grad(Ep)
            const double dep_on_dtheta = fe.first; //b.f_forces(theta);

            const Vec3d F1 = pa_r1 * (-dep_on_dtheta);
            const Vec3d F2 = pb_r2 * (-dep_on_dtheta);
            const Vec3d Fo = - ( F1 + F2 );

            buf.fx[a] += F1.x;
            buf.fy[a] += F1.y;
            buf.fz[a] += F1.z;
            buf.ep[a] += e / 3.;

            buf.fx[b] += Fo.x;
            buf.fy[b] += Fo.y;
            buf.fz[b] += Fo.z;
            buf.ep[b] += e / 3.;

            buf.fx[c] += F2.x;
            buf.fy[c] += F2.y;
            buf.fz[c] += F2.z;
            buf.ep[c] += e / 3.;

            if constexpr ( has_virial_field )
            {
              const Mat3d virial1 = tensor (F1,r1) * 0.5;
              const Mat3d virial2 = tensor (F2,r2) * 0.5;
              const Mat3d virialo = virial1 + virial2;
              buf.virial[a] += virial1;
              buf.virial[b] += virialo;
              buf.virial[c] += virial2;
            }
          }

          // 4. compute torsion intramolecular forces
          const unsigned int ntorsions = molecule_compute_parameters->m_molecules[mtype].m_nb_torsions;
          for(unsigned int i=0;i<ntorsions;i++)
          {
            const uint64_t x = molecule_compute_parameters->m_molecules[mtype].m_torsions[i];
            const unsigned int a = ( x >> 32 ) & 0xFF;
            const unsigned int b = ( x >> 24 ) & 0xFF;
            const unsigned int c = ( x >> 16 ) & 0xFF;
            const unsigned int d = ( x >> 8 ) & 0xFF;
            const unsigned int pidx = x & 0xFF;

            const Vec3d r0 = { buf.rx[a] , buf.ry[a] , buf.rz[a] };
            const Vec3d r1 = { buf.rx[b] , buf.ry[b] , buf.rz[b] };
            const Vec3d r2 = { buf.rx[c] , buf.ry[c] , buf.rz[c] };
            const Vec3d r3 = { buf.rx[d] , buf.ry[d] , buf.rz[d] };
            
            const Vec3d r10 = xform * periodic_r_delta( r1 , r0 , size_box , half_min_size_box );
            const Vec3d r12 = xform * periodic_r_delta( r1 , r2 , size_box , half_min_size_box );
            const Vec3d r23 = xform * periodic_r_delta( r2 , r3 , size_box , half_min_size_box );
            
            //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
            // Compute energy
            const Vec3d V = cross(r10,r12);
            const Vec3d W = cross(r23,r12);

            assert(V!=(Vec3d{0,0,0}));//case r10 and r12 aligned
            assert(W!=(Vec3d{0,0,0}));//case r12 and r23 aligned

            const double phi = signum(dot(cross(V,W), r12)) * angle(V, W);
            const auto force_energy = intramolecular_cos( molecule_compute_parameters->m_func_params[pidx] , phi );
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

            buf.fx[a] += Fa.x;
            buf.fy[a] += Fa.y;
            buf.fz[a] += Fa.z;
            buf.ep[a] += e;

            buf.fx[b] += Fo.x;
            buf.fy[b] += Fo.y;
            buf.fz[b] += Fo.z;
            buf.ep[b] += e;

            buf.fx[c] += Fc.x;
            buf.fy[c] += Fc.y;
            buf.fz[c] += Fc.z;
            buf.ep[c] += e;
            
            buf.fx[d] += Fd.x;
            buf.fy[d] += Fd.y;
            buf.fz[d] += Fd.z;
            buf.ep[d] += e;

            if constexpr ( has_virial_field )
            {
              const Mat3d vira = tensor (Fa,r10) * 0.5;
              const Mat3d virc = tensor (Fc,r12) * 0.5;
              const Mat3d vird = tensor (Fd,r23) * 0.5;
              
              const Mat3d virac = vira + virc;
              const Mat3d vircd = virc + vird;

              buf.virial[a] += vira;
              buf.virial[b] += virac;
              buf.virial[c] += vircd;
              buf.virial[d] += vird;
            }
          }
     
          // 5. compute imporper intramolecular forces
          const unsigned int nimpropers = molecule_compute_parameters->m_molecules[mtype].m_nb_impropers;
          for(unsigned int i=0;i<nimpropers;i++)
          {
            const uint64_t x = molecule_compute_parameters->m_molecules[mtype].m_impropers[i];
            const unsigned int a = ( x >> 32 ) & 0xFF;
            const unsigned int b = ( x >> 24 ) & 0xFF;
            const unsigned int c = ( x >> 16 ) & 0xFF;
            const unsigned int d = ( x >> 8 ) & 0xFF;
            const unsigned int pidx = x & 0xFF;

            const Vec3d r0 = { buf.rx[a] , buf.ry[a] , buf.rz[a] };
            const Vec3d r1 = { buf.rx[b] , buf.ry[b] , buf.rz[b] };
            const Vec3d r2 = { buf.rx[c] , buf.ry[c] , buf.rz[c] };
            const Vec3d r3 = { buf.rx[d] , buf.ry[d] , buf.rz[d] };

            const Vec3d r10 = xform * periodic_r_delta( r1 , r0 , size_box , half_min_size_box );
            const Vec3d r12 = xform * periodic_r_delta( r1 , r2 , size_box , half_min_size_box );
            const Vec3d r23 = xform * periodic_r_delta( r2 , r3 , size_box , half_min_size_box );

            //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
            // Compute energy
            const Vec3d V = cross(r10,r12);
            const Vec3d W = cross(r23,r12);

            assert(V!=(Vec3d{0,0,0}));//case r10 and r12 aligned
            assert(W!=(Vec3d{0,0,0}));//case r12 and r23 aligned

            const double phi = signum(dot(cross(V,W), r12)) * angle(V, W);
            const auto force_energy = intramolecular_cos( molecule_compute_parameters->m_func_params[pidx] , phi );
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

            const Vec3d Fo = - (Fa+Fc+Fd);

            buf.fx[a] += Fa.x;
            buf.fy[a] += Fa.y;
            buf.fz[a] += Fa.z;
            buf.ep[a] += e;

            buf.fx[b] += Fo.x;
            buf.fy[b] += Fo.y;
            buf.fz[b] += Fo.z;
            buf.ep[b] += e;

            buf.fx[c] += Fc.x;
            buf.fy[c] += Fc.y;
            buf.fz[c] += Fc.z;
            buf.ep[c] += e;
            
            buf.fx[d] += Fd.x;
            buf.fy[d] += Fd.y;
            buf.fz[d] += Fd.z;
            buf.ep[d] += e;

            if constexpr ( has_virial_field )
            {
              const Mat3d vira = tensor (Fa,r10) * 0.5;
              const Mat3d virc = tensor (Fc,r12) * 0.5;
              const Mat3d vird = tensor (Fd,r23) * 0.5;
              const Mat3d virac = vira + virc;
              const Mat3d vircd = virc + vird;
              buf.virial[a] += vira;
              buf.virial[b] += virac;
              buf.virial[c] += vircd;
              buf.virial[d] += vird;
            }
          }

          // 6. compute intramolecular pair forces
//          #if 0
          const unsigned int npairs = molecule_compute_parameters->m_molecules[mtype].m_nb_pairs;
          for(unsigned int i=0;i<npairs;i++)
          {
            const uint64_t x = molecule_compute_parameters->m_molecules[mtype].m_pairs[i];
            const unsigned int wi = ( x >> 32 ) & 0xFF;
            const unsigned int a = ( x >> 24 ) & 0xFF;
            const unsigned int b = ( x >> 16 ) & 0xFF;
            const unsigned int pidx_rf = ( x >> 8 ) & 0xFF;
            const unsigned int pidx_ljexp6 = x & 0xFF;

            const double weight = wi * 0.5;

            const Vec3d ra = { buf.rx[a] , buf.ry[a] , buf.rz[a] };
            const Vec3d rb = { buf.rx[b] , buf.ry[b] , buf.rz[b] };
            const Vec3d r = xform * periodic_r_delta( ra , rb , size_box , half_min_size_box );
            
            const double d2 = norm2(r);
            const double d = sqrt(d2);
            
            double ljexp6_e=0.0 , ljexp6_de=0.0;
            if( d <= molecule_compute_parameters->m_ljexp6_params[pidx_ljexp6].m_rcut )
            {
              molecule_compute_parameters->m_ljexp6_params[pidx_ljexp6].compute_energy( d , ljexp6_e , ljexp6_de );
            }

            double rf_e=0.0 , rf_de=0.0;
            if( d <= molecule_compute_parameters->m_rf_params[pidx_rf].m_param.rc )
            {
              const double c1 = molecule_compute_parameters->m_rf_params[pidx_rf].m_c1;
              const double c2 = molecule_compute_parameters->m_rf_params[pidx_rf].m_c2;
              reaction_field_compute_energy( molecule_compute_parameters->m_rf_params[pidx_rf].m_param, c1, c2, d, rf_e, rf_de );
//              rf_e -= molecule_compute_parameters->m_pair_params[pidx].m_ecut_rf;
            }
                        
            double e = ljexp6_e + rf_e ;
            double de = ljexp6_de + rf_de ;
            de /= d;
            e *= weight;
            de *= weight;

            /*if( weight != 0.0 )
            {
              size_t cell_a=0 , p_a=0;
              decode_cell_particle( molecule_list->atom_data(m,a) , cell_a, p_a );
              size_t cell_b=0 , p_b=0;
              decode_cell_particle( molecule_list->atom_data(m,b) , cell_b, p_b );
              double c1 = molecule_compute_parameters->m_pair_params[pidx].m_c1;
              double c2 = molecule_compute_parameters->m_pair_params[pidx].m_c2;
              if(c1>c2){double t=c1;c1=c2;c2=t;}
              printf("C%d.p%d / C%d.p%d : w=% .6e : r=% .6e , c1=% .6e , c2=% .6e , de=% .6e\n",int(cell_a),int(p_a),int(cell_b),int(p_b), weight, d,c1,c2, de);
            }*/
            
            // force = dE * dr->
            const Vec3d F = de * r;

            buf.fx[a] += F.x;
            buf.fy[a] += F.y;
            buf.fz[a] += F.z;
            buf.ep[a] += 0.5 * e;

            buf.fx[b] -= F.x;
            buf.fy[b] -= F.y;
            buf.fz[b] -= F.z;
            buf.ep[b] += 0.5 * e;

            if constexpr ( has_virial_field )
            {
              const Mat3d vir = tensor( F, r ) * -0.5;
              buf.virial[a] += vir;
              buf.virial[b] += vir;
            }
          }
//          #endif
        
          // 7. store results back to particles
          for(unsigned int a=0;a<natoms;a++)
          {
            size_t c=0 , p=0;
            decode_cell_particle( molecule_list->atom_data(m,a) , c, p );
            buf.rx[a] = cells[c][field::rx][p];
            buf.ry[a] = cells[c][field::ry][p];
            buf.rz[a] = cells[c][field::rz][p];
            cells[c][field::fx][p] += buf.fx[a];
            cells[c][field::fy][p] += buf.fy[a];
            cells[c][field::fz][p] += buf.fz[a];
            cells[c][field::ep][p] += buf.ep[a];
            if constexpr ( has_virial_field )
            {
              cells[c][field::virial][p] += buf.virial[a];
            }
          }
       
        } // end of loop on molecules
       
      } // end of parallel section

    }

  };

  template<class GridT> using IntramolecularForceComputeTmpl = IntramolecularForceCompute<GridT>;

  // === register factories ===
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory("intramolecular_compute", make_grid_variant_operator< IntramolecularForceComputeTmpl > );
  }

}
