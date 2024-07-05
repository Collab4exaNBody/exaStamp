#pragma once

#include <exaStamp/molecule/molecule_list.h>
#include <exaStamp/molecule/molecule_species.h>
#include <exaStamp/molecule/molecule_compute_buffer.h>
#include <exaStamp/molecule/periodic_r_delta.h>
#include <exaStamp/molecule/potential_functional.h>

#include <onika/cuda/cuda.h>
#include <onika/parallel/parallel_for.h>
#include <onika/cuda/ro_shallow_copy.h>

namespace exaStamp
{

  template<class CellsAccessorT, class RXT, class RYT, class RZT, class FXT, class FYT, class FZT, class ET, class VirialT = void >
  struct IntramolecularForceFunctor
  {    
    static inline constexpr bool has_virial_field = ! std::is_same_v< VirialT , void >;
    using VirialMemberT = std::conditional_t< has_virial_field , VirialT, bool >;
    using ComputeBufferT = std::conditional_t< has_virial_field , MoleculeVirialComputeBuffer, MoleculeComputeBuffer >;

    onika::cuda::ro_shallow_copy_t< MoleculeLists > m_molecule_list = {};
    onika::cuda::ro_shallow_copy_t< MoleculeSpeciesVector > m_molecules = {};
    onika::cuda::ro_shallow_copy_t< MoleculeSetComputeParams > molecule_compute_parameters = {};
    const Mat3d m_xform = { 1.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 1.0 };
    const Vec3d m_size_box = { 0.0 , 0.0 , 0.0 };
    const double m_half_min_size_box = 0.0;
    
    CellsAccessorT m_cells = {};
    RXT m_rx = {};
    RYT m_ry = {};
    RZT m_rz = {};
    FXT m_fx = {};
    FYT m_fy = {};
    FZT m_fz = {};
    ET m_ep = {};
    VirialMemberT m_virial = {};
    
    ONIKA_HOST_DEVICE_FUNC
    inline void operator () ( size_t m ) const
    {
      const unsigned int mtype = m_molecule_list.type(m);
//      const unsigned int natoms = m_molecules.m_molecules[mtype].m_nb_atoms;

      // 2. compute bonds intramolecular forces
      //    A--B
      const unsigned int nbonds = molecule_compute_parameters.m_molecules[mtype].m_nb_bonds;
      for(unsigned int i=0;i<nbonds;i++)
      {
        // bonds consist of 3 8-bits encoded values packed in a single 32-bit word (var. x), containing :
        // - atom A index in molecule (var. a), 
        // - atom B index in molecule (var. b),
        // - parameter pack index (var. pidx)
        const uint32_t x = molecule_compute_parameters.m_molecules[mtype].m_bonds[i];
        const unsigned int a = ( x >> 16 ) & 0xFF;
        const unsigned int b = ( x >> 8 ) & 0xFF;
        const unsigned int pidx = x & 0xFF;

        size_t ca=0 , pa=0;
        decode_cell_particle( m_molecule_list.atom_data(m,a) , ca, pa );
        const Vec3d ra = { m_cells[ca][m_rx][pa], m_cells[ca][m_ry][pa], m_cells[ca][m_rz][pa] };

        size_t cb=0 , pb=0;
        decode_cell_particle( m_molecule_list.atom_data(m,b) , cb, pb );
        const Vec3d rb = { m_cells[cb][m_rx][pb], m_cells[cb][m_ry][pb], m_cells[cb][m_rz][pb] };
        
        const Vec3d r = m_xform * periodic_r_delta( ra , rb , m_size_box , m_half_min_size_box );
        
        const double norm_r = norm(r);
        assert(norm_r>0);

        // compute the pair dEp/dr , E
        const auto [ dEp_dr , e ] = intramolecular_quar( molecule_compute_parameters.m_func_params[pidx] , norm_r );

        // Compute energy
        //const double e = fe.second;

        // Compute forces
        // F = dEp/dr nij
        //const double dEp_dr = fe.first;
        const Vec3d F = ( r * dEp_dr ) / norm_r;

        m_cells[ca][m_fx][pa] += F.x;
        m_cells[ca][m_fy][pa] += F.y;
        m_cells[ca][m_fz][pa] += F.z;
        m_cells[ca][m_ep][pa] += e * 0.5;

        m_cells[cb][m_fx][pb] -= F.x;
        m_cells[cb][m_fy][pb] -= F.y;
        m_cells[cb][m_fz][pb] -= F.z;
        m_cells[cb][m_ep][pb] += e * 0.5;

        if constexpr ( has_virial_field )
        {
          const auto vir = tensor(F,r) * 0.5;
          m_cells[ca][m_virial][pa] += vir;
          m_cells[cb][m_virial][pb] += vir;
        }
      }
      

      // 3. compute bends(angles) intramolecular forces
      //    A   C
      //    \  /
      //     B
      const unsigned int nbends = molecule_compute_parameters.m_molecules[mtype].m_nb_bends;
      for(unsigned int i=0;i<nbends;i++)
      {
        // bends (or angles) consist of 4 8-bits encoded values packed in a single 32-bit word (var. x), containing :
        // - atom A index in molecule (var. a), 
        // - atom B index in molecule (var. b),
        // - atom C index in molecule (var. c),
        // - parameter pack index (var. pidx)
        const uint32_t x = molecule_compute_parameters.m_molecules[mtype].m_bends[i];
        const unsigned int a = ( x >> 24 ) & 0xFF;
        const unsigned int b = ( x >> 16 ) & 0xFF;
        const unsigned int c = ( x >> 8 ) & 0xFF;
        const unsigned int pidx = x & 0xFF;

        size_t ca=0 , pa=0;
        decode_cell_particle( m_molecule_list.atom_data(m,a) , ca, pa );
        const Vec3d ra = { m_cells[ca][m_rx][pa], m_cells[ca][m_ry][pa], m_cells[ca][m_rz][pa] };

        size_t cb=0 , pb=0;
        decode_cell_particle( m_molecule_list.atom_data(m,b) , cb, pb );
        const Vec3d rb = { m_cells[cb][m_rx][pb], m_cells[cb][m_ry][pb], m_cells[cb][m_rz][pb] };

        size_t cc=0 , pc=0;
        decode_cell_particle( m_molecule_list.atom_data(m,c) , cc, pc );
        const Vec3d rc = { m_cells[cc][m_rx][pc], m_cells[cc][m_ry][pc], m_cells[cc][m_rz][pc] };
        
        // first harm
        const Vec3d r1 = m_xform * periodic_r_delta( rb , ra , m_size_box , m_half_min_size_box );

        // second harm
        const Vec3d r2 = m_xform * periodic_r_delta( rb , rc , m_size_box , m_half_min_size_box );

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
        const auto [ dep_on_dtheta , e ] = intramolecular_quar( molecule_compute_parameters.m_func_params[pidx] , theta );
        //const double e = fe.second; //b.f_energy(theta);
        //---------------------------------------------------------------------------------------------

        //---------------------------------------------------------------------------------------------
        // Compute forces
        // F = - grad(Ep)
        //const double dep_on_dtheta = fe.first; //b.f_forces(theta);

        const Vec3d F1 = pa_r1 * (-dep_on_dtheta);
        const Vec3d F2 = pb_r2 * (-dep_on_dtheta);
        const Vec3d Fo = - ( F1 + F2 );

        m_cells[ca][m_fx][pa] += F1.x;
        m_cells[ca][m_fy][pa] += F1.y;
        m_cells[ca][m_fz][pa] += F1.z;
        m_cells[ca][m_ep][pa] += e / 3.;

        m_cells[cb][m_fx][pb] += Fo.x;
        m_cells[cb][m_fy][pb] += Fo.y;
        m_cells[cb][m_fz][pb] += Fo.z;
        m_cells[cb][m_ep][pb] += e / 3.;

        m_cells[cc][m_fx][pc] += F2.x;
        m_cells[cc][m_fy][pc] += F2.y;
        m_cells[cc][m_fz][pc] += F2.z;
        m_cells[cc][m_ep][pc] += e / 3.;

        if constexpr ( has_virial_field )
        {
          const Mat3d virial1 = tensor (F1,r1) * 0.5;
          const Mat3d virial2 = tensor (F2,r2) * 0.5;
          const Mat3d virialo = virial1 + virial2;
          m_cells[ca][m_virial][pa] += virial1;
          m_cells[cb][m_virial][pb] += virialo;
          m_cells[cc][m_virial][pc] += virial2;
        }
      }


      // 4. compute torsion intramolecular forces
      //    A     D
      //    \    /
      //     B__C
      const unsigned int ntorsions = molecule_compute_parameters.m_molecules[mtype].m_nb_torsions;
      for(unsigned int i=0;i<ntorsions;i++)
      {
        // torsions consist of 5 8-bits encoded values packed in a single 64-bit word (var. x), containing :
        // - atom A index in molecule (var. a), 
        // - atom B index in molecule (var. b),
        // - atom C index in molecule (var. c),
        // - atom D index in molecule (var. d),
        // - parameter pack index (var. pidx)
        const uint64_t x = molecule_compute_parameters.m_molecules[mtype].m_torsions[i];
        const unsigned int a = ( x >> 32 ) & 0xFF;
        const unsigned int b = ( x >> 24 ) & 0xFF;
        const unsigned int c = ( x >> 16 ) & 0xFF;
        const unsigned int d = ( x >> 8 ) & 0xFF;
        const unsigned int pidx = x & 0xFF;

        size_t ca=0 , pa=0;
        decode_cell_particle( m_molecule_list.atom_data(m,a) , ca, pa );
        const Vec3d r0 = { m_cells[ca][m_rx][pa], m_cells[ca][m_ry][pa], m_cells[ca][m_rz][pa] };

        size_t cb=0 , pb=0;
        decode_cell_particle( m_molecule_list.atom_data(m,b) , cb, pb );
        const Vec3d r1 = { m_cells[cb][m_rx][pb], m_cells[cb][m_ry][pb], m_cells[cb][m_rz][pb] };

        size_t cc=0 , pc=0;
        decode_cell_particle( m_molecule_list.atom_data(m,c) , cc, pc );
        const Vec3d r2 = { m_cells[cc][m_rx][pc], m_cells[cc][m_ry][pc], m_cells[cc][m_rz][pc] };

        size_t cd=0 , pd=0;
        decode_cell_particle( m_molecule_list.atom_data(m,d) , cd, pd );
        const Vec3d r3 = { m_cells[cd][m_rx][pd], m_cells[cd][m_ry][pd], m_cells[cd][m_rz][pd] };
        
        const Vec3d r10 = m_xform * periodic_r_delta( r1 , r0 , m_size_box , m_half_min_size_box );
        const Vec3d r12 = m_xform * periodic_r_delta( r1 , r2 , m_size_box , m_half_min_size_box );
        const Vec3d r23 = m_xform * periodic_r_delta( r2 , r3 , m_size_box , m_half_min_size_box );
        
        //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
        // Compute energy
        const Vec3d V = cross(r10,r12);
        const Vec3d W = cross(r23,r12);

        assert(V!=(Vec3d{0,0,0}));//case r10 and r12 aligned
        assert(W!=(Vec3d{0,0,0}));//case r12 and r23 aligned

        const double phi = signum(dot(cross(V,W), r12)) * angle(V, W);
        const auto force_energy = intramolecular_cos( molecule_compute_parameters.m_func_params[pidx] , phi );
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

        m_cells[ca][m_fx][pa] += Fa.x;
        m_cells[ca][m_fy][pa] += Fa.y;
        m_cells[ca][m_fz][pa] += Fa.z;
        m_cells[ca][m_ep][pa] += e;

        m_cells[cb][m_fx][pb] += Fo.x;
        m_cells[cb][m_fy][pb] += Fo.y;
        m_cells[cb][m_fz][pb] += Fo.z;
        m_cells[cb][m_ep][pb] += e;

        m_cells[cc][m_fx][pc] += Fc.x;
        m_cells[cc][m_fy][pc] += Fc.y;
        m_cells[cc][m_fz][pc] += Fc.z;
        m_cells[cc][m_ep][pc] += e;
        
        m_cells[cd][m_fx][pd] += Fd.x;
        m_cells[cd][m_fy][pd] += Fd.y;
        m_cells[cd][m_fz][pd] += Fd.z;
        m_cells[cd][m_ep][pd] += e;

        if constexpr ( has_virial_field )
        {
          const Mat3d vira = tensor (Fa,r10) * 0.5;
          const Mat3d virc = tensor (Fc,r12) * 0.5;
          const Mat3d vird = tensor (Fd,r23) * 0.5;
          
          const Mat3d virac = vira + virc;
          const Mat3d vircd = virc + vird;

          m_cells[ca][m_virial][pa] += vira;
          m_cells[cb][m_virial][pb] += virac;
          m_cells[cc][m_virial][pc] += vircd;
          m_cells[cd][m_virial][pd] += vird;
        }
      }
 
      // 5. compute imporper intramolecular forces
      //    A  D
      //    \ /
      //     B--C
      const unsigned int nimpropers = molecule_compute_parameters.m_molecules[mtype].m_nb_impropers;
      for(unsigned int i=0;i<nimpropers;i++)
      {
        // impropers consist of 5 8-bits encoded values packed in a single 64-bit word (var. x), containing :
        // - atom A index in molecule (var. a), 
        // - atom B index in molecule (var. b),
        // - atom C index in molecule (var. c),
        // - atom D index in molecule (var. d),
        // - parameter pack index (var. pidx)
        const uint64_t x = molecule_compute_parameters.m_molecules[mtype].m_impropers[i];
        const unsigned int a = ( x >> 32 ) & 0xFF;
        const unsigned int b = ( x >> 24 ) & 0xFF;
        const unsigned int c = ( x >> 16 ) & 0xFF;
        const unsigned int d = ( x >> 8 ) & 0xFF;
        const unsigned int pidx = x & 0xFF;

        size_t ca=0 , pa=0;
        decode_cell_particle( m_molecule_list.atom_data(m,a) , ca, pa );
        const Vec3d r0 = { m_cells[ca][m_rx][pa], m_cells[ca][m_ry][pa], m_cells[ca][m_rz][pa] };

        size_t cb=0 , pb=0;
        decode_cell_particle( m_molecule_list.atom_data(m,b) , cb, pb );
        const Vec3d r1 = { m_cells[cb][m_rx][pb], m_cells[cb][m_ry][pb], m_cells[cb][m_rz][pb] };

        size_t cc=0 , pc=0;
        decode_cell_particle( m_molecule_list.atom_data(m,c) , cc, pc );
        const Vec3d r2 = { m_cells[cc][m_rx][pc], m_cells[cc][m_ry][pc], m_cells[cc][m_rz][pc] };

        size_t cd=0 , pd=0;
        decode_cell_particle( m_molecule_list.atom_data(m,d) , cd, pd );
        const Vec3d r3 = { m_cells[cd][m_rx][pd], m_cells[cd][m_ry][pd], m_cells[cd][m_rz][pd] };

        const Vec3d r10 = m_xform * periodic_r_delta( r1 , r0 , m_size_box , m_half_min_size_box );
        const Vec3d r12 = m_xform * periodic_r_delta( r1 , r2 , m_size_box , m_half_min_size_box );
        const Vec3d r23 = m_xform * periodic_r_delta( r2 , r3 , m_size_box , m_half_min_size_box );

        //-----------------------------COMPUTE ENERGY AND FORCE----------------------------------
        // Compute energy
        const Vec3d V = cross(r10,r12);
        const Vec3d W = cross(r23,r12);

        assert(V!=(Vec3d{0,0,0}));//case r10 and r12 aligned
        assert(W!=(Vec3d{0,0,0}));//case r12 and r23 aligned

        const double phi = signum(dot(cross(V,W), r12)) * angle(V, W);
        const auto force_energy = intramolecular_cos( molecule_compute_parameters.m_func_params[pidx] , phi );
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

        m_cells[ca][m_fx][pa] += Fa.x;
        m_cells[ca][m_fy][pa] += Fa.y;
        m_cells[ca][m_fz][pa] += Fa.z;
        m_cells[ca][m_ep][pa] += e;

        m_cells[cb][m_fx][pb] += Fo.x;
        m_cells[cb][m_fy][pb] += Fo.y;
        m_cells[cb][m_fz][pb] += Fo.z;
        m_cells[cb][m_ep][pb] += e;

        m_cells[cc][m_fx][pc] += Fc.x;
        m_cells[cc][m_fy][pc] += Fc.y;
        m_cells[cc][m_fz][pc] += Fc.z;
        m_cells[cc][m_ep][pc] += e;
        
        m_cells[cd][m_fx][pd] += Fd.x;
        m_cells[cd][m_fy][pd] += Fd.y;
        m_cells[cd][m_fz][pd] += Fd.z;
        m_cells[cd][m_ep][pd] += e;

        if constexpr ( has_virial_field )
        {
          const Mat3d vira = tensor (Fa,r10) * 0.5;
          const Mat3d virc = tensor (Fc,r12) * 0.5;
          const Mat3d vird = tensor (Fd,r23) * 0.5;
          const Mat3d virac = vira + virc;
          const Mat3d vircd = virc + vird;
          m_cells[ca][m_virial][pa] += vira;
          m_cells[cb][m_virial][pb] += virac;
          m_cells[cc][m_virial][pc] += vircd;
          m_cells[cd][m_virial][pd] += vird;
        }
      }

      // 6. compute intramolecular pair forces
      const unsigned int npairs = molecule_compute_parameters.m_molecules[mtype].m_nb_pairs;
      for(unsigned int i=0;i<npairs;i++)
      {
        const uint64_t x = molecule_compute_parameters.m_molecules[mtype].m_pairs[i];
        const unsigned int wi = ( x >> 32 ) & 0xFF;
        const unsigned int a = ( x >> 24 ) & 0xFF;
        const unsigned int b = ( x >> 16 ) & 0xFF;
        const unsigned int pidx_rf = ( x >> 8 ) & 0xFF;
        const unsigned int pidx_ljexp6 = x & 0xFF;

        const double weight = wi * 0.5;

        size_t ca=0 , pa=0;
        decode_cell_particle( m_molecule_list.atom_data(m,a) , ca, pa );
        const Vec3d ra = { m_cells[ca][m_rx][pa], m_cells[ca][m_ry][pa], m_cells[ca][m_rz][pa] };

        size_t cb=0 , pb=0;
        decode_cell_particle( m_molecule_list.atom_data(m,b) , cb, pb );
        const Vec3d rb = { m_cells[cb][m_rx][pb], m_cells[cb][m_ry][pb], m_cells[cb][m_rz][pb] };

        const Vec3d r = m_xform * periodic_r_delta( ra , rb , m_size_box , m_half_min_size_box );
        
        const double d2 = norm2(r);
        const double d = sqrt(d2);
        
        double ljexp6_e=0.0 , ljexp6_de=0.0;
        if( d <= molecule_compute_parameters.m_ljexp6_params[pidx_ljexp6].m_rcut )
        {
          molecule_compute_parameters.m_ljexp6_params[pidx_ljexp6].compute_energy( d , ljexp6_e , ljexp6_de );
        }

        double rf_e=0.0 , rf_de=0.0;
        if( d <= molecule_compute_parameters.m_rf_params[pidx_rf].m_param.rc )
        {
          const double c1 = molecule_compute_parameters.m_rf_params[pidx_rf].m_c1;
          const double c2 = molecule_compute_parameters.m_rf_params[pidx_rf].m_c2;
          reaction_field_compute_energy( molecule_compute_parameters.m_rf_params[pidx_rf].m_param, c1 * c2, d, rf_e, rf_de );
        }
                    
        double e = ljexp6_e + rf_e ;
        double de = ljexp6_de + rf_de ;
        de /= d;
        e *= weight;
        de *= weight;
        
        // force = dE * dr->
        const Vec3d F = de * r;

        m_cells[ca][m_fx][pa] += F.x;
        m_cells[ca][m_fy][pa] += F.y;
        m_cells[ca][m_fz][pa] += F.z;
        m_cells[ca][m_ep][pa] += 0.5 * e;

        m_cells[cb][m_fx][pb] -= F.x;
        m_cells[cb][m_fy][pb] -= F.y;
        m_cells[cb][m_fz][pb] -= F.z;
        m_cells[cb][m_ep][pb] += 0.5 * e;

        if constexpr ( has_virial_field )
        {
          const Mat3d vir = tensor( F, r ) * -0.5;
          m_cells[ca][m_virial][pa] += vir;
          m_cells[cb][m_virial][pb] += vir;
        }
      }

    }
    
  };

  template<class CellsAccessorT, class RXT, class RYT, class RZT, class FXT, class FYT, class FZT, class ET, class VirialT >
  inline
  IntramolecularForceFunctor<CellsAccessorT,RXT,RYT,RZT,FXT,FYT,FZT,ET,VirialT>
  make_intramolecular_functor( const MoleculeLists& ml, const MoleculeSpeciesVector& mols, const MoleculeSetComputeParams& mol_parms,
                               const Mat3d& xform , const Vec3d& size_box, double half_box_size ,
                               CellsAccessorT cells, const onika::FlatTuple<RXT,RYT,RZT,FXT,FYT,FZT,ET,VirialT>& fields )
  {
    return { ml , mols , mol_parms , xform , size_box, half_box_size , cells
           , fields.get(onika::tuple_index<0>)
           , fields.get(onika::tuple_index<1>)
           , fields.get(onika::tuple_index<2>)
           , fields.get(onika::tuple_index<3>)
           , fields.get(onika::tuple_index<4>)
           , fields.get(onika::tuple_index<5>)
           , fields.get(onika::tuple_index<6>)
           , fields.get(onika::tuple_index<7>)
           };
  }

  template<class CellsAccessorT, class RXT, class RYT, class RZT, class FXT, class FYT, class FZT, class ET >
  inline
  IntramolecularForceFunctor<CellsAccessorT,RXT,RYT,RZT,FXT,FYT,FZT,ET>
  make_intramolecular_functor( const MoleculeLists& ml, const MoleculeSpeciesVector& mols, const MoleculeSetComputeParams& mol_parms, 
                               const Mat3d& xform , const Vec3d& size_box, double half_box_size ,
                               CellsAccessorT cells, const onika::FlatTuple<RXT,RYT,RZT,FXT,FYT,FZT,ET>& fields )
  {
    return { ml , mols , mol_parms , xform , size_box , half_box_size , cells
           , fields.get(onika::tuple_index<0>)
           , fields.get(onika::tuple_index<1>)
           , fields.get(onika::tuple_index<2>)
           , fields.get(onika::tuple_index<3>)
           , fields.get(onika::tuple_index<4>)
           , fields.get(onika::tuple_index<5>)
           , fields.get(onika::tuple_index<6>)
           };
  }

}

namespace onika
{
  namespace parallel
  {  
    template<class CellsAccessorT, class RXT, class RYT, class RZT, class FXT, class FYT, class FZT, class ET, class VirialT >
    struct ParallelForFunctorTraits< exaStamp::IntramolecularForceFunctor<CellsAccessorT,RXT,RYT,RZT,FXT,FYT,FZT,ET,VirialT> >
    {      
      static inline constexpr bool CudaCompatible = true;
    };
  }
}


