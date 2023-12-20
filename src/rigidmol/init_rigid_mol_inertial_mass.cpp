#include <exaStamp/particle_species/particle_specie_yaml.h>
#include <exanb/core/basic_types_stream.h>

#include <exanb/core/operator.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/log.h>
#include <exanb/core/yaml_utils.h>

#include <exanb/core/basic_types_def.h>
#include <exanb/core/math_utils.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>

#ifdef XSTAMP_USE_EIGEN3
#include <Eigen/Dense>
#endif

namespace exaStamp
{

  using namespace exanb;

  struct InitRigidMolInertialMass : public OperatorNode
  {
    ADD_SLOT( ParticleSpecies , species , INPUT_OUTPUT , REQUIRED );

    inline void execute () override final
    {
      for(unsigned int a=0;a<species->size();a++)
      {
        ParticleSpecie& mol = species->at(a);
        if( mol.m_rigid_atom_count < 1 )
        {
          lerr << "rigid atom count must be >= 1"<<std::endl;
          std::abort();
        }
        /*if( mol.m_rigid_atom_count != mol.m_rigid_atom_names.size() )
        {
          lerr << "insconsistent rigid atom names count. count="<<mol.m_rigid_atom_count<<", names=";
          for(const auto& s:mol.m_rigid_atom_names) lerr<<s<<" ";
          lerr<<std::endl;
          std::abort();          
        }*/
        
        mol.m_minert = { 0. , 0. , 0. };
        for(unsigned int i=0;i<mol.m_rigid_atom_count;i++)
        {
          mol.m_rigid_atoms[0].m_pos = { 0., 0., 0. };
        }
        if( mol.m_rigid_atom_count > 1 )
        {
          Mat3d minert_tensor { 0.,0.,0.,  0.,0.,0.,  0.,0.,0. };
          for(unsigned int i=0;i<mol.m_rigid_atom_count;i++)
          {
            size_t site_atom_type = mol.m_rigid_atoms[i].m_atom_type;
            if( site_atom_type >= species->size() )
            {
              lerr << "bad rigid molecule site atom type "<<site_atom_type <<std::endl;
              std::abort();
            }
            double site_mass = species->at( site_atom_type ).m_mass;
            Vec3d site_relative_pos = mol.m_rigid_atoms[i].m_pos;
            //elements diagonaux du tenseur d'inertie
            minert_tensor.m11 += site_mass * ( site_relative_pos.y*site_relative_pos.y + site_relative_pos.z*site_relative_pos.z);
            minert_tensor.m22 += site_mass * ( site_relative_pos.x*site_relative_pos.x + site_relative_pos.z*site_relative_pos.z);
            minert_tensor.m33 += site_mass * ( site_relative_pos.y*site_relative_pos.y + site_relative_pos.x*site_relative_pos.x);
            //autres elements
            minert_tensor.m12 += -site_mass * site_relative_pos.x * site_relative_pos.y;
            minert_tensor.m13 += -site_mass * site_relative_pos.x * site_relative_pos.z;
            minert_tensor.m23 += -site_mass * site_relative_pos.y * site_relative_pos.z;
          }
          //on complete le tenseur d'inertie
          minert_tensor.m21 = minert_tensor.m12;
          minert_tensor.m31 = minert_tensor.m13;
          minert_tensor.m32 = minert_tensor.m23;
          //on recupere les moments d'inertie
          mol.m_minert = eigenvaluescompute(minert_tensor);
          //on recupere la matrice de passage de la base initiale a la base de vecteurs propres
          const Mat3d tmat = transformmatcompute(minert_tensor);
          //on calcule les positions des atomes dans la base de vecteurs propres
          for(unsigned int i=0;i<mol.m_rigid_atom_count;i++)
          {
            size_t site_atom_type = mol.m_rigid_atoms[i].m_atom_type;
            if( site_atom_type >= species->size() )
            {
              lerr << "bad rigid molecule site atom type "<<site_atom_type <<std::endl;
              std::abort();
            }
            Vec3d pos = mol.m_rigid_atoms[i].m_pos;
            mol.m_rigid_atoms[i].m_pos.x = tmat.m11 * pos.x + tmat.m12 * pos.y + tmat.m13 * pos.z;
            mol.m_rigid_atoms[i].m_pos.y = tmat.m21 * pos.x + tmat.m22 * pos.y + tmat.m23 * pos.z;
            mol.m_rigid_atoms[i].m_pos.z = tmat.m31 * pos.z + tmat.m32 * pos.y + tmat.m33 * pos.z;
          }
        }
      }
    }

    static inline Vec3d eigenvaluescompute( const Mat3d& M)
    {
#     ifdef XSTAMP_USE_EIGEN3
      using namespace Eigen;
      MatrixXd Mat(3,3);
      Mat(0,0) = M.m11;
      Mat(0,1) = M.m12;
      Mat(0,2) = M.m13;
      Mat(1,0) = M.m21;
      Mat(1,1) = M.m22;
      Mat(1,2) = M.m23;
      Mat(2,0) = M.m31;
      Mat(2,1) = M.m32;
      Mat(2,2) = M.m33;
      SelfAdjointEigenSolver<Matrix3d> eigensolver(Mat);
      auto eigenval = eigensolver.eigenvalues();
      return { eigenval(0) , eigenval(1) , eigenval(2) };
#     else
      Vec3d eigenvec[3];
      double eigenval[3];
      symmetric_matrix_eigensystem(M, eigenvec, eigenval);
      return { eigenval[0] , eigenval[1] , eigenval[2] };
#     endif
    }

    static inline Mat3d transformmatcompute( const Mat3d& M)
    {
#     ifdef XSTAMP_USE_EIGEN3
      using namespace Eigen;
      MatrixXd Mat(3,3);
      Mat(0,0) = M.m11;
      Mat(0,1) = M.m12;
      Mat(0,2) = M.m13;
      Mat(1,0) = M.m21;
      Mat(1,1) = M.m22;
      Mat(1,2) = M.m23;
      Mat(2,0) = M.m31;
      Mat(2,1) = M.m32;
      Mat(2,2) = M.m33;
      SelfAdjointEigenSolver<Matrix3d> eigensolver(Mat);
      auto eigenvec =  eigensolver.eigenvectors();
      auto matpass = eigenvec.inverse();
      Vec3d vecx {matpass.col(0)(0), matpass.col(0)(1), matpass.col(0)(2)};
      Vec3d vecy {matpass.col(1)(0), matpass.col(1)(1), matpass.col(1)(2)};
      Vec3d vecz {matpass.col(2)(0), matpass.col(2)(1), matpass.col(2)(2)};
      Mat3d matpassage;
      matpassage = make_mat3d(vecx, vecy, vecz); 
      return { matpassage };
#     else
      Vec3d Q[3];
      double w[3];
      symmetric_matrix_eigensystem(M, Q, w);
      return inverse( make_mat3d(Q[0], Q[1], Q[2]) );
#     endif
    }

  };

  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "init_rigid_mol_inertial_mass", make_simple_operator< InitRigidMolInertialMass > );
  }

}

