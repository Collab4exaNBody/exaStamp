#include <exanb/core/basic_types.h>

#include <exanb/defbox/deformation.h>
#include <exanb/defbox/deformation_math.h>

#include <exaStamp/mechanical/cell_particles_local_mechanical_metrics.h>
#include <exaStamp/mechanical/basic_algebra_for_mechanics.h>

#include <cmath>
#include <algorithm>

namespace exaStamp
{

  static inline int sign(const double& test) {return (test >= 0) - (test < 0);}

  // Quelques fonctions de pondération, à tester
  // Fonction d'atténuation de Wendland normalisée à 10
  static inline double Wendland_C2(const double& x) { return (1+2*x)*pow(1-x/2., 4.);  }
  // Fonction d'atténuation de Wendland, à tester  
  static inline double Wendland_weight(const double& dist, const double& DistWeight) {    
    double sigma = DistWeight/3.;    
    if(dist < 4*sigma) {    
      double alpha = 10*21/16./M_PI;
      double beta2 = 4/15.;
      return alpha*beta2/pow(sigma, 3.)*Wendland_C2(dist/sigma*sqrt(beta2));
    } else {
      return 0.;
    }
  } 
  // Fonction d'atténuation article TATB (P. Lafourcade, C. Denoual, J.-B. Maillet JPCC 2018)
  static inline double calc_weight(const double& DistCur, const double& DistWeight) {	  
    if(DistCur/DistWeight <= 0.5) {
      return 1.0 - 6.0 * pow((DistCur/(DistWeight)),2.0) + 6.0 * pow((DistCur/(DistWeight)),3.0);
    } else if ((DistCur/DistWeight > 0.5) && (DistCur/DistWeight < 1.0)) {
      return 2.0 - 6.0 * (DistCur/(DistWeight)) + 6.0 * pow((DistCur/(DistWeight)),2.0) - 2.0 * pow((DistCur/(DistWeight)),3.0);
    } else {
      return 0;
    }
  }
  
  template< class GridT >
  struct DeformationGradientComputeOp
  {
    using CellsT = decltype( GridT{}.cells() );
    const CellsT m_cells_t0;

    GridParticleLocalMechanicalMetrics & m_local_mechanical_data;

    bool compute_static_measures;
    bool compute_dynamic_measures;
    
    //    std::vector< std::vector<Mat3d> > & m_tensor_data;
    const Mat3d m_xform_t0;
    const Mat3d m_xform;
    const Mat3d m_lattice;
    const double m_rrDef;

    double (*weight_function)(const double&, const double&) = calc_weight;
    
    template<class ComputeBufferT, class CellParticlesT /*, class GridParticleLocksT, class ParticleLockT*/>
    inline void operator ()
    (
     // number of neighbors
     size_t n,

     // buffer holding attached data to particle's neighbors
     const ComputeBufferT& tab,
                
     // data and locks accessors for neighbors
     CellParticlesT cells
     //        , GridParticleLocksT locks
     //        , ParticleLockT & particle_lock
     ) const
    {
      GridParticleLocalMechanicalMetrics & local_mechanical_data = m_local_mechanical_data;
      //      std::vector< std::vector<Mat3d> > & tensor_data = m_tensor_data;

      // Static measures      
      Mat3d F,E,R,U,skewR;
      Vec3d mu;

      Mat3d tensorAF; // intitialized to all zero
      Mat3d invtensorAF; // intitialized to all zero        
      Mat3d tensorBF;      

      Vec3d deltaPosInit;
      Vec3d deltaPosCour;
      Vec3d deltaPosInitCris;
      Vec3d deltaPosCourCris;                
      Vec3d dtr;
      
      // Slip basis
      Mat3d A,Ltens,Ntens;
      Vec3d deltaCourInit,lvec,mvec,nvec,s;
      int slipped_neighbors_count=0;
      
      // Cinematic measures;
      Mat3d L,skewL;
      Vec3d phi;

      Mat3d tensorAL; // intitialized to all zero
      Mat3d invtensorAL; // intitialized to all zero        
      Mat3d tensorBL;      

      Vec3d deltaVitCour;
      
      Mat3d hht = m_xform * m_lattice;
      Mat3d hh0 = m_xform_t0 * m_lattice;

      Vec3d tmp;
      tmp.x = m_cells_t0[tab.cell][field::rx][tab.part];
      tmp.y = m_cells_t0[tab.cell][field::ry][tab.part];
      tmp.z = m_cells_t0[tab.cell][field::rz][tab.part];

      Vec3d PosInit = m_xform_t0 * tmp;

      Vec3d PosInitVois;
      Mat3d I = make_identity_matrix();
      
      // loop over neighborhood
      for(size_t i=0;i<n;i++)
	{

	  // Relative positions in the reference state
	  tmp.x = tab.ext.rx0[i];
	  tmp.y = tab.ext.ry0[i];
	  tmp.z = tab.ext.rz0[i];
	  PosInitVois = m_xform_t0 * tmp;
	  deltaPosInit = PosInitVois - PosInit;

	  // Relative positions in the current state	  
	  deltaPosCour.x = tab.drx[i];
	  deltaPosCour.y = tab.dry[i];
	  deltaPosCour.z = tab.drz[i];

	  // Dealing with P.B.C in the reference state	  
	  deltaPosInitCris = inverse(hh0) * deltaPosInit;
	  dtr.x = dtr.y = dtr.z = 0.0;
	  if(deltaPosInitCris.x > 0.5) {
	    deltaPosInitCris.x -= 1.0;
	    dtr.x -= 1.0;              
	  } else if(deltaPosInitCris.x < -0.5) {
	    deltaPosInitCris.x += 1.0;
	    dtr.x += 1.0;              
	  }

	  if(deltaPosInitCris.y > 0.5) {
	    deltaPosInitCris.y -= 1.0;
	    dtr.y -= 1.0;              
	  } else if(deltaPosInitCris.y < -0.5) {
	    deltaPosInitCris.y += 1.0;
	    dtr.y += 1.0;              
	  }

	  if(deltaPosInitCris.z > 0.5) {
	    deltaPosInitCris.z -= 1.0;
	    dtr.z -= 1.0;              
	  } else if(deltaPosInitCris.z < -0.5) {
	    deltaPosInitCris.z += 1.0;
	    dtr.z += 1.0;              
	  }
	  deltaPosInit = hh0 * deltaPosInitCris;


	  // Dealing with P.B.C in the current state	  	  
	  deltaPosCourCris = inverse(hht) * deltaPosCour;          
	  deltaPosCourCris.x += dtr.x;
	  deltaPosCourCris.y += dtr.y;
	  deltaPosCourCris.z += dtr.z;
	  
	  if(deltaPosCourCris.x > 0.5) {
	    deltaPosCourCris.x -= 1.0;
	  } else if(deltaPosCourCris.x < -0.5) {
	    deltaPosCourCris.x += 1.0;
	  }

	  if(deltaPosCourCris.y > 0.5) {
	    deltaPosCourCris.y -= 1.0;
	  } else if(deltaPosCourCris.y < -0.5) {
	    deltaPosCourCris.y += 1.0;
	  }

	  if(deltaPosCourCris.z > 0.5) {
	    deltaPosCourCris.z -= 1.0;
	  } else if(deltaPosCourCris.z < -0.5) {
	    deltaPosCourCris.z += 1.0;
	  }
	  deltaPosCour = hht * deltaPosCourCris;

	  // Choosing weight function for neighborhood ponderation
	  double poidsD;
	  double rrInit = norm(deltaPosInit);
	  double rrCour = norm(deltaPosCour);

	  // New version
          poidsD = weight_function(rrInit, m_rrDef);

	  if (compute_static_measures) {
	    if(rrInit <= m_rrDef) {
	      tensorAF += tensor(deltaPosInit, deltaPosInit) * poidsD * 1.0e20;
	      tensorBF += tensor(deltaPosCour, deltaPosInit) * poidsD * 1.0e20;	      
	    }

	    deltaCourInit = deltaPosCour - deltaPosInit;
	    if((norm2(deltaCourInit) >= m_rrDef*m_rrDef/100.) & (rrInit <= m_rrDef)) {	    
	      s += deltaCourInit;
	      slipped_neighbors_count++;
	    }
	  }

	  if (compute_dynamic_measures) {
	    
            if(rrCour <= m_rrDef) {
              deltaVitCour.x = tab.ext.vx[i] - cells[tab.cell][field::vx][tab.part];
              deltaVitCour.y = tab.ext.vy[i] - cells[tab.cell][field::vy][tab.part];
              deltaVitCour.z = tab.ext.vz[i] - cells[tab.cell][field::vz][tab.part];
              tensorAL += tensor(deltaPosCour, deltaPosCour) * poidsD * 1.0e20;
              tensorBL += tensor(deltaVitCour, deltaPosCour) * poidsD * 1.0e20;            
            }
	  }
	  
	}

      if(compute_static_measures) {
	
	invtensorAF = inverse(tensorAF);
	F = AikBkj(tensorBF, invtensorAF);

	// test on nan values
	if (std::isnan(F.m11) == true) F.m11 = 1.0;
	if (std::isnan(F.m12) == true) F.m12 = 0.0;
	if (std::isnan(F.m13) == true) F.m13 = 0.0;

	if (std::isnan(F.m21) == true) F.m21 = 0.0;
	if (std::isnan(F.m22) == true) F.m22 = 1.0;
	if (std::isnan(F.m23) == true) F.m23 = 0.0;

	if (std::isnan(F.m31) == true) F.m31 = 0.0;
	if (std::isnan(F.m32) == true) F.m32 = 0.0;
	if (std::isnan(F.m33) == true) F.m33 = 1.0;        

	// Assigning deformation gradient tensor to per-atom variable F
	local_mechanical_data[tab.cell].F[tab.part] = F;

	// Computing Green-Lagrange strain tensor
	E = 0.5 * ( F * transpose(F) - I);     
	// Assigning Green-Lagrange strain tensor to per-atom variable E      
	local_mechanical_data[tab.cell].E[tab.part] = E;

	// Computing pure rotation and pure stretch tebsir through polar decomposition
	RU_decomposition(F,R,U);
	// Assigning pure rotation and pure stretch tensors to per-atom variables R and U            
	local_mechanical_data[tab.cell].R[tab.part] = R;
	local_mechanical_data[tab.cell].U[tab.part] = U;

	// Computing microrotation vector
	skewR = 0.5 * (R - transpose(R));
	mu.x = 0.5*(skewR.m32-skewR.m23);
	mu.y = 0.5*(skewR.m13-skewR.m31);
	mu.z = 0.5*(skewR.m21-skewR.m12);
	// Assigning microrotation vector to per-atom variable mu      
	local_mechanical_data[tab.cell].mu[tab.part] = mu;

	// Computing slip vector
	if(slipped_neighbors_count > 0)
	  s = s/slipped_neighbors_count;
	// Assigning slip vector to per-atom variable s      
	local_mechanical_data[tab.cell].s[tab.part] = s;

	// Computing local slip orientation tripod
	A = F - I;
	Ltens = A * transpose(A);
	Ntens = transpose(A) * A;

	if((Ltens.m11 > Ltens.m22) && (Ltens.m11 > Ltens.m33)) {
	  lvec.x = sqrt(Ltens.m11);
	  lvec.y = sign(Ltens.m12)*sqrt(Ltens.m22);
	  lvec.z = sign(Ltens.m13)*sqrt(Ltens.m33);
	} else if((Ltens.m22 > Ltens.m11) && (Ltens.m22 > Ltens.m33)) {
	  lvec.x = sign(Ltens.m12)*sqrt(Ltens.m11);
	  lvec.y = sqrt(Ltens.m22);
	  lvec.z = sign(Ltens.m23)*sqrt(Ltens.m33);
	} else {
	  lvec.x = sign(Ltens.m13)*sqrt(Ltens.m11);
	  lvec.y = sign(Ltens.m23)*sqrt(Ltens.m22);
	  lvec.z = sqrt(Ltens.m33);
	}
	double tol=1.e-4;
	if(abs(lvec.z) > tol) {
	  lvec *= sign(lvec.z);
	}else if(abs(lvec.y) > tol) {
	  lvec *= sign(lvec.y);
	}
      
	if((Ntens.m11 > Ntens.m22) && (Ntens.m11 > Ntens.m33)) {
	  nvec.x = sqrt(Ntens.m11);
	  nvec.y = sign(Ntens.m12)*sqrt(Ntens.m22);
	  nvec.z = sign(Ntens.m13)*sqrt(Ntens.m33);
	} else if((Ntens.m22 > Ntens.m11) && (Ntens.m22 > Ntens.m33)) {
	  nvec.x = sign(Ntens.m12)*sqrt(Ntens.m11);
	  nvec.y = sqrt(Ntens.m22);
	  nvec.z = sign(Ntens.m23)*sqrt(Ntens.m33);
	} else {
	  nvec.x = sign(Ntens.m13)*sqrt(Ntens.m11);
	  nvec.y = sign(Ntens.m23)*sqrt(Ntens.m22);
	  nvec.z = sqrt(Ntens.m33);
	}           
	//forcer l'orientation cohérente entre les vecteurs        
	if(abs(nvec.z) > tol) {
	  nvec *= sign(nvec.z);
	} else if(abs(nvec.y) > tol) {
	  nvec *= sign(nvec.y);
	}

	mvec = cross(nvec, lvec);
	// Assigning local slip orientation tripod per-atom variables l,m,n
	local_mechanical_data[tab.cell].l[tab.part] = lvec;
	local_mechanical_data[tab.cell].m[tab.part] = mvec;
	local_mechanical_data[tab.cell].n[tab.part] = nvec;
      }

      if (compute_dynamic_measures) {
	
      	invtensorAL = inverse(tensorAL);
      	L = AikBkj(tensorBL, invtensorAL);

      	// test on nan values
      	if (std::isnan(L.m11) == true) L.m11 = 0.0;
      	if (std::isnan(L.m12) == true) L.m12 = 0.0;
      	if (std::isnan(L.m13) == true) L.m13 = 0.0;

      	if (std::isnan(L.m21) == true) L.m21 = 0.0;
      	if (std::isnan(L.m22) == true) L.m22 = 0.0;
      	if (std::isnan(L.m23) == true) L.m23 = 0.0;

      	if (std::isnan(L.m31) == true) L.m31 = 0.0;
      	if (std::isnan(L.m32) == true) L.m32 = 0.0;
      	if (std::isnan(L.m33) == true) L.m33 = 0.0;        

      	// Assigning deformation gradient tensor to per-atom variable L
	local_mechanical_data[tab.cell].L[tab.part] = L;

      	// Computing vorticity vector
      	skewL = 0.5 * (L - transpose(L));
      	phi.x = 0.5*(skewL.m32-skewL.m23);
      	phi.y = 0.5*(skewL.m13-skewL.m31);
      	phi.z = 0.5*(skewL.m21-skewL.m12);

      	// Assigning vorticity vector to per-atom variable phi
      	local_mechanical_data[tab.cell].phi[tab.part] = phi;
		
      }
    
    }
      
  };

    //Gradient of a vector from reference system
  template< typename GridT/*, typename SecondOrderAnalyzerStruct */>
    struct RefGradientComputeOp
    {
      using CellsT = decltype( GridT{}.cells() );
      const CellsT m_cells_t0;

      GridParticleLocalMechanicalMetrics & m_local_mechanical_data;

      const Mat3d m_xform_t0;
      const Mat3d m_xform;
      const Mat3d m_lattice;
      const double m_rrDef;

      double (*weight_function)(const double&, const double&) = calc_weight;

      //      SecondOrderAnalyzerStruct SecondOrderAnalyzer = false;

      template<class ComputeBufferT, class CellParticlesT /*, class GridParticleLocksT, class ParticleLockT*/>
      inline void operator ()
        (
        // number of neighbors
        size_t n,

        // buffer holding attached data to particle's neighbors
        const ComputeBufferT& tab,
                
        // data and locks accessors for neighbors
        CellParticlesT cells
	//	,GridParticleLocksT locks
	//	,ParticleLockT & particle_lock
        ) const
      {
	GridParticleLocalMechanicalMetrics & local_mechanical_data = m_local_mechanical_data;
	
	Mat3d vector_gradient_tensor;

	//Gradient constructors
	Vec3d deltaVec;
	size_t c_nbh;
	size_t p_nbh;
 
        Vec3d deltaPosInit;
        Vec3d deltaPosInitCris;

	Mat3d vec_tens_ref_position_tensor;
	Mat3d ref_tens_ref_position_tensor;
	Mat3d inv_ref_tens_ref_position_tensor;

        //Initial position
        Vec3d dtr;

        Mat3d hh0 = m_xform_t0 * m_lattice;	
	Vec3d tmp;
	tmp.x = m_cells_t0[tab.cell][field::rx][tab.part];
	tmp.y = m_cells_t0[tab.cell][field::ry][tab.part];
	tmp.z = m_cells_t0[tab.cell][field::rz][tab.part];	
	Vec3d PosInit = m_xform_t0 * tmp;
	
	Vec3d PosInitVois;
	
	// loop over neighborhood
        for(size_t i=0;i<n;i++)
        {
	  
          tmp.x = tab.ext.rx0[i];
          tmp.y = tab.ext.ry0[i];
          tmp.z = tab.ext.rz0[i];	  
	  PosInitVois = m_xform_t0 * tmp;
	  
	  deltaPosInit = PosInitVois - PosInit;

          dtr.x = dtr.y = dtr.z = 0.0;
	  deltaPosInitCris = inverse(hh0) * deltaPosInit;
	  dtr.x = dtr.y = dtr.z = 0.0;
	  if(deltaPosInitCris.x > 0.5) {
	    deltaPosInitCris.x -= 1.0;
	    dtr.x -= 1.0;              
	  } else if(deltaPosInitCris.x < -0.5) {
	    deltaPosInitCris.x += 1.0;
	    dtr.x += 1.0;              
	  }

	  if(deltaPosInitCris.y > 0.5) {
	    deltaPosInitCris.y -= 1.0;
	    dtr.y -= 1.0;              
	  } else if(deltaPosInitCris.y < -0.5) {
	    deltaPosInitCris.y += 1.0;
	    dtr.y += 1.0;              
	  }

	  if(deltaPosInitCris.z > 0.5) {
	    deltaPosInitCris.z -= 1.0;
	    dtr.z -= 1.0;              
	  } else if(deltaPosInitCris.z < -0.5) {
	    deltaPosInitCris.z += 1.0;
	    dtr.z += 1.0;              
	  }
	  deltaPosInit = hh0 * deltaPosInitCris;

          double rrInit = norm(deltaPosInit);
	  
          if(rrInit <= m_rrDef) {
            tab.nbh.get(i, c_nbh, p_nbh);
	    deltaVec = local_mechanical_data[c_nbh].mu[p_nbh] - local_mechanical_data[tab.cell].mu[tab.part];
          }
	  
          double poidsD;
	  
          poidsD = weight_function(rrInit, m_rrDef);
	  
          if(rrInit <= m_rrDef) {
	    
            ref_tens_ref_position_tensor += tensor(deltaPosInit, deltaPosInit) * poidsD * 1.0e20;
            vec_tens_ref_position_tensor += tensor(deltaVec, deltaPosInit) * poidsD * 1.0e20;            
          }       
        }

        inv_ref_tens_ref_position_tensor = inverse(ref_tens_ref_position_tensor);
        
        vector_gradient_tensor = AikBkj(vec_tens_ref_position_tensor, inv_ref_tens_ref_position_tensor);

        save_nan(vector_gradient_tensor);

	local_mechanical_data[tab.cell].vector_gradient_tensor[tab.part] = vector_gradient_tensor;

	Vec3d slip_base_vector1 = local_mechanical_data[tab.cell].l[tab.part];
	Vec3d slip_base_vector2 = local_mechanical_data[tab.cell].m[tab.part];
	Vec3d slip_base_vector3 = local_mechanical_data[tab.cell].n[tab.part];
	
	Mat3d TransferMatrix = make_mat3d(slip_base_vector1, slip_base_vector2, slip_base_vector3);
  
	Mat3d microrotation_gradient_proj_tensor = AikBkj(AikBkj(transpose(TransferMatrix), vector_gradient_tensor), TransferMatrix);
  
	double comp_coin = microrotation_gradient_proj_tensor.m21;
	double comp_vis = microrotation_gradient_proj_tensor.m22;

	local_mechanical_data[tab.cell].dislo_indic[tab.part] = sqrt(comp_coin*comp_coin + comp_vis*comp_vis);

	local_mechanical_data[tab.cell].vis[tab.part] = comp_vis;
	local_mechanical_data[tab.cell].coin[tab.part] = comp_coin;

	local_mechanical_data[tab.cell].dislo_line_ortho[tab.part] = comp_coin*slip_base_vector1 + comp_vis*slip_base_vector2;
	local_mechanical_data[tab.cell].dislo_line[tab.part] = cross(slip_base_vector3, local_mechanical_data[tab.cell].dislo_line_ortho[tab.part]);
      }
    };  

}
