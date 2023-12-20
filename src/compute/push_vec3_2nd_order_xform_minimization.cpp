#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/quantity.h>
#include <exanb/core/domain.h>

#include <onika/soatl/field_pointer_tuple.h>

#include <memory>

/*
      Vec3d d = inv_xform * ( Vec3d{dx[i],dy[i],dz[i]} * dt + Vec3d{d2x[i],d2y[i],d2z[i]} * dt2 );
      x[i] += d.x;
      y[i] += d.y;
      z[i] += d.z;
*/

namespace exaStamp
{
  using namespace exanb;

  // get particle mass from its type. assume the type index is 0 if particle hasn't type field
  template<bool has_atom_type> static inline double get_particle_mass(const double* __restrict__, const uint8_t* __restrict__, size_t);
  template<> inline double get_particle_mass<true>( const double* __restrict__ masses , const uint8_t* __restrict__ types, size_t j)
  {
    return masses[ types[j] ];
  }
  template<> inline double get_particle_mass<false>( const double* __restrict__ masses , const uint8_t* __restrict__ , size_t )
  {
    return masses[0];
  }
  
  template<
    class GridT,
    class Field_X, class Field_Y, class Field_Z,
    class Field_dX, class Field_dY, class Field_dZ,
    class Field_ddX, class Field_ddY, class Field_ddZ,
    class = AssertGridHasFields< GridT, Field_X, Field_Y, Field_Z, Field_dX, Field_dY, Field_dZ , Field_ddX, Field_ddY, Field_ddZ>
    >
  struct PushVec3SecondOrderXFormMIN : public OperatorNode
  {
    static constexpr onika::soatl::FieldId<Field_X>  f_X{};
    static constexpr onika::soatl::FieldId<Field_Y>  f_Y{};
    static constexpr onika::soatl::FieldId<Field_Z>  f_Z{};
    static constexpr onika::soatl::FieldId<Field_ddX> f_ddX{};
    static constexpr onika::soatl::FieldId<Field_ddY> f_ddY{};
    static constexpr onika::soatl::FieldId<Field_ddZ> f_ddZ{};
  
    static constexpr size_t SIMD_VECTOR_SIZE = GridT::CellParticles::ChunkSize ;

    // compile time constant indicating if grid has type field
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;
    
    using PointerTuple = onika::soatl::FieldPointerTuple< GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize , 
                          Field_X, Field_Y, Field_Z,
                          Field_dX, Field_dY, Field_dZ,
                          Field_ddX, Field_ddY, Field_ddZ >;
  
    ADD_SLOT( GridT  , grid       , INPUT_OUTPUT);
    ADD_SLOT( ParticleSpecies, species , INPUT , REQUIRED );    
    ADD_SLOT( double , dt         , INPUT , REQUIRED );
    ADD_SLOT( double , dt_scale   , INPUT , REQUIRED );
    ADD_SLOT( Domain , domain     , INPUT , REQUIRED );
    ADD_SLOT( double , tolF       , INPUT , 1.0e-12 );    

    inline void execute () override final
    {
      const double raw_dt = *(this->dt);
      const double scale = *dt_scale;

      const Mat3d inv_xform = domain->inv_xform();
      const Mat3d xform = domain->xform();
      const double ForceMAX = *(this->tolF);
      
      ldbg<<"PushVec3SecondOrder: dt="<<raw_dt<<", dt_scale="<<scale<<", inv_xform="<<inv_xform<<std::endl;
//      ldbg << "xform*inv_xform="<< inv_xform * domain->xform() << std::endl;

//      GridT& grid2              = *(this->grid);
      ParticleSpecies& species = *(this->species);

      // copy masses to aligned array
      size_t nSpecies = species.size();
      double masses[MAX_PARTICLE_SPECIES];
      for(size_t i=0;i<nSpecies;i++)
      {
        masses[i] = species[i].m_mass;
      }
      
      auto cells = grid->cells();
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();

      const double dt = raw_dt * scale;
      const double dtdemi = dt*0.5;      
      const double dt2 = dt*dt*0.5;

      if( domain->xform_is_identity() )
      {
      
#       pragma omp parallel
        {
          PointerTuple ptrs;   
          GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
          {
            size_t i = grid_ijk_to_index( dims , loc + gl );
            int n = cells[i].size();
            cells[i].capture_pointers( ptrs );
            
            auto* __restrict__ ptr_X = ptrs[ f_X ];
            auto* __restrict__ ptr_Y = ptrs[ f_Y ];
            auto* __restrict__ ptr_Z = ptrs[ f_Z ];
            
            const auto* __restrict__ ptr_ddX = ptrs[ f_ddX ];
            const auto* __restrict__ ptr_ddY = ptrs[ f_ddY ];
            const auto* __restrict__ ptr_ddZ = ptrs[ f_ddZ ];

#           pragma omp simd
            for(int k=0;k<n;k++)
            {
              ptr_X[k] += ptr_ddX[k] * dt2;
              ptr_Y[k] += ptr_ddY[k] * dt2;
              ptr_Z[k] += ptr_ddZ[k] * dt2;
            }
          }
          GRID_OMP_FOR_END
        }
        
      }
      else
      {
#       pragma omp parallel
        {
          PointerTuple ptrs;   
          GRID_OMP_FOR_BEGIN(dims-2*gl,_,loc, schedule(dynamic) )
          {
            size_t i = grid_ijk_to_index( dims , loc + gl );
            int n = cells[i].size();
            cells[i].capture_pointers( ptrs );
            
            auto* __restrict__ ptr_X = ptrs[ f_X ]; 
            auto* __restrict__ ptr_Y = ptrs[ f_Y ];
            auto* __restrict__ ptr_Z = ptrs[ f_Z ];
            
	    auto* __restrict__ ptr_ddX = ptrs[ f_ddX ];
            auto* __restrict__ ptr_ddY = ptrs[ f_ddY ];
            auto* __restrict__ ptr_ddZ = ptrs[ f_ddZ ];

	    const uint8_t* __restrict__ types = cells[i].field_pointer_or_null(field::type);    
	    //double Fmaxloc=0.;
#           pragma omp simd
            for(int k=0;k<n;k++)
            {

	      double masse = get_particle_mass<has_type_field>( masses, types, k );
              Vec3d accel = inv_xform *  Vec3d{ptr_ddX[k],ptr_ddY[k],ptr_ddZ[k]};
              Vec3d forces = accel;
	      
	      double FX = forces.x;
	      double FY = forces.y;
	      double FZ = forces.z;	      

	      //	      lout << "forces before = " << forces << std::endl;
	      
	      if(fabs(FX) > ForceMAX ) {
		if (FX>0.) FX = ForceMAX;
		if (FX<0.) FX = -1.0 * ForceMAX;
		forces.x = FX;
	      }

	      if(fabs(FY) > ForceMAX ) {
		if (FY>0.) FY = ForceMAX;
		if (FY<0.) FY = -1.0 * ForceMAX;
		forces.y = FY;		
	      }	      

	      if(fabs(FZ)> ForceMAX ) {
		if (FZ>0.) FZ = ForceMAX;
		if (FZ<0.) FZ = -1.0 * ForceMAX;
		forces.z = FZ;
	      }
	      
	      //	      lout << "forces after  = " << forces << std::endl;
	      accel = xform * forces;

	      ptr_ddX[k] = accel.x;
	      ptr_ddY[k] = accel.y;
	      ptr_ddZ[k] = accel.z;

	      double MasseI = dtdemi / masse;
              Vec3d v = inv_xform * (accel * MasseI );

              ptr_X[k] += v.x * dt;
              ptr_Y[k] += v.y * dt;
              ptr_Z[k] += v.z * dt;
            }
          }
          GRID_OMP_FOR_END
        }

      }
    }

  //private:
  //  onika::soatl::FieldArrays< GridT::CellParticles::Alignment , GridT::CellParticles::ChunkSize , field::_mass > m_masses;
  };

  template<class GridT> using PushAccelVelocityToPositionXFormMIN = PushVec3SecondOrderXFormMIN<GridT, field::_rx,field::_ry,field::_rz, field::_vx,field::_vy,field::_vz , field::_ax,field::_ay,field::_az >;
  
 // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
   OperatorNodeFactory::instance()->register_factory( "push_f_v_r_xform_minimization", make_grid_variant_operator< PushAccelVelocityToPositionXFormMIN > );
  }

}

