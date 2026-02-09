/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements. See the NOTICE file
distributed with this work for additional information
regarding copyright ownership. The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License. You may obtain a copy of the License at
  http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied. See the License for the
specific language governing permissions and limitations
under the License.
*/

#include <memory>

#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid_fields.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/physics/units.h>
#include <onika/physics/constants.h>
#include <exaStamp/unit_system.h>

#include <onika/memory/allocator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/source_term.h>

#include <mpi.h>
#include <iomanip>

#include "DeepPot.h"
#include <exaStamp/potential/mlips/deepmd/deepmd.h>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_ep, field::_fx, field::_fy, field::_fz >
    >
  class DeepMDForce : public OperatorNode
  {
    ADD_SLOT( MPI_Comm       , mpi          , INPUT , MPI_COMM_WORLD );
    ADD_SLOT( ParticleSpecies, species      , INPUT , REQUIRED );
    ADD_SLOT( double         , rcut_max     , INPUT_OUTPUT , 0.0 );
    ADD_SLOT( GridT  , grid         , INPUT , REQUIRED );
    ADD_SLOT( Domain , domain       , INPUT , REQUIRED );
    //    ADD_SLOT( std::string , model   , INPUT , REQUIRED );
    //    ADD_SLOT( std::vector<std::string> , coefs   , INPUT , REQUIRED );
    ADD_SLOT( long           , grid_subdiv  , INPUT , 3 );
    ADD_SLOT( GridCellValues , grid_cell_values      , INPUT_OUTPUT );
    ADD_SLOT( deepmd::DeepPot , deep_pot     , INPUT );
    
  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {

      *rcut_max = (*deep_pot).cutoff();
      
      if( grid->number_of_cells() == 0 ) { return; }

      lout << "DeepMD force computation" << std::endl;

      static constexpr double econv = EXASTAMP_CONST_QUANTITY( 1. * eV );
      //      double econv = ONIKA_QUANTITY( 1.0 * eV ).convert();
      static constexpr double fconv = EXASTAMP_CONST_QUANTITY( 1. * eV/ang );
      lout << "fconv = " << fconv << std::endl;
      static constexpr double vconv = EXASTAMP_CONST_QUANTITY( 1. * eV );
      
      int rank=0;
      MPI_Comm_rank(*mpi, &rank);

      GridT& grid = *(this->grid);
      // compile time constant indicating if grid has type field
      using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
      static constexpr has_type_field_t has_type_field{};
        
      IJK dims = grid.dimension();
      auto cells = grid.cells();
        
      int gl = grid.ghost_layers();
      IJK start = {gl,gl,gl};
      IJK end = dims-start;
      IJK local_grid_dims = end-start;
      AABB bnds = grid.grid_bounds();
      AABB bnds_no_ghosts = grid.grid_bounds_no_ghost();
        
      Mat3d xform = domain->xform();
      Mat3d lattice = diag_matrix( domain->extent() - domain->origin() );
      Mat3d lattice_loc = diag_matrix( bnds_no_ghosts.bmax - bnds_no_ghosts.bmin );
      Mat3d Hcur = transpose(xform * lattice_loc);

      int numberOfParticles = grid.number_of_particles();
      int numberOfCells = grid.number_of_cells();
      
      if (numberOfParticles < 0) lerr << "Error: Negative number of particles.\n" << std::endl;
      
      //      deepmd::DeepPot dp (*model);
      
      // std::vector<double > coord = {1., 0., 0., 0., 0., 1.5, 1. ,0. ,3.};
      // std::vector<double > cell = {10., 0., 0., 0., 10., 0., 0., 0., 10.};
      // std::vector<int > atype = {1, 0, 1};
      // double e;
      // std::vector<double > f, v;
      // dp.compute (e, f, v, coord, atype, cell);

      // Prepare your inputs
      // 3N coordinates [x1,y1,z1,x2,y2,z2,...]      
      std::vector<double> coord;
      coord.resize(3*numberOfParticles);
      // N atom types [0,1,0,1,...]      
      std::vector<int> atype;
      atype.resize(numberOfParticles);
      // 9 elements for 3x3 cell matrix
      std::vector<double> cell;
      cell.resize(9);
      cell[0] = Hcur.m11;
      cell[1] = Hcur.m12;
      cell[2] = Hcur.m13;      
      cell[3] = Hcur.m21;
      cell[4] = Hcur.m22;
      cell[5] = Hcur.m23;      
      cell[6] = Hcur.m31;
      cell[7] = Hcur.m32;
      cell[8] = Hcur.m33;
      // Output variables
      double energy;
      std::vector<double> force;
      std::vector<double> virial;

      std::vector<size_t> cell_offsets(numberOfCells + 1, 0);
      for(size_t i = 0; i < numberOfCells; i++) {
        cell_offsets[i + 1] = cell_offsets[i] + cells[i].size();
      }

      //	-------------------------------------------------------------------
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,i,cell_loc, schedule(dynamic) )
          {
            GridFieldSetPointerTuple< GridT, FieldSet<field::_rx,field::_ry,field::_rz,field::_type> > ptrs;
            cells[i].capture_pointers( ptrs );
            
            const auto* __restrict__ rx = ptrs[field::rx];
            const auto* __restrict__ ry = ptrs[field::ry];
            const auto* __restrict__ rz = ptrs[field::rz];
            const auto* __restrict__ pt = ptrs[field::type];
              
            const unsigned int n = cells[i].size();
              
            for(unsigned int j=0;j<n;j++)
              {
                size_t iloc = cell_offsets[i]+j;
                Vec3d r = { rx[j], ry[j], rz[j] };
                Vec3d rt = xform * r;                  
                coord[iloc + 0] = rt.x;
                coord[iloc + 1] = rt.y;
                coord[iloc + 2] = rt.z;
                atype[iloc] = pt[j];
              }
          }
        GRID_OMP_FOR_END
          }
      //	-------------------------------------------------------------------
      
      // Compute energy and forces
      std::string type_map_str;
      (*deep_pot).compute(energy, force, virial, coord, atype, cell);
      lout << "force 0 = " << force[0] << std::endl;
      //	-------------------------------------------------------------------
#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims,i,cell_loc, schedule(dynamic) )
          {
            GridFieldSetPointerTuple< GridT, FieldSet<field::_fx,field::_fy,field::_fz> > ptrs;
            cells[i].capture_pointers( ptrs );
            
            auto* fx = ptrs[field::fx];
            auto* fy = ptrs[field::fy];
            auto* fz = ptrs[field::fz];
              
            const unsigned int n = cells[i].size();
              
            for(unsigned int j=0;j<n;j++)
              {
                size_t iloc = cell_offsets[i]+j;
                fx[j] += force[iloc + 0] * fconv;
                fy[j] += force[iloc + 1] * fconv;
                fz[j] += force[iloc + 2] * fconv;
                //                ep += energy[iloc];
              }
          }
        GRID_OMP_FOR_END
          }
      //	-------------------------------------------------------------------
      
    }

  };

  template<class GridT> using DeepMDForceTmpl = DeepMDForce<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(deepmd_force)
  {
   OperatorNodeFactory::instance()->register_factory("deepmd_force", make_grid_variant_operator< DeepMDForceTmpl > );
  }

}
