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
#include <onika/physics/units.h>
#include <onika/memory/allocator.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/core/source_term.h>

#include <mpi.h>
#include <iomanip>

#include "MethodSphere.h"
#include <string>
#include <iostream>
#include <tuple>
#include <fstream>
#include "BoxConfiguration.h"
#include "calculateStress.h"
#include "Grid.h"
#include "typedef.h"
#include <regex>

namespace exaStamp
{
  using namespace exanb;

  template<
    class GridT ,
    class = AssertGridHasFields< GridT, field::_vx, field::_vy, field::_vz , field::_fx, field::_fy, field::_fz >
    >
  class StressGridMDStressLab : public OperatorNode
  {
    ADD_SLOT( MPI_Comm       , mpi          , INPUT , MPI_COMM_WORLD );
    
    ADD_SLOT( ParticleSpecies, species      , INPUT , REQUIRED );
    ADD_SLOT( double         , rcut_max     , INPUT_OUTPUT , 0.0 );  // neighborhood distance, in grid space

    ADD_SLOT( long           , grid_subdiv  , INPUT , 3 );
    ADD_SLOT( GridCellValues , grid_cell_values      , INPUT_OUTPUT );

    // Current configuration
    ADD_SLOT( GridT  , grid         , INPUT , REQUIRED );
    ADD_SLOT( Domain , domain       , INPUT , REQUIRED );
    
    // Reference configuration
    ADD_SLOT( GridT  , grid_t0                     , INPUT, OPTIONAL);
    ADD_SLOT( Mat3d  , xform_t0                    , INPUT, OPTIONAL);

    // MDStressLab++ specific input files
    ADD_SLOT( std::string , filename , INPUT, OPTIONAL );
    ADD_SLOT( std::vector<std::string> , modelnames , INPUT, REQUIRED );
    ADD_SLOT( int64_t         , timestep       , INPUT , REQUIRED );    
    ADD_SLOT( std::string     , mesg     , INPUT_OUTPUT , "DONE." );
    
  public:

    // -----------------------------------------------
    // -----------------------------------------------
    inline void execute ()  override final
    {

      if (filename.has_value()) {
        std::cout << "Compute stress field on grid of external file using MDStressLab++" << std::endl;
      
        // compile time constant indicating if grid has type field
        using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
        static constexpr has_type_field_t has_type_field{};

        GridT& grid = *(this->grid);
        
        //	-------------------------------------------------------------------
        // Get number of particles and create box configuration
        int numberOfParticles;
        int referenceAndFinal= true;
        std::string configFileName = *filename;
        std::ifstream file(configFileName);
        if(!file) MY_ERROR("ERROR: output.xyz could not be opened for reading!");
      
        file >> numberOfParticles;
        if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");
      
        // int numberOfParticles = grid.number_of_particles();
        // if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");
w
        BoxConfiguration body{numberOfParticles,referenceAndFinal};
        body.read(configFileName,referenceAndFinal);
        //	-------------------------------------------------------------------

        int nx,ny,nz;
        nx = 60;
        ny = 5;
        nz = 60;  
        Vector3d lowerLimit(0.,0.,0.);
        Vector3d upperLimit(217.8,16.5,198.0);
        ::Grid<Current> gridFromFile(lowerLimit,upperLimit,nx,ny,nz);
        //	-------------------------------------------------------------------

        /*![ComputeStress]*/
        MethodSphere hardy(6,"hardy");
        for (const auto modelname : *modelnames)
          {
            Kim kim(modelname);
            try {
              Stress<MethodSphere,Cauchy> hardyStress(hardy,&gridFromFile);
        
              calculateStress(body, kim,
                              std::tie(),
                              std::tie(hardyStress), true);
              hardyStress.write("project_hardy_" + modelname);
            }
            catch(const std::runtime_error& e){
              std::cout << e.what() << std::endl;
              std::cout << "Compute stress with projected forces failed. Moving on" << std::endl;
            }

            try{
              Stress<MethodSphere,Cauchy> hardyStress(hardy,&gridFromFile);

              // Calculate stress using the process_dedr, if possible
              calculateStress(body, kim,
                              std::tie(),
                              std::tie(hardyStress));
              hardyStress.write_voxel_grid("hardy_" + modelname,nx,ny,nz,lowerLimit,upperLimit);
            }
            catch(const std::runtime_error& e){
              std::cout << e.what() << std::endl;
              std::cout << "Compute stress with process_dedr failed. Moving on" << std::endl;
            }
          }

      } else {

        std::cout << "Compute stress field on internal exaStamp grid using MDStressLab++" << std::endl;

        int rank=0;
        MPI_Comm_rank(*mpi, &rank);

        GridT& grid = *(this->grid);
        GridT& grid_t0 = *(this->grid_t0);
        
        IJK dims = grid.dimension();
        auto cells = grid.cells();
        auto cells_t0 = grid_t0.cells();
        
        int gl = grid.ghost_layers();
        IJK start = {gl,gl,gl};
        IJK end = dims-start;
        IJK local_grid_dims = end-start;
        AABB bnds = grid.grid_bounds();
        AABB bnds_no_ghosts = grid.grid_bounds_no_ghost();
        
        ldbg << "dims = " << dims << std::endl;
        ldbg << "start = " << start << std::endl;
        ldbg << "end = " << end << std::endl;        
        ldbg << "gl = " << gl << std::endl;
        ldbg << "bounds = " << bnds << std::endl;
        ldbg << "bounds = " << bnds_no_ghosts << std::endl;
        
        Mat3d xform = domain->xform();
        Mat3d xform_t0 = *(this->xform_t0);
        Mat3d lattice = diag_matrix( domain->extent() - domain->origin() );
        Mat3d lattice_loc = diag_matrix( bnds_no_ghosts.bmax - bnds_no_ghosts.bmin );
        
        Mat3d Hcur = transpose(xform * lattice_loc);
        Mat3d Href = transpose(xform_t0 * lattice_loc);        

        //	-------------------------------------------------------------------
        // Get number of particles and create box configuration
        int numberOfParticles;
        int numberOfCells;
        int referenceAndFinal= true;
      
        numberOfParticles = grid.number_of_particles();
        numberOfCells = grid.number_of_cells();
        
        if (numberOfParticles < 0) MY_ERROR("Error: Negative number of particles.\n");
        BoxConfiguration body{numberOfParticles,referenceAndFinal};

        body.reference_box(0) = Href.m11;
        body.reference_box(1) = Href.m12;
        body.reference_box(2) = Href.m13;
        body.reference_box(3) = Href.m21;
        body.reference_box(4) = Href.m22;
        body.reference_box(5) = Href.m23;
        body.reference_box(6) = Href.m31;
        body.reference_box(7) = Href.m32;
        body.reference_box(8) = Href.m33;

        body.box(0) = Hcur.m11;
        body.box(1) = Hcur.m12;
        body.box(2) = Hcur.m13;
        body.box(3) = Hcur.m21;
        body.box(4) = Hcur.m22;
        body.box(5) = Hcur.m23;
        body.box(6) = Hcur.m31;
        body.box(7) = Hcur.m32;
        body.box(8) = Hcur.m33;        

        body.pbc(0) = 0;//domain->periodic_boundary_x();
        body.pbc(1) = 0;//domain->periodic_boundary_y();
        body.pbc(2) = 0;//domain->periodic_boundary_z();
        
        ldbg << std::endl;
        ldbg << "Box size = " << std::endl;
        ldbg << body.box << std::endl;
        ldbg << std::endl;
        ldbg << "Reference box size = " << std::endl;
        ldbg << body.reference_box << std::endl;
        ldbg << std::endl;
        ldbg << "Periodic boundary conditions = " << body.pbc << std::endl;

        const auto dom_dims = domain->grid_dimension();
        const auto dom_start = grid.offset();

        ldbg << "dom_dims = " << dom_dims << std::endl;
        ldbg << "dom_start = " << dom_start << std::endl;
        
        body.species.resize(numberOfParticles);
        ldbg << "body.species size = " << body.species.size() << std::endl;

        std::vector<size_t> cell_offsets(numberOfCells + 1, 0);
        for(size_t i = 0; i < numberOfCells; i++) {
          cell_offsets[i + 1] = cell_offsets[i] + cells[i].size();
        }
        
#     pragma omp parallel
        {
          GRID_OMP_FOR_BEGIN(dims,i,cell_loc, schedule(dynamic) )
            {
              GridFieldSetPointerTuple< GridT, FieldSet<field::_rx,field::_ry,field::_rz,field::_vx,field::_vy,field::_vz> > ptrs;

              GridFieldSetPointerTuple< GridT, FieldSet<field::_rx,field::_ry,field::_rz> > ptrs_t0;              
              cells[i].capture_pointers( ptrs );
              cells_t0[i].capture_pointers( ptrs_t0 );
              
              const auto* __restrict__ rx = ptrs[field::rx];
              const auto* __restrict__ ry = ptrs[field::ry];
              const auto* __restrict__ rz = ptrs[field::rz];

              const auto* __restrict__ vx = ptrs[field::vx];
              const auto* __restrict__ vy = ptrs[field::vy];
              const auto* __restrict__ vz = ptrs[field::vz];
              
              const auto* __restrict__ rx0 = ptrs_t0[field::rx];
              const auto* __restrict__ ry0 = ptrs_t0[field::ry];
              const auto* __restrict__ rz0 = ptrs_t0[field::rz];

              const unsigned int n = cells[i].size();
              
              for(unsigned int j=0;j<n;j++)
                {
                  size_t iloc = cell_offsets[i]+j;
                  Vec3d r = { rx[j], ry[j], rz[j] };
                  Vec3d rt = xform * r;                  
                  body.coordinates[Current](iloc,0) = rt.x;
                  body.coordinates[Current](iloc,1) = rt.y;
                  body.coordinates[Current](iloc,2) = rt.z;

                  body.velocities(iloc,0) = vx[j];
                  body.velocities(iloc,1) = vy[j];
                  body.velocities(iloc,2) = vz[j];
                  body.species[iloc] = "Ta";
                }

              const unsigned int n0 = cells_t0[i].size();
              for(unsigned int j=0;j<n;j++)
                {
                  size_t iloc = cell_offsets[i]+j;                  
                  Vec3d r0 = { rx0[j], ry0[j], rz0[j] };
                  Vec3d r0t = xform_t0 * r0;
                  body.coordinates[Reference](iloc,0) = r0t.x;
                  body.coordinates[Reference](iloc,1) = r0t.y;
                  body.coordinates[Reference](iloc,2) = r0t.z;
                }              
            }
          GRID_OMP_FOR_END
            }
        //	-------------------------------------------------------------------

        //	-------------------------------------------------------------------
        // Create MDStressLab grid to get the stress field on
        int nx = local_grid_dims.i;
        int ny = local_grid_dims.j;
        int nz = local_grid_dims.k;

        Vec3d origin_loc = xform * bnds_no_ghosts.bmin;
        Vec3d end_loc = xform * bnds_no_ghosts.bmax;
        
        Vector3d lowerLimit(origin_loc.x,origin_loc.y,origin_loc.z);
        Vector3d upperLimit(end_loc.x,end_loc.y,end_loc.z);
        
        ::Grid<Current> gridFromFile(lowerLimit,upperLimit,nx,ny,nz);
        //	-------------------------------------------------------------------

        MethodSphere hardy(6,"hardy");
        for (const auto modelname : *modelnames)
          {
            Kim kim(modelname);
            try{
              Stress<MethodSphere,Cauchy> hardyStress(hardy,&gridFromFile);

              // Calculate stress using the process_dedr, if possible
              calculateStress(body, kim,
                              std::tie(),
                              std::tie(hardyStress));
              std::string extension = "hardy_" + std::to_string(rank) + "_" + std::to_string(*timestep);
              hardyStress.write_voxel_grid(extension,nx,ny,nz,lowerLimit,upperLimit);
            }
            catch(const std::runtime_error& e){
              lerr << e.what() << std::endl;
              lerr << "Compute stress with process_dedr failed. Moving on" << std::endl;
            }
          }
        
      }
    }

  };

  template<class GridT> using StressGridMDStressLabTmpl = StressGridMDStressLab<GridT>;

  // === register factories ===
  ONIKA_AUTORUN_INIT(compute_stress_grid_mdstresslab)
  {
   OperatorNodeFactory::instance()->register_factory("compute_stress_grid_mdstresslab", make_grid_variant_operator< StressGridMDStressLabTmpl > );
  }

}
