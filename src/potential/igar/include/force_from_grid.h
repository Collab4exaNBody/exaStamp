/*
Licensed to the Apache Software Foundation (ASF) under one
or more contributor license agreements.  See the NOTICE file
distributed with this work for additional information
regarding copyright ownership.  The ASF licenses this file
to you under the Apache License, Version 2.0 (the
"License"); you may not use this file except in compliance
with the License.  You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing,
software distributed under the License is distributed on an
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
KIND, either express or implied.  See the License for the
specific language governing permissions and limitations
under the License.
*/
#pragma once

#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/grid_cell_particles/grid_cell_values.h>
#include <exanb/grid_cell_particles/grid_cell_values_utils.h>
#include <exanb/core/grid_particle_field_accessor.h>
#include <onika/math/quaternion_operators.h>
#include <onika/cuda/cuda.h>

namespace exanb
{

  namespace ParticleCellProjectionTools
  {
    using namespace GridCellValuesUtils;    

//     struct CreateCellValueField
//     {
//       GridCellValues& m_cell_values;
//       std::function<bool(const std::string&)> m_field_selector;
//       const ssize_t m_subdiv = 1;

//       template<class GridT, class FidT>
//       inline void operator () ( GridT& grid , const FidT& proj_field)
//       {
//         using field_type = typename FidT::value_type;
//         if constexpr ( type_support_weighting_v<field_type> )
//         {
//           unsigned int ncomp = 1;
//           if constexpr ( std::is_same_v<field_type,Vec3d> ) ncomp=3;
//           if constexpr ( std::is_same_v<field_type,IJK>   ) ncomp=3;
//           if constexpr ( std::is_same_v<field_type,Mat3d> ) ncomp=9;
//           const auto name = proj_field.short_name();
//           if( m_field_selector(name) )
//           {
//             if( ! m_cell_values.has_field(name) )
//             {
//               ldbg << "Add field "<<name<<std::endl;
//               m_cell_values.add_field(name, m_subdiv , ncomp );
//             }
//           }
//         }
//         else
//         {
//           ldbg << "Field " << proj_field.short_name() << " ignored (not projectable)" << std::endl;
//         }
//       }
//     };

//     struct CollectCellValueFieldToAdd
//     {
//       GridCellValues& m_cell_values;
//       std::vector< AddCellFieldInfo >& m_fields_to_add;
//       std::function<bool(const std::string&)> m_field_selector;
//       const ssize_t m_subdiv = 1;

//       template<class GridT, class FidT>
//       inline void operator () ( GridT& grid , const FidT& proj_field)
//       {
//         using field_type = typename FidT::value_type;
//         if constexpr ( type_support_weighting_v<field_type> )
//         {
//           unsigned int ncomp = 1;
//           if constexpr ( std::is_same_v<field_type,Vec3d> ) ncomp=3;
//           if constexpr ( std::is_same_v<field_type,IJK>   ) ncomp=3;
//           if constexpr ( std::is_same_v<field_type,Mat3d> ) ncomp=9;
//           const auto name = proj_field.short_name();
//           if( m_field_selector(name) )
//           {
//             if( ! m_cell_values.has_field(name) )
//             {
//               ldbg << "Add field "<<name<<std::endl;
//               //m_cell_values.add_field(name, m_subdiv , ncomp );
//               m_fields_to_add.push_back( { name, static_cast<size_t>(m_subdiv), ncomp } );
//             }
//           }
//         }
//         else
//         {
//           ldbg << "Field " << proj_field.short_name() << " ignored (not projectable)" << std::endl;
//         }
//       }
//     };


//     template<class CellsT>
//     struct ProjectCellValueField
//     {
//       const CellsT m_cells_acc;
//       GridCellValues& m_cell_values;
//       std::function<bool(const std::string&)> m_field_selector;
//       const double m_splat_size = 1.0;
//       const ssize_t m_subdiv = 1;
      
//       template<class GridT, class FidT >
//       inline void operator () ( GridT& grid , const FidT& proj_field )
//       {
//         using field_type = typename FidT::value_type;
//         if constexpr ( type_support_weighting_v<field_type> )
//         {
//           using const_ncomp = std::integral_constant<unsigned int
//                 , ( std::is_same_v<field_type,Vec3d> || std::is_same_v<field_type,IJK> ) ? 3
//                 : ( std::is_same_v<field_type,Mat3d> ? 9
//                 : ( std::is_same_v<field_type,Quaternion> ? 4
//                 : ( std::is_arithmetic_v<field_type> ? 1 : 0 ) ) ) >;
//           static constexpr unsigned int ncomp = const_ncomp::value;
          
//           auto cells = grid.cells();
//           const IJK dims = grid.dimension();
//           const ssize_t gl = grid.ghost_layers();
//           const IJK dims_no_gl = dims - 2*gl;
//           //const ssize_t gl = grid.ghost_layers();      
//           const double cell_size = grid.cell_size();
//           const double subcell_size = cell_size / m_subdiv;
//   #       ifndef NDEBUG
//           const ssize_t n_cells = grid.number_of_cells();
//           const ssize_t n_subcells = m_subdiv * m_subdiv * m_subdiv;
//   #       endif

//           if( m_splat_size >= subcell_size )
//           {
//             lerr << "Warning: splat_size ("<<m_splat_size<<") >= subcell_size("<<subcell_size<<"), projected values wil be truncated"<<std::endl;
//           }

//           const auto name = proj_field.short_name();

//           if( m_field_selector(name) )
//           {
//             ldbg << "Project field "<<name<<std::endl;
//             double * __restrict__ ptr = nullptr;
//             assert( size_t(m_subdiv) == m_cell_values.field(name).m_subdiv );
//             assert( size_t(m_subdiv * m_subdiv * m_subdiv * ncomp) == m_cell_values.field(name).m_components );
//             auto accessor = m_cell_values.field_data(name);
//             ptr = accessor.m_data_ptr;
//             const size_t stride = accessor.m_stride;

//             if( ncomp == 0 )
//             {
//               lerr << "Warning: Selected field "<<name<<" cannot be projected to cell values, it will be set to 0" << std::endl;
//             }
            
//             // computes per cell Mass, per cell Mass*Velocity^2, per cell x-component of Momentum and per cell x-component or Virial
//   #         pragma omp parallel
//             {
//   #           pragma omp for schedule(static)
//               for(ssize_t i = 0; i < dims.i*dims.j*dims.k; ++i)
//               {
//                 for(ssize_t j = 0; j < m_subdiv*m_subdiv*m_subdiv; ++j)
//                 {
// 		              size_t index = i * stride + j * ncomp;
// 	                for(unsigned int k=0;k<ncomp;k++) ptr[ index + k ] = 0.0;
//                 }
//               }

//               GRID_OMP_FOR_BEGIN(dims_no_gl,i_no_gl,cell_loc_no_gl, schedule(dynamic) )
//               {
//                 const IJK cell_loc = cell_loc_no_gl + gl;
//                 const size_t i = grid_ijk_to_index( dims , cell_loc );
//                 const Vec3d cell_origin = grid.cell_position( cell_loc );
//                 const auto* __restrict__ rx = cells[i][field::rx];
//                 const auto* __restrict__ ry = cells[i][field::ry];
//                 const auto* __restrict__ rz = cells[i][field::rz];
//                 const unsigned int n = cells[i].size();
//                 for(unsigned int j=0;j<n;j++)
//                 {
//                   Vec3d r = { rx[j] , ry[j] , rz[j] };
//                   IJK center_cell_loc;
//                   IJK center_subcell_loc;
//                   Vec3d rco = r - cell_origin;
//                   localize_subcell( rco, cell_size, subcell_size, m_subdiv, center_cell_loc, center_subcell_loc );
//                   center_cell_loc += cell_loc;

//                   for(int ck=-1;ck<=1;ck++)
//                   for(int cj=-1;cj<=1;cj++)
//                   for(int ci=-1;ci<=1;ci++)
//                   {
//                     IJK nbh_cell_loc;
//                     IJK nbh_subcell_loc;
//                     subcell_neighbor( center_cell_loc, center_subcell_loc, m_subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc );
//                     if( grid.contains(nbh_cell_loc) )
//                     {
//                       ssize_t nbh_cell_i = grid_ijk_to_index( dims , nbh_cell_loc );
//                       ssize_t nbh_subcell_i = grid_ijk_to_index( IJK{m_subdiv,m_subdiv,m_subdiv} , nbh_subcell_loc );
//                       assert( nbh_cell_i>=0 && nbh_cell_i<n_cells );
//                       assert( nbh_subcell_i>=0 && nbh_subcell_i<n_subcells );

//                       // compute weighted contribution of particle to sub cell
//                       Vec3d nbh_cell_origin = grid.cell_position(nbh_cell_loc);
//                       AABB subcell_box = { nbh_cell_origin + nbh_subcell_loc*subcell_size , nbh_cell_origin + (nbh_subcell_loc+1)*subcell_size };
//                       const double w = particle_weight(r, m_splat_size, subcell_box);
//                       using FieldContribType = decltype( field_type{} * w );
//                       [[maybe_unused]] FieldContribType field_contrib;
//                       if constexpr ( ncomp > 0 )
//                       {
//                         field_contrib = m_cells_acc[i][proj_field][j] * w;
//                       }
//                       if constexpr ( ncomp == 1 )
//                       {
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i ] += field_contrib;
//                       }
//                       if constexpr ( ncomp == 3 )
//                       {
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 3 + 0 ] += field_contrib.x;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 3 + 1 ] += field_contrib.y;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 3 + 2 ] += field_contrib.z;
//                       }
//                       if constexpr ( ncomp == 4 )
//                       {
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 4 + 0 ] += field_contrib.w;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 4 + 1 ] += field_contrib.x;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 4 + 2 ] += field_contrib.y;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 4 + 3 ] += field_contrib.z;
//                       }
//                       if constexpr ( ncomp == 9 )
//                       {
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 0 ] += field_contrib.m11;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 1 ] += field_contrib.m12;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 2 ] += field_contrib.m13;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 3 ] += field_contrib.m21;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 4 ] += field_contrib.m22;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 5 ] += field_contrib.m23;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 6 ] += field_contrib.m31;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 7 ] += field_contrib.m32;
//   #                     pragma omp atomic update
//                         ptr[ nbh_cell_i * stride + nbh_subcell_i * 9 + 8 ] += field_contrib.m33;
//                       }
//                     }
//                   }
//                 }
//               }
//               GRID_OMP_FOR_END;
//             }
//           }
//         }// if type support weighting        
//       }
//     };

//     template<class LDBGT, class GridT, class FieldAccessorTupleT>
//     inline void get_particle_force_from_grid(
//         LDBGT& ldbg
//       , GridT& grid
//       , GridCellValues& grid_cell_values
//       , long grid_subdiv
//       , double splat_size
//       , std::function<bool(const std::string&)> field_selector 
//       , const FieldAccessorTupleT & grid_fields )
//     {
    
//       ldbg << "Available fields for particle / grid rojection : ";
//       print_field_tuple(ldbg,grid_fields);
//       ldbg << std::endl;

//       // create cell value fields
//       std::vector<AddCellFieldInfo> fields_to_add;
//       CollectCellValueFieldToAdd collect_fields = {grid_cell_values,fields_to_add,field_selector,grid_subdiv};
//       apply_grid_field_tuple( grid, collect_fields , grid_fields );
//       ldbg << "add "<< fields_to_add.size() << " cell fields :"<<std::endl;
//       for(const auto& f:fields_to_add) ldbg << "\t" << f.m_name<<std::endl;
//       grid_cell_values.add_fields( fields_to_add );

// //      using CellParticleAcessor = GridParticleFieldAccessor<typename GridT::CellParticles * const>;
//       auto cells = grid.cells_accessor(); //{ grid.cells() };

//       ProjectCellValueField<decltype(cells)/*CellParticleAcessor*/> fields_projector = { cells , grid_cell_values,field_selector,splat_size,grid_subdiv};
//       apply_grid_field_tuple( grid, fields_projector , grid_fields );
//     }

    template<class LDBGT, class GridT, class FieldTupleT>
    inline void get_particle_force_from_grid(
        LDBGT& ldbg
      , GridT& grid
      , GridCellValues& grid_cell_values
      , [[maybe_unused]] const FieldTupleT& particle_fields
      , double energy_factor )
    {
      static constexpr const char* energy_field_name = "igar_ep";

      if( !grid_cell_values.has_field(energy_field_name) )
      {
        ldbg << "get_particle_force_from_grid: no '" << energy_field_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }

      ldbg << "get_particle_force_from_grid: subcell-level triquadratic Lagrange interpolation" << std::endl;

      const IJK dims = grid.dimension();
      const ssize_t gl = grid.ghost_layers();
      const IJK dims_no_gl = dims - 2*gl;
      const double cell_size = grid.cell_size();

      auto accessor = grid_cell_values.field_data(energy_field_name);
      const double* __restrict__ energy_ptr = accessor.m_data_ptr;
      const size_t energy_stride = accessor.m_stride;
      const ssize_t subdiv = static_cast<ssize_t>( grid_cell_values.field(energy_field_name).m_subdiv );
      const double subcell_size = cell_size / subdiv;
      const double inv_subcell_size = 1.0 / subcell_size;

      auto cells = grid.cells();

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_gl, _cell_flat_i, cell_loc_no_gl, schedule(static))
        {
          const IJK cell_loc = cell_loc_no_gl + gl;
          const size_t cell_i = grid_ijk_to_index(dims, cell_loc);
          const Vec3d cell_origin = grid.cell_position(cell_loc);
          const unsigned int np = cells[cell_i].size();
          const auto* __restrict__ rx = cells[cell_i][field::rx];
          const auto* __restrict__ ry = cells[cell_i][field::ry];
          const auto* __restrict__ rz = cells[cell_i][field::rz];
          auto* __restrict__ fx = cells[cell_i][field::fx];
          auto* __restrict__ fy = cells[cell_i][field::fy];
          auto* __restrict__ fz = cells[cell_i][field::fz];
          auto* __restrict__ ep = cells[cell_i][field::ep];

          for(unsigned int p = 0; p < np; p++)
          {
            // Position relative to cell corner
            const Vec3d rco = { rx[p] - cell_origin.x, ry[p] - cell_origin.y, rz[p] - cell_origin.z };

            // Locate the subcell the particle falls in
            IJK center_cell_loc, center_subcell_loc;
            localize_subcell(rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc);
            center_cell_loc = center_cell_loc + cell_loc;  // relative → absolute

            // Normalized offset from particle's subcell center: u,v,w ∈ [-0.5, 0.5].
            // Lagrange nodes (neighboring subcell centers) are at -1, 0, +1 in these units.
            const double u = rco.x * inv_subcell_size - (center_subcell_loc.i + 0.5);
            const double v = rco.y * inv_subcell_size - (center_subcell_loc.j + 0.5);
            const double w = rco.z * inv_subcell_size - (center_subcell_loc.k + 0.5);
            ldbg << "u,v,w = (" << u << "," << v << "," << w << ")" << std::endl;
            ldbg << "center_cell_loc = " << center_cell_loc << std::endl;
            // Build 3x3x3 stencil of subcell energy values centered on particle's subcell.
            // subcell_neighbor handles cross-cell-boundary indexing transparently.
            // stencil[ck+1][cj+1][ci+1], center at stencil[1][1][1]. Out-of-domain = 0.
            double stencil[3][3][3] = {};
            for(int ck=-1; ck<=1; ck++)
            for(int cj=-1; cj<=1; cj++)
            for(int ci=-1; ci<=1; ci++)
            {
              IJK nbh_cell_loc, nbh_subcell_loc;
              subcell_neighbor(center_cell_loc, center_subcell_loc, subdiv, IJK{ci,cj,ck}, nbh_cell_loc, nbh_subcell_loc);
              if( grid.contains(nbh_cell_loc) )
              {
                const ssize_t nbh_cell_i    = grid_ijk_to_index(dims, nbh_cell_loc);
                const ssize_t nbh_subcell_i = grid_ijk_to_index(IJK{subdiv,subdiv,subdiv}, nbh_subcell_loc);
                stencil[ck+1][cj+1][ci+1]  = energy_factor * energy_ptr[ nbh_cell_i * energy_stride + nbh_subcell_i ];
              }
            }
            ldbg << "stencil 1 1 1 = " << stencil[1][1][1] << std::endl;
            // // --- DEBUG ---                                                                                                                                                                                                    
            // if( p == 0 )                                                                                                                                                                                                        
            // {
            //   lout << "[IGAR DEBUG] cell_loc=(" << cell_loc.i << "," << cell_loc.j << "," << cell_loc.k << ")"                                                                                                                  
            //       << " subdiv=" << subdiv                                                                                                                                                                                      
            //       << " energy_stride=" << energy_stride                                                                                                                                                                        
            //       << " energy_factor=" << energy_factor                                                                                                                                                                        
            //       << " center_subcell_loc=(" << center_subcell_loc.i << "," << center_subcell_loc.j << "," << center_subcell_loc.k << ")"                                                                                      
            //       << " u=" << u << " v=" << v << " w=" << w                                                                                                                                                                    
            //       << std::endl;
                                                                                                                                                                                                                                
            //   // center stencil entry raw (before energy_factor)
            //   const ssize_t dbg_cell_i    = grid_ijk_to_index(dims, center_cell_loc);                                                                                                                                           
            //   const ssize_t dbg_subcell_i = grid_ijk_to_index(IJK{subdiv,subdiv,subdiv}, center_subcell_loc);                                                                                                                   
            //   lout << "[IGAR DEBUG] raw energy_ptr[" << dbg_cell_i << "*" << energy_stride << "+" << dbg_subcell_i << "]"                                                                                                       
            //       << " = " << energy_ptr[ dbg_cell_i * energy_stride + dbg_subcell_i ]                                                                                                                                         
            //       << "  stencil[1][1][1]=" << stencil[1][1][1]                                                                                                                                                                 
            //       << std::endl;                                                                                                                                                                                                
            // }        
            // // --- END DEBUG --- 

            // Lagrange basis:  L0(t)=t(t-1)/2  L1(t)=1-t^2  L2(t)=t(t+1)/2
            // Derivatives:    dL0=t-1/2        dL1=-2t       dL2=t+1/2
            const double Lu[3]  = { u*(u-1.0)*0.5,  1.0-u*u,  u*(u+1.0)*0.5 };
            const double Lv[3]  = { v*(v-1.0)*0.5,  1.0-v*v,  v*(v+1.0)*0.5 };
            const double Lw[3]  = { w*(w-1.0)*0.5,  1.0-w*w,  w*(w+1.0)*0.5 };
            const double dLu[3] = { u-0.5,  -2.0*u,  u+0.5 };
            const double dLv[3] = { v-0.5,  -2.0*v,  v+0.5 };
            const double dLw[3] = { w-0.5,  -2.0*w,  w+0.5 };

            double ep_interp = 0.0;
            double dEdx = 0.0, dEdy = 0.0, dEdz = 0.0;
            for(int k=0; k<3; k++)
            for(int j=0; j<3; j++)
            for(int i=0; i<3; i++)
            {
              const double e = stencil[k][j][i];
              ep_interp += Lu[i]  * Lv[j]  * Lw[k] * e;
              dEdx      += dLu[i] * Lv[j]  * Lw[k] * e;
              dEdy      += Lu[i]  * dLv[j] * Lw[k] * e;
              dEdz      += Lu[i]  * Lv[j]  * dLw[k] * e;
            }

            ep[p] = ep_interp;
            fx[p] += -dEdx * inv_subcell_size;
            fy[p] += -dEdy * inv_subcell_size;
            fz[p] += -dEdz * inv_subcell_size;
          }
        }
        GRID_OMP_FOR_END;
      }
    }

    template<class LDBGT, class GridT>
    inline void get_particle_force_from_gradient_grid(
        LDBGT& ldbg
      , GridT& grid
      , const GridCellValues& grid_cell_values
      , double force_factor )
    {
      static constexpr const char* dEdx_name = "igar_dEdx";
      static constexpr const char* dEdy_name = "igar_dEdy";
      static constexpr const char* dEdz_name = "igar_dEdz";
      static constexpr const char* e_name = "igar_ep";

      if( !grid_cell_values.has_field(dEdx_name) )
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << dEdx_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }
      if( !grid_cell_values.has_field(dEdy_name) )
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << dEdy_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }
      if( !grid_cell_values.has_field(dEdz_name) )
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << dEdz_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }
      if( !grid_cell_values.has_field(e_name) )
      {
        ldbg << "get_particle_force_from_gradient_grid: no '" << e_name << "' in grid_cell_values, skip" << std::endl;
        return;
      }
      const GridCellField& e_field = grid_cell_values.field(e_name);
      const GridCellField& dEdx_field = grid_cell_values.field(dEdx_name);
      const GridCellField& dEdy_field = grid_cell_values.field(dEdy_name);
      const GridCellField& dEdz_field = grid_cell_values.field(dEdz_name);
      const ssize_t subdiv = static_cast<ssize_t>( dEdx_field.m_subdiv );
      const ssize_t subcell_count = subdiv * subdiv * subdiv;
      assert( subdiv > 0 );
      assert( dEdy_field.m_subdiv == static_cast<size_t>(subdiv) );
      assert( dEdz_field.m_subdiv == static_cast<size_t>(subdiv) );
      assert( dEdx_field.m_components >= static_cast<size_t>(subcell_count) );
      assert( dEdy_field.m_components >= static_cast<size_t>(subcell_count) );
      assert( dEdz_field.m_components >= static_cast<size_t>(subcell_count) );
      assert( dEdx_field.m_components % static_cast<size_t>(subcell_count) == 0 );
      assert( dEdy_field.m_components % static_cast<size_t>(subcell_count) == 0 );
      assert( dEdz_field.m_components % static_cast<size_t>(subcell_count) == 0 );
      const size_t e_comps = e_field.m_components / static_cast<size_t>(subcell_count);
      const size_t dEdx_comps = dEdx_field.m_components / static_cast<size_t>(subcell_count);
      const size_t dEdy_comps = dEdy_field.m_components / static_cast<size_t>(subcell_count);
      const size_t dEdz_comps = dEdz_field.m_components / static_cast<size_t>(subcell_count);

      const IJK dims = grid.dimension();
      const ssize_t gl = grid.ghost_layers();
      const IJK dims_no_gl = dims - 2*gl;
      const double cell_size = grid.cell_size();
      const double subcell_size = cell_size / subdiv;

      auto e_accessor = grid_cell_values.field_data(e_name);
      auto dEdx_accessor = grid_cell_values.field_data(dEdx_name);
      auto dEdy_accessor = grid_cell_values.field_data(dEdy_name);
      auto dEdz_accessor = grid_cell_values.field_data(dEdz_name);
      const double* __restrict__ e_ptr = e_accessor.m_data_ptr;
      const double* __restrict__ dEdx_ptr = dEdx_accessor.m_data_ptr;
      const double* __restrict__ dEdy_ptr = dEdy_accessor.m_data_ptr;
      const double* __restrict__ dEdz_ptr = dEdz_accessor.m_data_ptr;
      const size_t e_stride = e_accessor.m_stride;
      const size_t dEdx_stride = dEdx_accessor.m_stride;
      const size_t dEdy_stride = dEdy_accessor.m_stride;
      const size_t dEdz_stride = dEdz_accessor.m_stride;

      auto cells = grid.cells();

      ldbg << "get_particle_force_from_gradient_grid: apply force from igar_dEdx/igar_dEdy/igar_dEdz"
           << " subdiv=" << subdiv << " force_factor=" << force_factor << std::endl;

#     pragma omp parallel
      {
        GRID_OMP_FOR_BEGIN(dims_no_gl, _cell_flat_i, cell_loc_no_gl, schedule(static))
        {
          const IJK cell_loc = cell_loc_no_gl + gl;
          const size_t cell_i = grid_ijk_to_index(dims, cell_loc);
          const Vec3d cell_origin = grid.cell_position(cell_loc);
          const unsigned int np = cells[cell_i].size();
          const auto* __restrict__ rx = cells[cell_i][field::rx];
          const auto* __restrict__ ry = cells[cell_i][field::ry];
          const auto* __restrict__ rz = cells[cell_i][field::rz];
          auto* __restrict__ e = cells[cell_i][field::ep];
          auto* __restrict__ fx = cells[cell_i][field::fx];
          auto* __restrict__ fy = cells[cell_i][field::fy];
          auto* __restrict__ fz = cells[cell_i][field::fz];

          for(unsigned int p = 0; p < np; p++)
          {
            const Vec3d rco = { rx[p] - cell_origin.x, ry[p] - cell_origin.y, rz[p] - cell_origin.z };

            IJK center_cell_loc, center_subcell_loc;
            localize_subcell(rco, cell_size, subcell_size, subdiv, center_cell_loc, center_subcell_loc);
            center_cell_loc = center_cell_loc + cell_loc;

            if( grid.contains(center_cell_loc) )
            {
              const ssize_t center_cell_i = grid_ijk_to_index(dims, center_cell_loc);
              const ssize_t center_subcell_i = grid_ijk_to_index(IJK{subdiv,subdiv,subdiv}, center_subcell_loc);

              fx[p] += -force_factor * dEdx_ptr[center_cell_i * dEdx_stride + center_subcell_i * dEdx_comps];
              fy[p] += -force_factor * dEdy_ptr[center_cell_i * dEdy_stride + center_subcell_i * dEdy_comps];
              fz[p] += -force_factor * dEdz_ptr[center_cell_i * dEdz_stride + center_subcell_i * dEdz_comps];
              e[p] += e_ptr[center_cell_i * e_stride + center_subcell_i * e_comps];
            }
          }
        }
        GRID_OMP_FOR_END;
      }
    }

  } // ParticleCellProjectionTools

} // ecanb
