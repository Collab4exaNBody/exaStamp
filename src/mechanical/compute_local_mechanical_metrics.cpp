#include <exanb/core/basic_types_yaml.h>
#include <exanb/core/basic_types.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/basic_types_stream.h>
#include <exanb/core/log.h>
#include <exanb/core/domain.h>

#include <exanb/fields.h>

#include <exanb/defbox/deformation.h>
#include <exanb/defbox/deformation_stream.h>
#include <exanb/defbox/deformation_yaml.h>
#include <exanb/defbox/deformation_math.h>

#include <exaStamp/mechanical/cell_particles_local_mechanical_metrics.h>
#include <exaStamp/mechanical/compute_local_mechanical_metrics.h>

#include <exanb/core/config.h>
#include <exanb/particle_neighbors/chunk_neighbors.h>
#include <exanb/compute/compute_cell_particle_pairs.h>

#include <cmath>

namespace exaStamp
{
  using namespace exanb;
  
  template< class GridT >
  class ComputeLocalMechanicalMetricsOperator : public OperatorNode
  {
    //field id
    using has_id_field_t = typename GridT::CellParticles::template HasField < field::_id > ;
    static constexpr bool has_id_field = has_id_field_t::value;

    //field type
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    //field rxf
    using has_rxf_field_t = typename GridT::CellParticles::template HasField < field::_rxf > ;
    static constexpr bool has_rxf_field = has_rxf_field_t::value;    

    //field ryf
    using has_ryf_field_t = typename GridT::CellParticles::template HasField < field::_ryf > ;
    static constexpr bool has_ryf_field = has_ryf_field_t::value;    

    //field rzf
    using has_rzf_field_t = typename GridT::CellParticles::template HasField < field::_rzf > ;
    static constexpr bool has_rzf_field = has_rzf_field_t::value;

    //field vx
    using has_vx_field_t = typename GridT::CellParticles::template HasField < field::_vx > ;
    static constexpr bool has_vx_field = has_vx_field_t::value;

    //field vy
    using has_vy_field_t = typename GridT::CellParticles::template HasField < field::_vy > ;
    static constexpr bool has_vy_field = has_vy_field_t::value;

    //field vz
    using has_vz_field_t = typename GridT::CellParticles::template HasField < field::_vz > ;
    static constexpr bool has_vz_field = has_vz_field_t::value;
    
    using VariablesVec = std::vector<std::string>;
    
    ADD_SLOT( GridT                             , grid                  , INPUT );
    ADD_SLOT( Domain                            , domain                , INPUT );
    ADD_SLOT( Deformation                       , defbox                , INPUT );    
    ADD_SLOT( long                              , timestep              , INPUT );
    ADD_SLOT( bool                              , is_ghosts             , INPUT , false);
    ADD_SLOT( GridParticleLocalMechanicalMetrics, local_mechanical_data , OUTPUT );
    
    ADD_SLOT( bool                              , use_filtered_positions      , INPUT , false);
    ADD_SLOT( bool                              , compute_static_measures     , INPUT , true);
    ADD_SLOT( bool                              , compute_dynamic_measures    , INPUT , false);
    ADD_SLOT( bool                              , perform_dislocation_analysis, INPUT , false);    
    ADD_SLOT( std::string                       , weight_function             , INPUT , "calc_weight");    
    ADD_SLOT( GridT                             , grid_t0                     , INPUT, OPTIONAL);
    ADD_SLOT( Mat3d                             , xform_t0                    , INPUT, OPTIONAL);
    ADD_SLOT( double                            , rcut_max                    , INPUT , 0.0);    
    ADD_SLOT( double                            , rcut_localdef               , INPUT , REQUIRED);
    
    // Neighbors
    ADD_SLOT( exanb::GridChunkNeighbors                , chunk_neighbors             , INPUT , exanb::GridChunkNeighbors{} ); // New neighbor data structure

  public:
  
    inline void execute () override final
    {
      assert( chunk_neighbors->number_of_cells() == grid->number_of_cells() );
      bool has_chunk_neighbors = chunk_neighbors.has_value();
      if( !has_chunk_neighbors )
      {
        lerr << "No neighbors input data available" << std::endl;
        std::abort();
      }
      using std::vector;
      
      GridT& grid = *(this->grid);
      //      bool is_ghosts = *(this->is_ghosts);
      
      size_t n_cells = grid.number_of_cells();
      auto cells = grid.cells();      

      GridParticleLocalMechanicalMetrics& local_mechanical_data = *(this->local_mechanical_data);
      local_mechanical_data.resize( n_cells );

      if( ! grid_t0.has_value() )
        {
          lerr << "Lagrangian measures requested, but no grid_t0 given as input" << std::endl;
          std::abort();
        }

      bool use_filtered_positions = (*(this->use_filtered_positions) && has_rxf_field && has_ryf_field && has_rzf_field);
      bool compute_static_measures = *(this->compute_static_measures);
      bool compute_dynamic_measures = (*(this->compute_static_measures) && has_vx_field && has_vy_field && has_vz_field);
      bool perform_dislocation_analysis = *(this->perform_dislocation_analysis);
      
      double (*weight_function_mode)(const double&, const double&);
      if(*(this->weight_function) == "wendland") {weight_function_mode = Wendland_weight;} 
      else {weight_function_mode = calc_weight;}
      
      Mat3d xform = domain->xform();
      Mat3d xform_t0 = *(this->xform_t0);      
      double rrDef = *(this->rcut_localdef);
      Mat3d lattice = diag_matrix(domain->extent()-domain->origin());
      
      LinearXForm cp_xform { xform };
      ComputePairOptionalLocks<false> cp_locks {};

      if (compute_static_measures) {
	for(size_t c=0; c<n_cells;++c)
	  {
	    local_mechanical_data[c].F.resize( cells[c].size() );
	    local_mechanical_data[c].E.resize( cells[c].size() );	  
	    local_mechanical_data[c].R.resize( cells[c].size() );
	    local_mechanical_data[c].U.resize( cells[c].size() );
	    local_mechanical_data[c].mu.resize( cells[c].size() );
	    local_mechanical_data[c].s.resize( cells[c].size() );
	    local_mechanical_data[c].l.resize( cells[c].size() );
	    local_mechanical_data[c].m.resize( cells[c].size() );
	    local_mechanical_data[c].n.resize( cells[c].size() );
	  }
      }

      if (compute_dynamic_measures) {
	for(size_t c=0; c<n_cells;++c)
	  {
	    local_mechanical_data[c].L.resize( cells[c].size() );
	    local_mechanical_data[c].phi.resize( cells[c].size() );	  

	  }
      }

      if (perform_dislocation_analysis) {
	for(size_t c=0; c<n_cells;++c)
	  {
	    local_mechanical_data[c].vector_gradient_tensor.resize( cells[c].size() );
	    local_mechanical_data[c].dislo_indic.resize( cells[c].size() );
	    local_mechanical_data[c].vis.resize( cells[c].size() );
	    local_mechanical_data[c].coin.resize( cells[c].size() );
	    local_mechanical_data[c].dislo_line.resize( cells[c].size() );
	    local_mechanical_data[c].dislo_line_ortho.resize( cells[c].size() );
	  }
      }
      
      using CPBufT = ComputePairBuffer2<false,true,R0VExtraStorage,CopyParticleR0V>;
      std::function<void(CPBufT&)> cp_init_func = [this]( CPBufT & cpbuf )
      	{
      	  cpbuf.process_neighbor.m_cells_t0 = grid_t0->cells();
      	};
      auto cp_force_buf = make_compute_pair_buffer( cp_init_func );

      DeformationGradientComputeOp<GridT> deformation_gradient_compute_op { grid_t0->cells() , local_mechanical_data, compute_static_measures, compute_dynamic_measures, xform_t0, xform, lattice, rrDef, weight_function_mode};
	
      exanb::GridChunkNeighborsLightWeightIt<false> nbh_it{ *chunk_neighbors };
      auto optional = make_compute_pair_optional_args( nbh_it, ComputePairNullWeightIterator{} , cp_xform, cp_locks );
      
      if (use_filtered_positions)
      {	
	      PosititionFields< field::_rxf , field::_ryf, field::_rzf> FilteredPositionFields;
	      compute_cell_particle_pairs( grid, *rcut_max , false, optional, cp_force_buf, deformation_gradient_compute_op , FieldSet<>{}, parallel_execution_context() , FilteredPositionFields);
	      if (perform_dislocation_analysis) 
	      {
	        RefGradientComputeOp<GridT> dislocation_measure_compute_op {grid_t0->cells(), local_mechanical_data, xform_t0, xform, lattice, rrDef, weight_function_mode};	  
	        compute_cell_particle_pairs( grid, *rcut_max , false, optional, cp_force_buf, dislocation_measure_compute_op, FieldSet<>{}, parallel_execution_context() , FilteredPositionFields);	  
	      }	
      }
      else
      {
	      compute_cell_particle_pairs( grid, *rcut_max , false, optional, cp_force_buf, deformation_gradient_compute_op , FieldSet<>{}, parallel_execution_context() );
	      if (perform_dislocation_analysis)
	      {
	        RefGradientComputeOp<GridT> dislocation_measure_compute_op {grid_t0->cells(), local_mechanical_data, xform_t0, xform, lattice, rrDef, weight_function_mode};	  
	        compute_cell_particle_pairs( grid, *rcut_max , false, optional, cp_force_buf, dislocation_measure_compute_op, FieldSet<>{}, parallel_execution_context() );	  
	      }
      }         
    }

    struct alignas(onika::memory::DEFAULT_ALIGNMENT) R0VExtraStorage
    {
      // temporary arrays
      alignas(onika::memory::DEFAULT_ALIGNMENT) double rx0[exanb::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double ry0[exanb::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double rz0[exanb::MAX_PARTICLE_NEIGHBORS];
      
      alignas(onika::memory::DEFAULT_ALIGNMENT) double vx[exanb::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double vy[exanb::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double vz[exanb::MAX_PARTICLE_NEIGHBORS];
    };

    // functor that populate compute buffer's extended storage for particle charges
    struct CopyParticleR0V
    {
      using CellsT = decltype( GridT{}.cells() );
      CellsT m_cells_t0 = nullptr;

      template<class ComputeBufferT, class FieldArraysT, class NbhDataT >
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, const FieldArraysT * cells, size_t cell_b, size_t p_b, const NbhDataT& nbh_data) const noexcept
      {
        const auto i = tab.count;
        
        tab.ext.rx0[i] = m_cells_t0[cell_b][field::rx][p_b];
        tab.ext.ry0[i] = m_cells_t0[cell_b][field::ry][p_b];
        tab.ext.rz0[i] = m_cells_t0[cell_b][field::rz][p_b];
      
        tab.ext.vx[i] = cells[cell_b][field::vx][p_b];
        tab.ext.vy[i] = cells[cell_b][field::vy][p_b];
        tab.ext.vz[i] = cells[cell_b][field::vz][p_b];

        DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, nbh_data );
      }
    };

    // private:
    struct alignas(onika::memory::DEFAULT_ALIGNMENT) R0ExtraStorage
    {
      // temporary arrays
      alignas(onika::memory::DEFAULT_ALIGNMENT) double rx0[exanb::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double ry0[exanb::MAX_PARTICLE_NEIGHBORS];
      alignas(onika::memory::DEFAULT_ALIGNMENT) double rz0[exanb::MAX_PARTICLE_NEIGHBORS];
    };
  
    // functor that populate compute buffer's extended storage for particle charges
    struct CopyParticleR0
    {
      using CellsT = decltype( ((GridT*)nullptr)->cells() );
      CellsT m_cells_t0 = nullptr;

      template<class ComputeBufferT, class FieldArraysT, class NbhDataT >
      ONIKA_HOST_DEVICE_FUNC inline void operator () (ComputeBufferT& tab, const Vec3d& dr, double d2, const FieldArraysT * cells, size_t cell_b, size_t p_b, const NbhDataT& nbh_data) const noexcept
      {
        const auto i = tab.count;
        
        tab.ext.rx0[i] = m_cells_t0[cell_b][field::rx][p_b];
        tab.ext.ry0[i] = m_cells_t0[cell_b][field::ry][p_b];
        tab.ext.rz0[i] = m_cells_t0[cell_b][field::rz][p_b];
 
        DefaultComputePairBufferAppendFunc{} ( tab, dr, d2, cells, cell_b, p_b, nbh_data );
     }
   };

    //POST-TREATMENT FOR GRAD PHI PROJECTION ON SLIP PLANES AND CARACTERISATION OF DISLOCATION
    struct DislocationAnalyzer
    {
      GridParticleLocalMechanicalMetrics & m_local_mechanical_data;
  
      inline void operator ()
        ( 
          const Mat3d& microrotation_gradient_tensor,
        
          size_t cell_index,
          size_t part_index
        )
      {
        GridParticleLocalMechanicalMetrics & local_mechanical_data = m_local_mechanical_data;

        Vec3d slip_base_vector1 = m_local_mechanical_data[cell_index].l[part_index];
        Vec3d slip_base_vector2 = m_local_mechanical_data[cell_index].m[part_index];
        Vec3d slip_base_vector3 = m_local_mechanical_data[cell_index].n[part_index];

        Mat3d TransferMatrix = make_mat3d(slip_base_vector1, slip_base_vector2, slip_base_vector3);

        Mat3d microrotation_gradient_proj_tensor = AikBkj(AikBkj(transpose(TransferMatrix), microrotation_gradient_tensor), TransferMatrix);

        double comp_coin = microrotation_gradient_proj_tensor.m21;
        double comp_vis = microrotation_gradient_proj_tensor.m22;

        local_mechanical_data[cell_index].dislo_indic[part_index] = sqrt(comp_coin*comp_coin + comp_vis*comp_vis);

        local_mechanical_data[cell_index].vis[part_index] = comp_vis;
        local_mechanical_data[cell_index].coin[part_index] = comp_coin;

        local_mechanical_data[cell_index].dislo_line_ortho[part_index] = comp_coin*slip_base_vector1 + comp_vis*slip_base_vector2;
        local_mechanical_data[cell_index].dislo_line[part_index] = cross(slip_base_vector3, local_mechanical_data[cell_index].dislo_line_ortho[part_index]);
      }
    };
   
};

  template<class GridT> using ComputeLocalMechanicalMetricsOperatorTmpl = ComputeLocalMechanicalMetricsOperator<GridT>;
  
  // === register factories ===  
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "compute_local_mechanical_metrics", make_grid_variant_operator< ComputeLocalMechanicalMetricsOperatorTmpl > );
  }

}
