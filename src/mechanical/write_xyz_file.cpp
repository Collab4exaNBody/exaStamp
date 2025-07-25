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

#include <onika/math/basic_types_yaml.h>
#include <onika/math/basic_types.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/parallel_grid_algorithm.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <onika/math/basic_types_stream.h>
#include <exanb/core/domain.h>
#include <onika/string_utils.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/grid_fields.h>

#include <exanb/defbox/deformation.h>
#include <exanb/defbox/deformation_stream.h>
#include <exanb/defbox/deformation_yaml.h>
#include <exanb/defbox/deformation_math.h>
#include <onika/physics/constants.h>

#include <exaStamp/mechanical/cell_particles_local_metrics.h>
#include <exaStamp/mechanical/cell_particles_local_mechanical_metrics.h>
#include <exaStamp/mechanical/cell_particles_local_structural_metrics.h>

#include <onika/soatl/packed_field_arrays.h>
#include <memory>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <string>
#include <iomanip>
#include <experimental/filesystem>
#include <math.h>

namespace exaStamp
{
  
  using namespace exanb;
  
  template< class GridT
/*           ,class = AssertGridHasFields< GridT , field::_rx0, field::_ry0, field::_rz0 > */
           >
  class WriteXYZfilesOperator : public OperatorNode
  {
    //field type
    using has_type_field_t = typename GridT::CellParticles::template HasField < field::_type > ;
    static constexpr bool has_type_field = has_type_field_t::value;

    //field id
    using has_id_field_t = typename GridT::CellParticles::template HasField < field::_id > ;
    static constexpr bool has_id_field = has_id_field_t::value;

    //field id mol
    using has_id_mol_field_t = typename GridT::CellParticles::template HasField < field::_idmol > ;
    static constexpr bool has_id_mol_field = has_id_mol_field_t::value;

    //field rxf
    using has_rxf_field_t = typename GridT::CellParticles::template HasField < field::_rxf > ;
    static constexpr bool has_rxf_field = has_rxf_field_t::value;    
    //field ryf
    using has_ryf_field_t = typename GridT::CellParticles::template HasField < field::_ryf > ;
    static constexpr bool has_ryf_field = has_ryf_field_t::value;    
    //field rzf
    using has_rzf_field_t = typename GridT::CellParticles::template HasField < field::_rzf > ;
    static constexpr bool has_rzf_field = has_rzf_field_t::value;
    
    using VariablesVec = std::vector<std::string>;
        
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT );
    ADD_SLOT( GridT              , grid                , INPUT );
    ADD_SLOT( Domain             , domain              , INPUT );
    ADD_SLOT( Deformation        , defbox              , INPUT );    
    ADD_SLOT( std::string        , filename            , INPUT , REQUIRED               , DocString{"File name"} );
    ADD_SLOT( long               , timestep            , INPUT);
    ADD_SLOT( bool               , is_ghosts           , INPUT , false);
    ADD_SLOT( ParticleSpecies    , species             , INPUT , REQUIRED );

    // Structure that enables the output of per-atom variables computed through the operator compute_local_metrics to be called before writting out atoms if needed
    ADD_SLOT( GridParticleLocalMetrics          , local_data            , INPUT, OPTIONAL);
    ADD_SLOT( GridParticleLocalMechanicalMetrics, local_mechanical_data , INPUT, OPTIONAL);
    ADD_SLOT( GridParticleLocalStructuralMetrics, local_structural_data , INPUT, OPTIONAL);        
    ADD_SLOT( VariablesVec                      , per_atom_data         , INPUT, OPTIONAL);
    ADD_SLOT( bool                              , use_filtered_positions, INPUT , false);

  public:
  
    inline void execute () override final
    {
      namespace fs = std::experimental::filesystem;
      using std::string;
      using std::vector;
      using std::ostringstream;

      static const double conv_energy = 1.e4 * onika::physics::atomicMass / onika::physics::elementaryCharge;	// internal units to eV

      GridT& grid = *(this->grid);
      bool is_ghosts = *(this->is_ghosts);
      Mat3d xform = domain->xform();
      bool use_filtered_positions = (*(this->use_filtered_positions) && has_rxf_field && has_ryf_field && has_rzf_field);
      
      GridParticleLocalMetrics local_metrics;
      GridParticleLocalMechanicalMetrics local_mechanical_metrics;
      GridParticleLocalStructuralMetrics local_structural_metrics;            
      VariablesVec per_atom_variables;
      
      if ( (local_data.has_value() or local_mechanical_data.has_value() or local_structural_data.has_value()) and per_atom_data.has_value() ) {
	if ( local_data.has_value() ) {
	  local_metrics = *(this->local_data);
	}
	if ( local_mechanical_data.has_value() ) {
	  local_mechanical_metrics = *(this->local_mechanical_data);
	}
	if ( local_structural_data.has_value() ) {
	  local_structural_metrics = *(this->local_structural_data);
	}	
      	per_atom_variables = *(this->per_atom_data);
      	std::sort(per_atom_variables.begin(), per_atom_variables.end());	
      }
      
      size_t n_cells = grid.number_of_cells();
      auto cells = grid.cells();
      
      unsigned long nb_particles = grid.number_of_particles();
      if(!is_ghosts) { nb_particles -= grid.number_of_ghost_particles(); }
      
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

//      // initialisation : remove and recreate folder where .xyz will be stored
//      ldbg << "create xyz_datas dir" << std::endl;
//      static bool init = true;
//      if(rank==0 && init)
//        {
//          init = false;
//          fs::remove_all("xyz_datas");
//          std::error_code ec;      
//          fs::create_directory("xyz_datas", ec);
//        }
//
//      ostringstream filename_xyz;
//
//      // define the structure for the .xyz file name
//      filename_xyz << std::string("xyz_datas/atomic_structure_");
//      filename_xyz << std::setfill('0') << std::setw(15) << std::to_string(*timestep);
//      filename_xyz << std::string(".xyz");

      ostringstream filename_xyz;
      filename_xyz << std::string(filename->c_str());
//      filename_xyz << std::string(filename);

      // we don't want a proc try to write in a folder that doesn't exist
      MPI_Barrier(*mpi);

      // structure for file opening/writing in mpi
      MPI_File mpiFile;
      MPI_Status status;

      // all the processors open the .xyz file
      MPI_File_open(MPI_COMM_WORLD, filename_xyz.str().c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiFile);

      string header, xyz_positions;
      unsigned long offset, offset_header, nb_particles_offset, nb_particles_cell;

      // compute the sum of particle of above processors me included
      MPI_Scan(&nb_particles, &nb_particles_offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, *mpi);

      // retrieve my number of particles to compute the right offset
      nb_particles_offset = nb_particles_offset - nb_particles;
      nb_particles_cell = 0;

      // count total number of particles among all processors
      unsigned long total_particle_number = 0;
      for(size_t c=0; c<n_cells;++c)
      {
        int np = 0;
        if( !grid.is_ghost_cell(c) || is_ghosts )
        {
          np = cells[c].size();
        }
        total_particle_number += np;
      }
      MPI_Allreduce(MPI_IN_PLACE,&total_particle_number,1,MPI_UNSIGNED_LONG,MPI_SUM,*mpi);

      //matrix defining the shape of the box. It is written in the header of the xyz file
      Mat3d lattice = diag_matrix(domain->extent()-domain->origin());
      Mat3d lot = transpose(xform * lattice);
      if ( (local_data.has_value() or local_mechanical_data.has_value() or local_structural_data.has_value()) and per_atom_data.has_value() ) {// Writting out atoms type, position as well as per-atom metrics computed through compute_local_metrics operator
      	header = onika::format_string("%ld\nLattice=\"%10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e\" ",total_particle_number, lot.m11, lot.m12, lot.m13, lot.m21, lot.m22, lot.m23, lot.m31, lot.m32, lot.m33);
      	for (size_t var=0; var<per_atom_variables.size(); ++var)
      	  {
      	    header += onika::format_string("\t %s ", per_atom_variables[var]);
      	  }
      	header += onika::format_string("\n");	
      } else {
	header = onika::format_string("%ld\nLattice=\"%10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e\"\n",total_particle_number, lot.m11, lot.m12, lot.m13, lot.m21, lot.m22, lot.m23, lot.m31, lot.m32, lot.m33);
      }

      offset_header = strlen(header.c_str());
      offset = 0;

      // only processor 0 writes the header
      if (rank==0){
      	MPI_File_write(mpiFile, header.data() , header.length() , MPI_CHAR , &status);
      }

      if ( (local_data.has_value() or local_mechanical_data.has_value() or local_structural_data.has_value()) and per_atom_data.has_value() ) {// Writting out atoms type, position as well as per-atom metrics computed through compute_local_metrics operator
      
	// routine d'écriture MPI des particules dans le fichier .xyz
	for(size_t c=0; c<n_cells;++c)
	  {
	    int np = 0;
	    if( !grid.is_ghost_cell(c) || is_ghosts )
	      {
		const auto * __restrict__ rx = cells[c][field::rx];
		const auto * __restrict__ ry = cells[c][field::ry];
		const auto * __restrict__ rz = cells[c][field::rz];
		if (use_filtered_positions) {
		  rx = cells[c][field::rxf];
		  ry = cells[c][field::ryf];
		  rz = cells[c][field::rzf];
		}
		const uint8_t* __restrict__ types = cells[c].field_pointer_or_null(field::type); ONIKA_ASSUME_ALIGNED(types);
	  
		np = cells[c].size();

		for(int pos=0;pos<np;++pos)
		  {
		    Vec3d pos_vec = {rx[pos],ry[pos],rz[pos]};
		    const char* type_name = "XX";
		    int atom_type = 0;
		    if( types != nullptr ) atom_type = types[pos];
		    type_name = species->at(atom_type).m_name;
		    double mass = species->at(atom_type).m_mass;
		    double conv_force = mass*conv_energy;
		    pos_vec = xform * pos_vec;

		    string local_metrics_chain;		
		    local_metrics_chain.clear();
		    // Ce bloc doit peut-être être mis dans une fonction qui renvoie une chaine de caractère contenant toutes les variables à sortir ? 
		    for (size_t var=0; var<per_atom_variables.size(); ++var)
		      {
			if (std::binary_search(per_atom_variables.begin(), per_atom_variables.end(), per_atom_variables[var]) and (per_atom_variables[var] == "pe"))
			  {
			    local_metrics_chain += onika::format_string("\t % .10e", conv_energy*local_metrics[c].pe[pos]);
			  }
			if (std::binary_search(per_atom_variables.begin(), per_atom_variables.end(), per_atom_variables[var]) and (per_atom_variables[var] == "f"))
			  {
			    local_metrics_chain += onika::format_string("\t % .10e \t % .10e \t % .10e", conv_force*local_metrics[c].f[pos].x, conv_force*local_metrics[c].f[pos].y, conv_force*local_metrics[c].f[pos].z);
			  }
			if (std::binary_search(per_atom_variables.begin(), per_atom_variables.end(), per_atom_variables[var]) and (per_atom_variables[var] == "v"))
			  {
			    local_metrics_chain += onika::format_string("\t % .10e \t % .10e \t % .10e", local_metrics[c].v[pos].x, local_metrics[c].v[pos].y, local_metrics[c].v[pos].z);
			  }
			if (std::binary_search(per_atom_variables.begin(), per_atom_variables.end(), per_atom_variables[var]) and (per_atom_variables[var] == "F"))
			  {
			    //			    lout << "Looking for per-atom deformation gradient tensor " << std::endl;
			    local_metrics_chain += onika::format_string("\t % .10e % .10e % .10e % .10e % .10e % .10e % .10e % .10e % .10e",
								 local_mechanical_metrics[c].F[pos].m11,
								 local_mechanical_metrics[c].F[pos].m12,
								 local_mechanical_metrics[c].F[pos].m13,
								 local_mechanical_metrics[c].F[pos].m21,
								 local_mechanical_metrics[c].F[pos].m22,
								 local_mechanical_metrics[c].F[pos].m23,
								 local_mechanical_metrics[c].F[pos].m31,
								 local_mechanical_metrics[c].F[pos].m32,
								 local_mechanical_metrics[c].F[pos].m33);
			  }
			if (std::binary_search(per_atom_variables.begin(), per_atom_variables.end(), per_atom_variables[var]) and (per_atom_variables[var] == "L"))
			  {
			    //			    lout << "Looking for per-atom velocity gradient tensor " << std::endl;
			    local_metrics_chain += onika::format_string("\t % .10e % .10e % .10e % .10e % .10e % .10e % .10e % .10e % .10e",
								 local_mechanical_metrics[c].L[pos].m11,
								 local_mechanical_metrics[c].L[pos].m12,
								 local_mechanical_metrics[c].L[pos].m13,
								 local_mechanical_metrics[c].L[pos].m21,
								 local_mechanical_metrics[c].L[pos].m22,
								 local_mechanical_metrics[c].L[pos].m23,
								 local_mechanical_metrics[c].L[pos].m31,
								 local_mechanical_metrics[c].L[pos].m32,
								 local_mechanical_metrics[c].L[pos].m33);
			  }
			if (std::binary_search(per_atom_variables.begin(), per_atom_variables.end(), per_atom_variables[var]) and (per_atom_variables[var] == "bispectrum"))
			  {
			    //			    lout << "size bs : " << local_structural_metrics[c].bispectrum[pos].size() << std::endl;
			    for (size_t indbi=0; indbi < local_structural_metrics[c].bispectrum[pos].size(); indbi++)
			      {
				local_metrics_chain += onika::format_string(" % .10e",local_structural_metrics[c].bispectrum[pos][indbi]);
			      }
			  }
		      }
				  
		    xyz_positions = onika::format_string("%-10s \t % .10e \t % .10e \t % .10e", type_name, pos_vec.x, pos_vec.y, pos_vec.z) + local_metrics_chain + onika::format_string(" \n");
		
		    offset = offset_header + (pos + nb_particles_cell + nb_particles_offset) * strlen(xyz_positions.c_str());
		    MPI_File_write_at( mpiFile, offset, xyz_positions.data(), xyz_positions.length() , MPI_CHAR , &status );
		  }

	      }
	    if (np != 0)
	      {
		nb_particles_cell += np;
	      }
	  }

      } else {
      
	// routine d'écriture MPI des particules dans le fichier .xyz
	for(size_t c=0; c<n_cells;++c)
	  {
	    int np = 0;
	    if( !grid.is_ghost_cell(c) || is_ghosts )
	      {
		const auto * __restrict__ rx = cells[c][field::rx];
		const auto * __restrict__ ry = cells[c][field::ry];
		const auto * __restrict__ rz = cells[c][field::rz];
		if (use_filtered_positions) {
		  rx = cells[c][field::rxf];
		  ry = cells[c][field::ryf];
		  rz = cells[c][field::rzf];
		}		
		const uint8_t* __restrict__ types = cells[c].field_pointer_or_null(field::type); ONIKA_ASSUME_ALIGNED(types);
	  
		np = cells[c].size();

		for(int pos=0;pos<np;++pos)
		  {
		    Vec3d pos_vec = {rx[pos],ry[pos],rz[pos]};
		    const char* type_name = "XX";
		    int atom_type = 0;
		    if( types != nullptr ) atom_type = types[pos];
		    type_name = species->at(atom_type).m_name;
		    pos_vec = xform * pos_vec;
		    xyz_positions = onika::format_string("%-10s \t % .10e \t % .10e \t % .10e \n", type_name, pos_vec.x, pos_vec.y, pos_vec.z);
		    offset = offset_header + (pos + nb_particles_cell + nb_particles_offset) * strlen(xyz_positions.c_str());
		    MPI_File_write_at( mpiFile, offset, xyz_positions.data(), xyz_positions.length() , MPI_CHAR , &status );
		  }

	      }
	    if (np != 0)
	      {
		nb_particles_cell += np;
	      }
	  }
	
      }
      
      MPI_File_close(&mpiFile);
    }
    
};

  template<class GridT> using WriteXYZfilesOperatorTmpl = WriteXYZfilesOperator<GridT>;
  
  // === register factories ===  
  ONIKA_AUTORUN_INIT(write_xyz_file)
  {
    OperatorNodeFactory::instance()->register_factory( "write_xyz_file", make_grid_variant_operator< WriteXYZfilesOperatorTmpl > );
  }

}
