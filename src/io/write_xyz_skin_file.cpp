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
#include <exanb/core/quantity.h>
#include <exanb/core/string_utils.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/grid_fields.h>

#include <exanb/defbox/deformation.h>
#include <exanb/defbox/deformation_stream.h>
#include <exanb/defbox/deformation_yaml.h>
#include <exanb/defbox/deformation_math.h>

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
  class WriteXYZSkinfilesOperator : public OperatorNode
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
    
    ADD_SLOT( MPI_Comm           , mpi                 , INPUT );
    ADD_SLOT( GridT              , grid                , INPUT );
    ADD_SLOT( Domain             , domain              , INPUT );
    ADD_SLOT( Deformation        , defbox              , INPUT );    
    ADD_SLOT( long               , timestep            , INPUT);
    ADD_SLOT( bool               , is_ghosts           , INPUT , false);
    ADD_SLOT( ParticleSpecies    , species    , INPUT , REQUIRED );
    ADD_SLOT( GridT              , grid_t0    , INPUT, OPTIONAL );

    // rcut and reference positions for deformation gradient tensor computation
    ADD_SLOT( double             , skin_tol            , INPUT , 0.1);        
    //    ADD_SLOT( PositionLongTermBackup , backup_r_lt , INPUT);

  public:
  
    inline void execute () override final
    {
      namespace fs = std::experimental::filesystem;
      using std::string;
      using std::vector;
      using std::ostringstream;

      GridT& grid = *(this->grid);
      GridT& grid_t0 = *(this->grid_t0);
      
      ParticleSpecies& species = *(this->species);
      bool is_ghosts = *(this->is_ghosts);
      Mat3d xform = domain->xform();

      double skin_low = *skin_tol;
      double skin_high = 1. - (*skin_tol);
      
      size_t n_cells = grid.number_of_cells();
      auto cells = grid.cells();
      auto cells_t0 = grid_t0.cells();      
      
      size_t nb_particles = grid.number_of_particles();
      if(!is_ghosts) { nb_particles -= grid.number_of_ghost_particles(); }
      
      int rank=0, np=1;
      MPI_Comm_rank(*mpi, &rank);
      MPI_Comm_size(*mpi, &np);

      // initialisation : remove and recreate folder where .xyz will be stored
      ldbg << "create xyz_datas_skin dir" << std::endl;
      static bool init = true;
      if(rank==0 && init)
        {
          init = false;
          fs::remove_all("xyz_datas_skin");
          std::error_code ec;
          fs::create_directory("xyz_datas_skin", ec);
        }

      ostringstream filename_xyz;

      // define the structure for the .xyz file name
      filename_xyz << std::string("xyz_datas_skin/atomic_structure_");
      filename_xyz << std::setfill('0') << std::setw(15) << std::to_string(*timestep);
      filename_xyz << std::string(".xyz");

      // we don't want a proc try to write in a folder that doesn't exist
      MPI_Barrier(*mpi);

      // structure for file opening/writing in mpi
      MPI_File mpiFile;
      MPI_Status status;

      // all the processors open the .xyz file
      MPI_File_open(MPI_COMM_WORLD, filename_xyz.str().c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &mpiFile);

      string header, xyz_positions;
      int offset, offset_header, nb_particles_cell;

      // // compute the sum of particle of above processors me included
      // MPI_Scan(&nb_particles, &nb_particles_offset, 1, MPI_INT, MPI_SUM, comm);
      // // retrieve my number of particles to compute the right offset
      // nb_particles_offset = nb_particles_offset - nb_particles;
      nb_particles_cell = 0;

      Mat3d Id33 = make_identity_matrix();      
      Mat3d lattice = diag_matrix(domain->extent()-domain->origin());

      Mat3d lot = transpose(xform * lattice);
      Mat3d lot_t0 = transpose(Id33 * lattice);      

      lout << "lattice : " << lattice << std::endl;            
      lout << "lot     : " << lot << std::endl;
      lout << "lot_t0  : " << lot_t0 << std::endl;      
      // compute the sum of particles within a given skin per proc
      unsigned long nb_particles_skin_per_proc = 0;
      unsigned long nb_particles_skin_per_cell = 0;
      for(size_t c=0; c<n_cells;++c)
      {
        if( !grid.is_ghost_cell(c) || is_ghosts )
        {
          const auto * __restrict__ rx0 = cells_t0[c][field::rx];
          const auto * __restrict__ ry0 = cells_t0[c][field::ry];
          const auto * __restrict__ rz0 = cells_t0[c][field::rz];
          np = cells[c].size();
          for(int pos=0;pos<np;++pos)
          {
            Vec3d pos_vec_t0 = {rx0[pos],ry0[pos],rz0[pos]};
            //      lout << pos_vec_t0 << std::endl;
            pos_vec_t0 = inverse(lot_t0) * pos_vec_t0;
            if ( ((pos_vec_t0.x < skin_low) || (pos_vec_t0.x > skin_high)) || ((pos_vec_t0.y < skin_low) || (pos_vec_t0.y > skin_high)) || ((pos_vec_t0.z < skin_low) || (pos_vec_t0.z > skin_high)) ) {
              nb_particles_skin_per_proc += 1;
            }
            pos_vec_t0 = lot_t0 * pos_vec_t0;       
          }

        }
        if (nb_particles_skin_per_cell != 0)
        {
          nb_particles_skin_per_cell += nb_particles_skin_per_cell;
        }
      }
      // offset of particles for writing correctly in mpi in the xyz file      
      unsigned long nb_particles_skin_offset=0;
      MPI_Scan( &nb_particles_skin_per_proc , &nb_particles_skin_offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, *mpi);
      nb_particles_skin_offset = nb_particles_skin_offset - nb_particles_skin_per_proc;
      
      // count total number of particles among all processors
      unsigned long nb_particles_skin_total=0;
      nb_particles_skin_total = nb_particles_skin_per_proc;
      MPI_Allreduce(MPI_IN_PLACE,&nb_particles_skin_total,1,MPI_UNSIGNED_LONG,MPI_SUM,*mpi);

      // header of xyz file containing the total number of skin particles
      header = format_string("%ld\nLattice=\"%10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e %10.12e\"\n",nb_particles_skin_total, lot.m11, lot.m12, lot.m13, lot.m21, lot.m22, lot.m23, lot.m31, lot.m32, lot.m33);
      offset_header = strlen(header.c_str());
      offset = 0;

      // only processor 0 writes the header
      if (rank==0){
        MPI_File_write( mpiFile, header.data(), header.length() , MPI_CHAR , &status );
      }
      
      // routine d'écriture MPI des particules dans le fichier .xyz
      for(size_t c=0; c<n_cells;++c)
      {
        int np = 0;
        int pos_skin = 0;       
        if( !grid.is_ghost_cell(c) || is_ghosts )
        {
          const auto * __restrict__ rx = cells[c][field::rx];
          const auto * __restrict__ ry = cells[c][field::ry];
          const auto * __restrict__ rz = cells[c][field::rz];

          const auto * __restrict__ rx0 = cells_t0[c][field::rx];
          const auto * __restrict__ ry0 = cells_t0[c][field::ry];
          const auto * __restrict__ rz0 = cells_t0[c][field::rz];
          
          const uint8_t* __restrict__ types = cells[c].field_pointer_or_null(field::type); ONIKA_ASSUME_ALIGNED(types);
          
          np = cells[c].size();
          int pos_cell = 0;               
          for(int pos=0;pos<np;++pos)
          {
            Vec3d pos_vec = {rx[pos],ry[pos],rz[pos]};      
            Vec3d pos_vec_t0 = {rx0[pos],ry0[pos],rz0[pos]};        
            pos_vec_t0 = inverse(lot_t0) * pos_vec_t0;
            bool keep_particle = false;
            if ( ((pos_vec_t0.x < skin_low) || (pos_vec_t0.x > skin_high)) || ((pos_vec_t0.y < skin_low) || (pos_vec_t0.y > skin_high)) || ((pos_vec_t0.z < skin_low) || (pos_vec_t0.z > skin_high)) ) {
              nb_particles_skin_per_proc += 1;
              keep_particle = true;
            }
            pos_vec_t0 = lot_t0 * pos_vec_t0;

            pos_vec = xform * pos_vec;

            const char* type_name = "XX";
            if( types != nullptr ) { type_name = species[types[pos]].m_name; }
            
            if (keep_particle == true) {
              // lout << "offset due to header  : " << offset_header << std::endl;
              // lout << "offset due to proc    : " << nb_particles_skin_offset << std::endl;         
              // lout << "offset due to cell    : " << nb_particles_cell << std::endl;        
              // lout << "offset due to in cell : " << pos_cell << std::endl;         
              xyz_positions = format_string("%-10s \t % .10e \t % .10e \t % .10e \n", type_name, pos_vec.x, pos_vec.y, pos_vec.z);
              offset = offset_header + (pos_cell + nb_particles_cell + nb_particles_skin_offset) * strlen(xyz_positions.c_str());
              MPI_File_write_at( mpiFile, offset, xyz_positions.data() , xyz_positions.length() , MPI_CHAR , &status );
              pos_skin += 1;
              pos_cell += 1;          
            }

          }

        }
        if (pos_skin != 0)
        {
          nb_particles_cell += pos_skin;
        }
      }

      MPI_File_close(&mpiFile);
    }
    
};

  template<class GridT> using WriteXYZSkinfilesOperatorTmpl = WriteXYZSkinfilesOperator<GridT>;
  
  // === register factories ===  
  __attribute__((constructor)) static void register_factories()
  {
    OperatorNodeFactory::instance()->register_factory( "write_xyz_skin_file", make_grid_variant_operator< WriteXYZSkinfilesOperatorTmpl > );
  }

}
