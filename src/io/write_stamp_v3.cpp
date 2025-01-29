#include <onika/math/basic_types_yaml.h>
#include <onika/scg/operator.h>
#include <onika/scg/operator_slot.h>
#include <onika/scg/operator_factory.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/grid.h>
#include <exanb/core/domain.h>
#include <onika/math/basic_types_stream.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <onika/log.h>
#include <onika/physics/units.h>

#include <exaStamp/io/StampV3LegacyIOStructures.h>
#include <exaStamp/compute/unit_system.h>

#include <mpi.h>
#include <fstream>

namespace exaStamp
{
  using namespace exanb;

  template<typename GridT>
  class WriteStampV3Operator : public OperatorNode
  {

    ADD_SLOT(MPI_Comm       , mpi        , INPUT , MPI_COMM_WORLD , DocString{"MPI communicator"} );
    ADD_SLOT(std::string    , filename   , INPUT , REQUIRED       , DocString{"File name"} );
    ADD_SLOT(Domain         , domain     , INPUT , REQUIRED       , DocString{"Domain"} );
    ADD_SLOT(GridT          , grid       , INPUT , REQUIRED       , DocString{"Particle grid"} );
    ADD_SLOT(long           , timestep   , INPUT , REQUIRED               , DocString{"Iteration number"} );
    ADD_SLOT(double         , dt         , INPUT , REQUIRED               , DocString{"Time step"} );
    ADD_SLOT(ParticleSpecies, species    , INPUT , OPTIONAL         , DocString{"Particle species data block"} );
    ADD_SLOT(long           , block_size , INPUT , 1048576                , DocString{"Block size for collective writing"} );
    ADD_SLOT(double         , physical_time, INPUT , REQUIRED , DocString{"Physical time written to file"} );
    ADD_SLOT(bool           , gen_user_input, INPUT , false         , DocString{"print user input to include in order to read file"} );

  public:
    // -----------------------------------------------
    // ----------- Operator documentation ------------
    inline std::string documentation() const override final
    {
      return R"EOF(
        Creates a protection file in StampV3 format:
        this format is fully compatible with the Stamp3 MD code (read/write)
        )EOF";
    }

    inline void execute () override final
    {

      ldbg << "-> entering stampv3 dump routine" << std::endl;

      using has_field_id_t = typename GridT:: template HasField <field::_id>;
      static constexpr bool has_field_id = has_field_id_t();

      using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, field::_vx,field::_vy,field::_vz, field::_fx,field::_fy,field::_fz , field::_virial, field::_id, field::_type>;
      static constexpr double coord_conv = EXASTAMP_CONST_QUANTITY(1.0 * m); // conversion factor between stampv3 unit system and ExaStampV2 internal units
      static constexpr double vel_conv = EXASTAMP_CONST_QUANTITY(1.0 * m/s);
 
      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();
      IJK gstart { gl, gl, gl };
      IJK gend = dims - IJK{ gl, gl, gl };
      auto cells = grid->cells();
      Mat3d xform = domain->xform();

      unsigned long long number_of_particles = 0;
      GRID_BLOCK_FOR_BEGIN(dims,gstart,gend,i,_)
      {
        number_of_particles += cells[i].size();
      }
      GRID_BLOCK_FOR_END

      MPI_Comm comm = *mpi;
      int np=1, rank=0;
      MPI_Comm_size(comm,&np);
      MPI_Comm_rank(comm,&rank);
      
      std::vector<unsigned long long> proc_particles(np,0);
      MPI_Allgather(&number_of_particles,1,MPI_UNSIGNED_LONG_LONG,proc_particles.data(),1,MPI_UNSIGNED_LONG_LONG,comm);
      size_t all_number_of_particles = 0;
      size_t particle_offset = 0;
      for(int i=0;i<np;i++)
      {
        all_number_of_particles += proc_particles[i];
        if( i < rank ) { particle_offset += proc_particles[i]; }
      }
 
      ldbg << "\t-> writing header section" << std::endl;
      // --------- header -------------
      LegacyHeaderIOStruct header;
      header.iterationNumber = *timestep;
      header.time = *physical_time;
      header.xmax = domain->extent().x / coord_conv;
      header.ymax = domain->extent().y / coord_conv;
      header.zmax = domain->extent().z / coord_conv;
      header.xmin = domain->origin().x / coord_conv;
      header.ymin = domain->origin().y / coord_conv;
      header.zmin = domain->origin().z / coord_conv;
      header.particlesTotalNumber = all_number_of_particles;
      header.domainNumber = 0;
      ldbg << "\t\tdomain: ("<<header.xmin<<","<<header.ymin<<","<<header.zmin<<") -> ("<<header.xmax<<","<<header.ymax<<","<<header.zmax<<") : size ("
           <<header.xmax-header.xmin<<","<<header.ymax-header.ymin<<","<<header.zmax-header.zmin<<")"  << std::endl;
 
      bool write_xsv2_extension = false;

      if( !domain->xform_is_identity() && (rank==0 || rank==(np-1)) )
      {
        ldbg <<"\t-> enable XSV2 extension"<<std::endl;
        header.domainNumber = XsV2ExtensionMarker.x;
        write_xsv2_extension = true;
      }

      LegacySystemIOFile dumpFile;
      dumpFile.open(filename->c_str(),"w");
      
      if( rank == 0 ) { dumpFile.writeHeader(header); }
      else { dumpFile.skipHeader(); }
      
      dumpFile.incrementCurrentOffset( particle_offset );
      ldbg << "\t-> header section written" << std::endl;
      ldbg << "\t-> writing particles section" << std::endl;
      
      // --------- particle data ------     
      LegacyParticleIOStruct * particlesArray = new LegacyParticleIOStruct[*block_size];
      unsigned long long pcounter = 0;
      unsigned long long tcounter = 0;
      size_t id_counter = particle_offset + 1;
      
      GRID_BLOCK_FOR_BEGIN(dims,gstart,gend,i,_)
      {
        size_t n = cells[i].size();
        for(size_t j=0;j<n;j++)
        {
          ParticleTupleIO pt = cells[i][j];

          // atom type
          int atom_type = pt[field::type];
          particlesArray[pcounter].particleType = atom_type;

          // position r
          Vec3d r = Vec3d{pt[field::rx],pt[field::ry],pt[field::rz]};
          domain_periodic_location( *domain , r );
          assert( is_inside( domain->bounds() , r ) );

          // r = ( xform * r ) / coord_conv;
          r = r / coord_conv;
          
          particlesArray[pcounter].coordinates[0] = r.x;
          particlesArray[pcounter].coordinates[1] = r.y;
          particlesArray[pcounter].coordinates[2] = r.z;

          // velocity v
          Vec3d v = Vec3d{pt[field::vx],pt[field::vy],pt[field::vz]} / vel_conv;
          particlesArray[pcounter].velocity[0] = v.x;
          particlesArray[pcounter].velocity[1] = v.y;
          particlesArray[pcounter].velocity[2] = v.z;

          // virial diagonal component
          Mat3d virial = pt[field::virial];
          particlesArray[pcounter].johnDoe[0] = virial.m11;
          particlesArray[pcounter].johnDoe[1] = virial.m22;
          particlesArray[pcounter].johnDoe[2] = virial.m33;

          // mass
          double mass = 1.0;
          if( species.has_value() )
          {
            mass = species->at(atom_type).m_mass;
          }
          particlesArray[pcounter].quaternion[0] = mass;
  
          // kinetic energy
          particlesArray[pcounter].quaternion[1] = mass * norm2(v);          

          // velocity norm
          particlesArray[pcounter].quaternion[2] = norm(v);        

          // force norm
          Vec3d f = Vec3d{pt[field::fx],pt[field::fy],pt[field::fz]};
          particlesArray[pcounter].quaternion[3] = norm(f);

          if( has_field_id )
          {
            particlesArray[pcounter].particleID = pt[field::id] + 1;
          }
          else
          {
            particlesArray[pcounter].particleID = id_counter++ ;
          }
          
          ++ pcounter;
          if( pcounter >= size_t(*block_size) )
          {
            dumpFile.writeArrayOfParticles(particlesArray, pcounter);
            tcounter += pcounter;
            pcounter = 0;
          }
        }
      }
      GRID_BLOCK_FOR_END

      if( pcounter > 0 )
      {
        dumpFile.writeArrayOfParticles(particlesArray, pcounter);
        tcounter += pcounter;
        pcounter = 0;
      }
      ldbg << "\t-> particle section written" << std::endl;

      // free write buffer
      delete[]  particlesArray;

      // extension block if needed
      if( write_xsv2_extension && rank==(np-1) )
      {
        dumpFile.writeExtendedBlock( &xform , sizeof(xform) );
        ldbg << "\t-> XSV2 extension : xform = " <<xform << std::endl;
      }
      
      // close file
      dumpFile.close();

      // sum write counts from all processors
      MPI_Allreduce(MPI_IN_PLACE,&tcounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,comm);

      if( tcounter != static_cast<size_t>(header.particlesTotalNumber) )
      {
        lerr << "Wrong number of particles in output: expected "<<header.particlesTotalNumber<<" but "<<tcounter<<" has been written"<<std::endl;
        std::abort();
      }      

      ldbg << "-> leaving stampv3 dump routine" << std::endl;

      if( *gen_user_input )
      {
        lout << std::endl << std::endl << "# include the following in user input file" << std::endl << std::endl << "# Domain properties (some will be overwritten by data reader)" << std::endl;
        exanb::print_user_input( *domain , lout );
        lout << std::endl << "# Species description" << std::endl << "species:"<< std::endl ;
        exaStamp::print_user_input( *species , lout , 1 );
      }

    }

    inline void yaml_initialize(const YAML::Node& node) override final
    {
      YAML::Node tmp;
      if( node.IsScalar() )
      {
        tmp["filename"] = node;
      }
      else { tmp = node; }
      this->OperatorNode::yaml_initialize( tmp );
    }

  };

  // === register factories ===  
  ONIKA_AUTORUN_INIT(write_stamp_v3)
  {
    OperatorNodeFactory::instance()->register_factory( "write_stamp_v3", make_grid_variant_operator< WriteStampV3Operator > );
  }

}
