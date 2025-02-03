#pragma once

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
#include <onika/thread.h>
#include <onika/physics/units.h>

#include <exaStamp/io/StampV3LegacyIOStructures.h>
#include <exaStamp/unit_system.h>
#include <onika/mpi/all_value_equal.h>

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
#include <limits>

// uncomment the folling line to enable read timing prints
// #define STAMP_V3_PROFILING 1

namespace exaStamp
{
  using namespace exanb;

  struct StampV3InputParticle
  {
    Vec3d r , v;
    uint64_t id;
  };

  template<typename GridT>
  static inline void read_stamp_v3(
    MPI_Comm comm,
    const std::string& file_name,
    double enlarge_bounds,
    ReadBoundsSelectionMode bounds_mode,
    GridT& grid,
    Domain& domain, // domain is IN/OUT (.m_cell_size and .m_grid_dims are used as input hints)
    int64_t& iteration_number,
    double& phys_time,
    bool pbc_adjust_xform
    )
  {
    using std::endl;
    using std::string;
    using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz,field::_vx,field::_vy,field::_vz, field::_id, field::_type>;
    using ParticleTuple = decltype( grid.cells()[0][0] );

    // ThC: checks that grid is empty and zero initialized
    grid.reset();

    // read data by block to avoid bugs in MPI/IO ...
    static constexpr size_t io_chunk_size = 1048576;

    // conversion constants
    static constexpr double coord_conv = EXASTAMP_CONST_QUANTITY( 1.0 * m ); // conversion factor between input unity to exaStamp internal unity
    static constexpr double vel_conv = EXASTAMP_CONST_QUANTITY( 1.0 * m/s );

    // MPI Initialization
    int rank=0, np=1;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    // extract name without path
    std::string basename;
    std::string::size_type p = file_name.rfind("/");
    if( p != std::string::npos ) basename = file_name.substr(p+1);
    else basename=file_name;

    // read Stamp V3 file
    LegacyHeaderIOStruct entete ;
    LegacySystemIOFile mpiioDumpFile ;
    lout << "============ "<< basename<<" ============" << endl << std::flush ;

    mpiioDumpFile.open(file_name.c_str());
    mpiioDumpFile.readHeader(entete);

    assert( onika::mpi::all_value_equal(comm,entete) );

    iteration_number = entete.iterationNumber;
    AABB file_bounds = { Vec3d{ entete.xmin , entete.ymin , entete.zmin } * coord_conv , Vec3d{ entete.xmax , entete.ymax , entete.zmax } * coord_conv };
    phys_time = entete.time; //UnityConverterHelper::convert(entete.time, "s"); // physical time from file
    
    size_t n_particles = entete.particlesTotalNumber;
    size_t start = ( n_particles * rank ) / np;
    size_t end = ( n_particles * (rank+1) ) / np;
    size_t count = end - start;

    lout << "Iteration Number = " << iteration_number << endl;
    lout << "Total Particles  = " << n_particles << endl;
    lout << "File bounds      = " << file_bounds << endl;
    lout << "Simulation Time  = " << entete.time << endl;
    lout << "MPI Proc. (PE)   = " << np << endl;
    lout << "Particles / PE   = " << count << endl;
    lout << "Bounds mode      = " << bounds_mode << endl;
    lout << std::flush ;
    
    MPI_Barrier(comm);

    const double dmax = std::numeric_limits<double>::max();
    const double dmin = std::numeric_limits<double>::lowest();
    AABB local_bounds = { Vec3d{ dmax , dmax , dmax } , Vec3d{ dmin , dmin , dmin } };
    long long id_min = std::numeric_limits<long long>::max();
    long long id_max = -1;
    long long nb_zero_r = 0;
    
    StampV3InputParticle* particleBuf = new StampV3InputParticle[count];
    LegacyParticleIOStruct* particlesArray = new LegacyParticleIOStruct[io_chunk_size];
    mpiioDumpFile.incrementCurrentOffset( start );
    
    ldbg << "Read "<<count<<" particles " << std::flush;
    
    bool read_xsv2_extension = ( entete.domainNumber == XsV2ExtensionMarker.x );
    Mat3d xsv2ext_xform = make_identity_matrix();

    size_t read_cursor = 0;
    size_t to_read = count;
    
    while( to_read > 0 )
    {
      ldbg << '.' << std::flush;

      size_t read_n = std::min( to_read , io_chunk_size );
      assert( read_n >= 1 );
      
      std::memset( particlesArray , 0 , read_n * sizeof(LegacyParticleIOStruct) );
      mpiioDumpFile.readArrayOfParticles(particlesArray, read_n);
      
#     pragma omp parallel
      {
        AABB thread_local_bounds = local_bounds;
#       pragma omp for schedule(static) reduction(min:id_min) reduction(max:id_max) reduction(+:nb_zero_r)
        for(size_t i=0;i<read_n;i++)
        {
          Vec3d r = make_vec3d(particlesArray[i].coordinates) * coord_conv ;
          int64_t id = particlesArray[i].particleID;
          if( (r.x==0.0 && r.y==0.0 && r.z==0.0) || std::isinf(r.x) || std::isinf(r.y) || std::isinf(r.z) || std::isnan(r.x) || std::isnan(r.y) || std::isnan(r.z) )
          {
            ldbg << "zero particle #"<<id<<std::endl;
            ++ nb_zero_r;
          }
          id_min = std::min( id_min, static_cast<long long>(id) );
          id_max = std::max( id_max, static_cast<long long>(id) );
          thread_local_bounds = extend( thread_local_bounds , r );

          size_t read_buf_idx = read_cursor + i;
          particleBuf[read_buf_idx].r = r;
          particleBuf[read_buf_idx].v = make_vec3d(particlesArray[i].velocity) * vel_conv;
          assert( id>=0 && id < (1ll<<(63-8)) );
          particleBuf[read_buf_idx].id = id<<8;
          particleBuf[read_buf_idx].id |= (particlesArray[i].particleType & 0xFF);
        }
#       pragma omp critical
        {
          local_bounds = extend( local_bounds , thread_local_bounds );
        }
      }
      
      to_read -= read_n;
      read_cursor += read_n;
    }
    delete [] particlesArray; particlesArray=nullptr;
    ldbg << "Local bounds (1) = " << local_bounds << endl;

    // broadcast transform matrix
    if( read_xsv2_extension )
    {
      static_assert( sizeof(xsv2ext_xform) == 9*sizeof(double) , "Mat3d has different size than exepected" );
      if(rank==(np-1))
      {
        mpiioDumpFile.readExtendedBlock( &xsv2ext_xform, sizeof(xsv2ext_xform) );
      }
      MPI_Bcast((double*) &xsv2ext_xform,9,MPI_DOUBLE,np-1,comm);
      ldbg <<   "Stampv3ext XForm = "<< xsv2ext_xform <<endl;
    }

    mpiioDumpFile.close();

    // if ids are not numbered from 0, shift them
    MPI_Allreduce(MPI_IN_PLACE,&id_min,1,MPI_LONG_LONG,MPI_MIN,comm);
    MPI_Allreduce(MPI_IN_PLACE,&id_max,1,MPI_LONG_LONG,MPI_MAX,comm);
    MPI_Allreduce(MPI_IN_PLACE,&nb_zero_r,1,MPI_LONG_LONG,MPI_SUM,comm);
    ldbg << "Id min           = "<<id_min << endl;
    ldbg << "Id max           = "<<id_max << endl;
    ldbg << "Nb zero particle = "<<nb_zero_r << endl;
    
    // after this points, coordinates in particleArray avec been converted to code length units.
    
    // reduce bounds across all processors (find global bounding box )
    double tmp[6] = { -local_bounds.bmin.x, -local_bounds.bmin.y, -local_bounds.bmin.z , local_bounds.bmax.x, local_bounds.bmax.y, local_bounds.bmax.z };
    MPI_Allreduce(MPI_IN_PLACE,tmp,6,MPI_DOUBLE,MPI_MAX,comm);
    AABB all_bounds = { {-tmp[0],-tmp[1],-tmp[2]} , {tmp[3],tmp[4],tmp[5]} };
    lout << "World bounds     = " << all_bounds << endl;
    lout << "Enlarge bounds   = "<<enlarge_bounds<<endl;

    // adjust domain settings to fit user requirements and file data
    if( pbc_adjust_xform && !domain.xform_is_identity() )
    {
      lerr << "WARNING: stampv3 reader expected XForm to be identity, resetting XForm"<<std::endl;
      domain.set_xform( make_identity_matrix() );
    }
    
    compute_domain_bounds(domain, bounds_mode,enlarge_bounds,file_bounds,all_bounds, pbc_adjust_xform );
    size_t n_domain_cells = grid_cell_count(domain.grid_dimension()) ;
    lout << "Domain XForm     = "<< xsv2ext_xform * domain.xform() <<endl;
    lout << "Domain bounds    = "<<domain.bounds()<<endl;
    lout << "Domain size      = "<<bounds_size(domain.bounds()) <<endl;
    lout << "Real size        = "<<bounds_size(domain.bounds()) * Vec3d{domain.xform().m11,domain.xform().m22,domain.xform().m33} <<endl;
    lout << "Cell size        = "<<domain.cell_size()<<endl;
    lout << "Grid dimensions  = "<<domain.grid_dimension()<<" ("<<n_domain_cells<<" cells)"<< endl;

    // now, we start building the local processor's grid
    grid.set_origin( domain.bounds().bmin );
    grid.set_cell_size( domain.cell_size() ); 

    if( pbc_adjust_xform && !domain.xform_is_identity() )
    {
      Mat3d inv_xform = domain.inv_xform();
      ldbg << "domain xform adjusted with respect to PBC , inv_xform = " << inv_xform << std::endl;
      local_bounds = AABB{ Vec3d{ dmax , dmax , dmax } , Vec3d{ dmin , dmin , dmin } };
      for (size_t i = 0;i<count; i++)
      {
        particleBuf[i].r = inv_xform * particleBuf[i].r;
        local_bounds = extend( local_bounds , particleBuf[i].r );
      }
      ldbg << "local bounds = " << local_bounds << std::endl;  
    }
    
    // take into account matrix from file
    Mat3d defmat = xsv2ext_xform * domain.xform();
    domain.set_xform( defmat );

    // compute local processor's grid size and offset
    assert( local_bounds.bmin.x >= domain.bounds().bmin.x && local_bounds.bmin.y >= domain.bounds().bmin.y && local_bounds.bmin.z >= domain.bounds().bmin.z );
    IJK local_offset = make_ijk( ( local_bounds.bmin - domain.bounds().bmin ) / domain.cell_size() ); // make_ijk( Vec3d ) uses a floor operation on each of x,y and z. this is what we want.
    IJK local_extent = make_ijk( ceil( ( local_bounds.bmax - domain.origin() ) / domain.cell_size() ) );    
    IJK local_dims = local_extent - local_offset;
    AABB agjusted_local_bounds = { local_offset*domain.cell_size() + domain.origin() ,  local_extent*domain.cell_size() + domain.origin() };
    ldbg << "Local offset     = " << local_offset << endl;
    ldbg << "Local dimensions = " << local_dims << endl;
    ldbg << "Local bounds     = " << agjusted_local_bounds << endl;
    assert( is_inside(agjusted_local_bounds,local_bounds) );

    // adapt local processor's grid dimensions and offset to just fit particles it reads    
    grid.set_dimension( local_dims );
    grid.set_offset( local_offset );

//    MPI_Barrier(comm);
//    auto T0 = std::chrono::high_resolution_clock::now();
    size_t n_otb_particles = 0;    
    //spin_mutex_array cell_locks;
    //size_t number_of_cells = grid.number_of_cells();
    //cell_locks.resize( number_of_cells );
    auto cells = grid.cells();

//#   pragma omp parallel for schedule(static)
    for (size_t i = 0;i<count; i++)
    {
      //if( i!=0 && i%(1024*1024)==0 ) std::cout << '.' << std::flush ;
      uint8_t typ = particleBuf[i].id & 0xFF ;
      int64_t id = ( particleBuf[i].id >> 8 ) - id_min;

      double vx = particleBuf[i].v.x;
      double vy = particleBuf[i].v.y;
      double vz = particleBuf[i].v.z;

      Vec3d r = particleBuf[i].r;
      IJK dom_loc = domain_periodic_location( domain, r );

      // find domain's grid location and modify r w.r.t periodic conditions
      IJK loc = dom_loc - local_offset;
      if( grid.contains(loc) )
      {
        assert( loc == grid.locate_cell( r ) );
        size_t cell_index = grid_ijk_to_index(local_dims, loc);
        assert( cell_index < grid.number_of_cells() );
        ParticleTuple t = ParticleTupleIO(r.x,r.y,r.z, vx,vy,vz, id, typ);
        //cell_locks[cell_index].lock();
        cells[cell_index].push_back( t , grid.cell_allocator() );
        //cell_locks[cell_index].unlock();
      }
      else
      {
       ++ n_otb_particles;
      }
    }
    delete [] particleBuf; particleBuf=nullptr;

    grid.rebuild_particle_offsets();

    if( n_otb_particles > 0 )
    {
      lerr << "Warning: " << n_otb_particles << " particles outside local bounds have been ignored" << std::endl;
    }

    lout << "========================================" << endl << endl;
  }

}

