#pragma once

#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/grid_fields.h>
#include <onika/math/basic_types_stream.h>
#include <onika/log.h>
#include <onika/physics/units.h>
#include <exaStamp/unit_system.h>
#include <exanb/core/particle_id_codec.h>

#include <exanb/io/mpi_file_io.h>
#include <exaStamp/molecule/stampv4_io.h>
#include <onika/mpi/all_value_equal.h>
#include <onika/oarray.h>

#include <exaStamp/molecule/molecule_species.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
//#include <chrono>
#include <string>
#include <limits>
#include <memory>
#include <map>

// uncomment the folling line to enable read timing prints
//#define STAMP_V4_VERBOSE_DEBUG 1

namespace exaStamp
{
  using namespace exanb;

  template<class LDBG, class ARRAY>
  static inline void set_local_bound(LDBG& ldbg, size_t count, const ARRAY& array, AABB& local_bounds, const Vec3d pos_scaling = {1.,1.,1.} )
  {

    if( count > 0 )
      {
        Vec3d p0 = {array[0].Position[0], array[0].Position[1], array[0].Position[2]};
        p0 = p0 * pos_scaling;

        local_bounds = AABB{ p0 , p0 };
      }

#   pragma omp parallel
    {
      AABB thread_local_bounds = local_bounds;
#     pragma omp for schedule(static)
      for (size_t i = 1;i<count; i++)
        {
          Vec3d p_i = {array[i].Position[0], array[i].Position[1], array[i].Position[2]};
          p_i = p_i * pos_scaling;
          thread_local_bounds = extend( thread_local_bounds , p_i );
        }
#     pragma omp critical
      {
        local_bounds = extend( local_bounds , thread_local_bounds );
      }
    }//end parallel

    ldbg << "Local bounds     = " << local_bounds << std::endl;

  }

  template< class _Version = std::integral_constant<unsigned int,41> , class _UseIOMolRigid42 = std::integral_constant<bool,(_Version::value==42)> >
  struct StampV4Options
  {
    static inline constexpr unsigned int Version = _Version::value;
    static inline constexpr bool UseIOMolRigid42 = _UseIOMolRigid42::value;
  };

  template<class LDBG, class GridT, class StampV4OptionsT = StampV4Options<> >
  static inline void read_stamp_v4(
                       exanb::MpiIO& file,
                       LDBG& ldbg,
                       MPI_Comm comm,
                       const std::string& file_name,
                       double enlarge_bounds,
                       ReadBoundsSelectionMode bounds_mode,
                       GridT& grid,
                       Domain& domain,
                       int64_t& iteration_number,
                       double& dt,
                       double& phystime,
                       ParticleSpecies& species,
                       bool enable_xform_scale,
                       bool pbc_adjust_xform,
                       StampV4OptionsT options = {}
                       )
  {
    ldbg << "Begin read_stamp_v4  " << std::endl << std::flush ;
    //using namespace stampv4;
    using IndexT = std::conditional_t< options.Version==41 , int32_t , int64_t >;
    using stampv4::IOVersion;
    //using IOGraines = stampv4::IOGraines;
    using IOMolecules = stampv4::IOMolecules;
    using IOEntete = stampv4::IOEntete<IndexT>;
    using IOAtomes = stampv4::IOAtomes<IndexT>;
    using IOMolRig = std::conditional_t< options.UseIOMolRigid42 , stampv4::IOMolRigV4_2 , stampv4::IOMolRigV4_1 >;

    // MPI Initialization
    int rank=0, np=1;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &np);

    //==========================================================================================
    //==========================================================================================
    //==========================================================================================
    std::string basename;
    std::string::size_type p = file_name.rfind("/");
    if( p != std::string::npos ) basename = file_name.substr(p+1);
    else basename=file_name;
    lout << "============ "<< basename <<" ============" << std::endl ;

    ldbg << "\t--> Start reading the _entete_ of the MPIIO file - offset = " << file.current_offset() << "\n" << std::endl ;
    // ------------------------------
    // --------- entete -------------
    std::unique_ptr<IOEntete> entete = std::make_unique<IOEntete>();
    file.read( entete.get() );
    file.increment_offset( entete.get() );
    assert( exanb::all_value_equal(comm, *entete.get()) );

    //AABB file_bounds = { Vec3d{ entete.xmin , entete.ymin , entete.zmin } * coord_conv , Vec3d{ entete.xmax , entete.ymax , entete.zmax } * coord_conv };

    lout << "Format            = Stamp_v4\n";
    lout << "Version           = " << options.Version << "\n";
    lout << "Number of atoms   = " << entete->Natomes << "\n";
    lout << "Cell length x     = " << std::setw(14) << entete->long_a                << "\n";
    lout << "Cell length y     = " << std::setw(14) << entete->long_b                << "\n";
    lout << "Cell length z     = " << std::setw(14) << entete->long_c                << "\n";
    lout << "Cell angle x      = " << std::setw(14) << entete->angle_a               << "\n";
    lout << "Cell angle y      = " << std::setw(14) << entete->angle_b               << "\n";
    lout << "Cell angle z      = " << std::setw(14) << entete->angle_g               << "\n";
    lout << "Matrix            = " << std::setw(14) << entete->MatriceCR[0][0] << " " 
                                   << std::setw(14) << entete->MatriceCR[0][1] << " " 
                                   << std::setw(14) << entete->MatriceCR[0][2] << " ; "
                                   << std::setw(14) << entete->MatriceCR[1][0] << " " 
                                   << std::setw(14) << entete->MatriceCR[1][1] << " " 
                                   << std::setw(14) << entete->MatriceCR[1][2] << " ; "
                                   << std::setw(14) << entete->MatriceCR[2][0] << " " 
                                   << std::setw(14) << entete->MatriceCR[2][1] << " " 
                                   << std::setw(14) << entete->MatriceCR[2][2] << "\n";
    lout << "Iteration number  = " << entete->NumeroIterationAbsolu << "\n";
    lout << "Physical time     = " << std::setw(14) << entete->tempsPhysique         << "\n";
    ldbg << "current time step = " << std::setw(14) << entete->dt_adaptatif          << "\n";
    ldbg << "\t\t*****\n";
    ldbg << "\t\tenergie totale      = " << std::setw(14) << entete->EnergieTotale         << "\n";
    ldbg << "\t\tenergie potentielle = " << std::setw(14) << entete->EnergiePotentielle    << "\n";
    ldbg << "\t\tenergie cinetique   = " << std::setw(14) << entete->EnergieCinetique      << "\n";
    ldbg << "\t\t*****\n";
    ldbg << "\t\tbloc molecule       = " << std::setw( 6) << entete->bloc_molecules        << "\n";
    ldbg << "\t\tbloc DPD            = " << std::setw( 6) << entete->bloc_dpd              << "\n";
    ldbg << "\t\tbloc graines        = " << std::setw( 6) << entete->bloc_graines          << "\n";
    ldbg << "\t\tbloc molrig         = " << std::setw( 6) << entete->bloc_molrig           << "\n";
    ldbg << "\t\tbloc monomeres      = " << std::setw( 6) << entete->bloc_monomeres        << "\n";
    ldbg << "\t\tbloc polymerisation = " << std::setw( 6) << entete->bloc_polymerisation   << "\n";
    ldbg << "\t\tbloc posfiltre      = " << std::setw( 6) << entete->bloc_posfiltre        << "\n";
    ldbg << std::endl;

    size_t n_particles           = entete->Natomes;
    size_t particle_offset_start = ( n_particles * rank ) / np;
    size_t particle_offset_end   = ( n_particles * (rank+1) ) / np;
    size_t count                 = particle_offset_end-particle_offset_start;

    // count /= 256;

    lout << "Number of PE      = " << np           << "\n";
    lout << "Particles / PE    = " << count        << "\n";
    lout << "Bounds mode       = " << bounds_mode  << "\n";
    ldbg << std::endl;


    // conversion of the cell parameters to the deformation matrix
    // https://en.wikipedia.org/wiki/Fractional_coordinates
    double alpha = EXASTAMP_QUANTITY( entete->angle_a * degree );
    double beta  = EXASTAMP_QUANTITY( entete->angle_b * degree );
    double gamma = EXASTAMP_QUANTITY( entete->angle_g * degree );

    double ma = EXASTAMP_QUANTITY( entete->long_a * m );
    double mb = EXASTAMP_QUANTITY( entete->long_b * m );
    double mc = EXASTAMP_QUANTITY( entete->long_c * m );

    Vec3d a{ma, 0, 0};
    Vec3d b{mb * std::cos(gamma), mb * std::sin(gamma), 0};

    double cx = mc * std::cos(beta);
    double cy = mc * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) /std::sin(gamma);
    double cz = sqrt(mc * mc - cx * cx - cy * cy);

    Vec3d c{cx, cy, cz};
    
    bool scaling_matrix = false;
    Vec3d position_scaling { 1., 1., 1. };
    if( entete->angle_a==90. && entete->angle_b==90. && entete->angle_g==90. )
    {
      ldbg << "\t--> Orthogonal matrix" << std::endl;
      a = Vec3d{ ma , 0. , 0. };
      b = Vec3d{ 0. , mb , 0. };
      c = Vec3d{ 0. , 0. , mc };
      scaling_matrix = true;
    }

    Mat3d xform = make_identity_matrix();

    // FIXME: XSv2 extension not compatible with StampV4 format changes
    if( entete->has_xsv2_extension() )
    {
      xform.m11 = entete->MatriceCR[0][0];
      xform.m21 = entete->MatriceCR[1][0];
      xform.m31 = entete->MatriceCR[2][0];
      xform.m12 = entete->MatriceCR[0][1];
      xform.m22 = entete->MatriceCR[1][1];
      xform.m32 = entete->MatriceCR[2][1];
      xform.m13 = entete->MatriceCR[0][2];
      xform.m23 = entete->MatriceCR[1][2];
      xform.m33 = entete->MatriceCR[2][2];
      ldbg << "\t--> extended xsv2 mode, matrix = " << xform << "\n" << std::endl;
      //scaling_matrix = is_diagonal_matrix( xform );
    }
    else 
    {
      xform.m11 = entete->MatriceCR[0][0];
      xform.m21 = entete->MatriceCR[1][0];
      xform.m31 = entete->MatriceCR[2][0];
      xform.m12 = entete->MatriceCR[0][1];
      xform.m22 = entete->MatriceCR[1][1];
      xform.m32 = entete->MatriceCR[2][1];
      xform.m13 = entete->MatriceCR[0][2];
      xform.m23 = entete->MatriceCR[1][2];
      xform.m33 = entete->MatriceCR[2][2];
      xform = xform * make_mat3d(a,b,c);
      ldbg << "\t--> original StampV4 matrix = " << xform << "\n" << std::endl;
    }

    // get general informations about the dynamics
    iteration_number = entete->NumeroIterationAbsolu ;
    phystime         = EXASTAMP_QUANTITY( entete->tempsPhysique * s );
    dt               = EXASTAMP_QUANTITY( entete->dt_adaptatif * s );

    MPI_Barrier(comm);
    // --------- entete -------------
    // ------------------------------
    ldbg << "\t--> End reading the _entete_ of the MPIIO file - offset = " << file.current_offset() << "\n" << std::endl;

    //==========================================================================================
    //==========================================================================================
    //==========================================================================================

    //init data for local bound
    const double dmax = std::numeric_limits<double>::max();
    const double dmin = std::numeric_limits<double>::lowest();
    AABB local_bounds = { Vec3d{ dmax , dmax , dmax } , Vec3d{ dmin , dmin , dmin } };

    ldbg << "\t--> Start reading atoms from the MPIIO file - offset = " << file.current_offset() << "\n" << std::endl;
    
    // ------------------------------
    // --------- atomes -------------
    ldbg << "\t--> Allocate read buffer for "<<count<<" atoms" << std::endl;
    std::unique_ptr<IOAtomes[]> atomes = std::make_unique<IOAtomes[]>(count);
    file.increment_offset( atomes.get(), particle_offset_start );           //displace the offset to the beginning of the _atomes_ block assigned to this process
    file.read( atomes.get(), count);
    file.increment_offset( atomes.get(), n_particles - particle_offset_start ); //displace the offset to the end of the global _atomes_ block

    // safety, avoid non null terminated strings
    const std::string REPLACEMOLNAME="????";
    //const char* TARGETMOLNAME="H2O";
    for (size_t i = 0;i<count; i++)
    {
      atomes[i].typeAtome[15] = '\0';
      //if( REPLACEMOLNAME == atomes[i].typeAtome ) strncpy(atomes[i].typeAtome,TARGETMOLNAME,16);
      //atomes[i].typeAtome[15] = '\0';
    }

    ldbg << "\t--> lecture types ok"<< std::endl;

    // init domain bounds    
    if( entete->has_xsv2_extension() )
    {
      ldbg << "\t--> Extension XSV2" << std::endl;
      AABB bounds = { {entete->bloc_Rjohndoe01,entete->bloc_Rjohndoe02,entete->bloc_Rjohndoe03} , {entete->bloc_Rjohndoe04,entete->bloc_Rjohndoe05,entete->bloc_Rjohndoe06} };
      double cell_size = entete->bloc_Rjohndoe07;
      domain.set_cell_size( cell_size );
      domain.set_bounds( bounds );
    }
    else 
    {
      // In stamp, reduce coordinates are between -0.5 and 0.5
      Vec3d min_dom {- 0.5, - 0.5, - 0.5};
      Vec3d max_dom {  0.5,   0.5,   0.5};
      if( scaling_matrix && enable_xform_scale && is_diagonal(xform) )
      {
        position_scaling = diagonal(xform);
        min_dom = min_dom * position_scaling;
        max_dom = max_dom * position_scaling;
        xform = make_identity_matrix();
        ldbg << "\t--> Scale positions, reset matrix to identity, scaling="<<position_scaling << std::endl;
      }
      else
      {
        ldbg << "\t--> keep original positions, stored matrix = "<< xform << std::endl;

        // if domain is not roughly square, and we use matrix from file, adjust it so that cells are as cubic shaped as possible
        position_scaling = Vec3d { norm(column1(xform)) , norm(column2(xform)) , norm(column3(xform)) };
        double minar = std::min( std::min( position_scaling.x , position_scaling.y ) , position_scaling.z );
        if( minar <= 1.e-16 )
        {
          lerr << "Bad matrix aspect ratio, too distorted : minar="<<minar<<" ar="<<position_scaling <<" , xform="<<xform<<std::endl;
          std::abort();
        }
        //dom_aspect_ratio = dom_aspect_ratio / minar;
        double maxar = std::max( std::max( position_scaling.x , position_scaling.y ) , position_scaling.z );
        ldbg << "\t--> dom aspect ratio = "<<  position_scaling << std::endl;
        if( maxar/minar > 1.5 )
        {
          Vec3d cell_shape = position_scaling / domain.grid_dimension();
          double cell_size = std::min( std::min( cell_shape.x , cell_shape.y ) , cell_shape.z );
          IJK gdims = make_ijk( position_scaling / cell_size );
          //domain.set_cell_size( cell_size );
          domain.set_grid_dimension( gdims );
          min_dom = min_dom * position_scaling;
          max_dom = max_dom * position_scaling;
          xform = xform * diag_matrix( reciprocal(position_scaling) );
          ldbg << "\t--> AR adjust : gdims="<< gdims << std::endl;
          ldbg << "\t--> AR adjust : xform="<< xform << std::endl;      
          ldbg << "\t--> AR adjust : min_dom="<< min_dom << std::endl;      
          ldbg << "\t--> AR adjust : max_dom="<< max_dom << std::endl;
          if( ! pbc_adjust_xform )
          {
            ldbg << "\t--> AR adjust : pbc_adjust_xform forced ON" << std::endl;
            pbc_adjust_xform = true;
          }
        }
        else
        {
          position_scaling = Vec3d { 1. , 1. , 1. };
        }
      }
      domain.set_bounds( AABB{ min_dom, max_dom} );
    }
      
    ldbg << "\t--> initial bounds = " << domain.bounds() << std::endl;

    // define local_bound in mpi proc
    set_local_bound(ldbg, count, atomes, local_bounds, position_scaling);
    Vec3d position_offset = { 0. , 0. , 0. };
    if( ! is_inside( domain.bounds() , local_bounds ) )
    {
      position_offset = ( domain.bounds().bmin - local_bounds.bmin ) + ( ( bounds_size(domain.bounds()) - bounds_size(local_bounds) ) * 0.5 );
      ldbg << "automatic position offset = "<< position_offset<<std::endl;
      local_bounds = AABB { local_bounds.bmin + position_offset , local_bounds.bmax + position_offset };
    }
    if( ! is_inside( domain.bounds() , local_bounds ) )
    {
      lerr << "Local particle bounds outside domain bounds\n"
           << "local bounds = "<<local_bounds<<" , size="<<bounds_size(local_bounds)<< "\n"
           << "domain bounds = "<<domain.bounds()<<" , size="<<bounds_size(domain.bounds())<<std::endl<<std::flush;
      std::abort();
    }
    ldbg << "\t--> Particle bounds = " << local_bounds << std::endl;
    ldbg << "\t--> cell_size = " << domain.cell_size() << std::endl;

    // compute domain size and enlarge it as requested 
    domain.set_xform( make_identity_matrix() );
    ldbg << "\t*** compute_domain_bounds"<< std::endl;
    exanb::compute_domain_bounds(domain, ReadBoundsSelectionMode::DOMAIN_BOUNDS, enlarge_bounds, AABB{}, AABB{}, pbc_adjust_xform);    
    ldbg << "\t--> pbc_adjust_xform = "<<std::boolalpha<<pbc_adjust_xform << std::endl;
    ldbg << "\t--> bounds = "<<domain.bounds() << std::endl;
    ldbg << "\t--> generated xform = "<<domain.xform() << std::endl;
    Mat3d position_inv_xform = domain.inv_xform();
    xform = domain.xform() * xform;
    domain.set_xform( xform );
    ldbg << "\t--> final xform = "<< domain.xform() << std::endl;

    // now, we start building the local processor's grid
    lout << "Cell size         = " << domain.cell_size() << std::endl;
    lout << "Domain bounds     = " << domain.bounds() << std::endl;
    grid.set_origin( domain.bounds().bmin );
    grid.set_cell_size( domain.cell_size() );

    // compute local processor's grid size and offset
    //  -> make_ijk( Vec3d ) uses a floor operation on each of x,y and z. this is what we want.
    IJK local_offset = make_ijk( ( local_bounds.bmin - domain.bounds().bmin ) / domain.cell_size() ); 
    IJK local_extent = make_ijk( ceil( ( local_bounds.bmax - domain.bounds().bmin ) / domain.cell_size() ) );
    IJK local_dims = local_extent - local_offset;
    AABB adjusted_local_bounds = { local_offset*domain.cell_size() + domain.bounds().bmin ,  local_extent*domain.cell_size() + domain.bounds().bmin };
    ldbg << "Local dimensions = " << local_dims << std::endl;
    ldbg << "Local offset = " << local_offset << std::endl;
    ldbg << "adj. loc. bounds = " << adjusted_local_bounds << std::endl;
    assert( is_inside(adjusted_local_bounds,local_bounds) );

    // adapt local processor's grid dimensions and offset to just fit particles it reads
    grid.set_dimension( local_dims );
    grid.set_offset( local_offset );

    //MPI_Barrier(comm);
    //auto T0 = std::chrono::high_resolution_clock::now();

    auto type_name_to_index = [&species]( const char* typeAtome ) -> int 
    {
      int nSpecies = species.size();
      int at_type = 0;
      for(at_type = 0; at_type < nSpecies ; ++at_type)
      {
        if( species[at_type].name() == typeAtome ) return at_type;
      }
      lerr << "Atom of type " << typeAtome << " is unknown, the "<<nSpecies<<" available types are :" << std::endl;
      for(at_type = 0; at_type < nSpecies ; ++at_type) lerr << "\t" << species[at_type].name()<<std::endl << std::flush;
      std::abort();
      return -1;
    };


    //MPI_Barrier(comm);

    // --------- atomes -------------
    // ------------------------------
    ldbg << "\t--> End reading atoms from MPIIO file - offset = " << file.current_offset() << "\n" << std::endl;

    //==========================================================================================
    //==========================================================================================
    //==========================================================================================

    // ---------------------------------
    // --------- molecules -------------
    std::unique_ptr<IOMolecules[]> molecules = nullptr;
    if( entete->bloc_molecules )
    {
      ldbg << "\t--> Start reading molecules from the MPIIO file - offset = " << file.current_offset() << "\n" << std::endl;

      molecules = std::make_unique<IOMolecules[]>(count);;
      file.increment_offset( molecules.get(), particle_offset_start );           //displace the offset to the beginning of the _molecules_ block assigned to this process
      file.read( molecules.get(), count);
      file.increment_offset( molecules.get(), n_particles - particle_offset_start ); //displace the offset to the end of the global _molecules_ block
          
      std::map<std::string, uint8_t> convert_name_to_id;

      //Impression des donnees moleculaires
#     if !defined(NDEBUG) && defined(STAMP_V4_VERBOSE_DEBUG)
      for (size_t i = 0;i<count; i++)
      {
        ldbg 
             << "\t\tatom id "       << std::setw( 6) << atomes[i].numeroAtome 
             << " - typeFF "         << std::setw( 6) << molecules[i].typeAtomeFF 
             << " - molecule id "    << std::setw( 6) << molecules[i].numeroMolecule 
             << " - molecule type "  << std::setw( 6) << molecules[i].typeMolecule 
             << " - connectivity ["  << std::setw( 6) << molecules[i].Connectivite[0] << " ; "
                                     << std::setw( 6) << molecules[i].Connectivite[1] << " ; "
                                     << std::setw( 6) << molecules[i].Connectivite[2] << " ; "
                                     << std::setw( 6) << molecules[i].Connectivite[3] << " ; "
                                     << std::setw( 6) << molecules[i].Connectivite[4] << "] "
             << std::endl;
      }
      ldbg << std::endl;
#     endif

      ldbg << "\t--> End reading molecules from the MPIIO file - offset = " << file.current_offset() << "\n" << std::endl << std::flush;
    }

    if( entete->bloc_dpd )
    {
      lerr << "StampV4 read error: DPD optional block not supported\n";
      std::abort();
    }

    if( entete->bloc_graines )
    {
      lerr << "StampV4 read error: 'graines' optional block not supported\n";
      std::abort();
    }

    std::unique_ptr<IOMolRig[]> rigidmols = nullptr;
    if( entete->bloc_molrig )
    {
      ldbg << "\t--> Start reading rigid molecules from the MPIIO file - offset = " << file.current_offset() << "\n" << std::endl;

      rigidmols = std::make_unique<IOMolRig[]>(count);
      file.increment_offset( rigidmols.get(), particle_offset_start );           //displace the offset to the beginning of the _molecules_ block assigned to this process
      file.read( rigidmols.get(), count);
      file.increment_offset( rigidmols.get(), n_particles - particle_offset_start ); //displace the offset to the end of the global _molecules_ block
          
      //Impression des donnees moleculaires
#     if !defined(NDEBUG) && defined(STAMP_V4_VERBOSE_DEBUG)
      for (size_t i = 0;i<count; i++)
      {
        ldbg 
             << "\t\tatom id " << std::setw(6) << atomes[i].numeroAtome 
             << " quat=("   <<std::setw(6)<< rigidmols[i].quaternion     [0] <<","<<std::setw(6)<< rigidmols[i].quaternion     [1] <<","<< std::setw(6)<< rigidmols[i].quaternion     [2]<<","<<std::setw(6)<< rigidmols[i].quaternion[3]<<")"
             << " angmom=(" <<std::setw(6)<< rigidmols[i].momentangulaire[0] <<","<<std::setw(6)<< rigidmols[i].momentangulaire[1] <<","<< std::setw(6)<< rigidmols[i].momentangulaire[2]<<")"
             //<< " orient=(" <<std::setw(6)<< rigidmols[i].orientation    [0] <<","<<std::setw(6)<< rigidmols[i].orientation    [1] <<","<< std::setw(6)<< rigidmols[i].orientation    [2]<<")"
             << std::endl;
      }
      ldbg << std::endl;
#     endif

      ldbg << "\t--> End reading rigid molecules from the MPIIO file - offset = " << file.current_offset() << "\n" << std::endl << std::flush;
    }

    //==========================================================================================
    //==========================================================================================
    //==========================================================================================

    // ------------------------------------------------
    // --------- comptage des types de molecule -------
    std::vector<int> at_mol_place( count , 0 );


    // ------------------------------------------------
    // --------- transmission des donnees -------------
    using MoleculeTupleIO = onika::soatl::FieldTuple<field::_rx, field::_ry, field::_rz, 
                                              field::_vx, field::_vy, field::_vz, 
                                              field::_id, field::_type, 
                                              field::_orient , field::_angmom,
                                              field::_idmol, field::_cmol>;
    using MoleculeTuple = decltype( grid.cells()[0][0] );
    using RigidMoleculeOpt = onika::soatl::FieldTuple< field::_orient , field::_angmom >;
    using MoleculeConOpt = onika::soatl::FieldTuple< field::_idmol, field::_cmol >;
    static constexpr bool has_mol_con_fields = GridT::template HasField<field::_idmol> ::value && GridT::template HasField<field::_cmol> ::value;
    static constexpr bool has_rigid_mol_fields = GridT::template HasField<field::_orient> ::value && GridT::template HasField<field::_angmom> ::value;
    
    std::vector<RigidMoleculeOpt> rigid_mol_opt;
    std::vector<MoleculeConOpt> mol_con_opt;
    std::vector<uint64_t> particle_optional_place;

    static constexpr double vel_conv = EXASTAMP_CONST_QUANTITY( 1.0 * m/s );
    Vec3d SumVelocity = {0.,0.,0.};
    //Gestion des blocs de donnees, mise en commun dans la grid
    for (size_t i = 0;i<count; i++)
    {
      //atom absolute ID
      int64_t at_id  = atomes[i].numeroAtome;
      if( at_id < 0 )
      {
        fatal_error() << "read_stamp_v4 : Warning: bad id ("<<at_id<<") found" << std::endl;
      }

      //atom type ID
      int at_type_i = type_name_to_index(atomes[i].typeAtome);
      if( at_type_i == -1 )
      {
        lerr << "\t\tatom id="<< std::setw( 6) << atomes[i].numeroAtome<<" index="<<i
           << " type="<< std::setw( 6) << atomes[i].typeAtome << "("<< at_type_i <<")"
           << " pos=("<< std::setw(14) << std::setprecision(6) << atomes[i].Position[0] << "," 
                      << std::setw(14) << std::setprecision(6) << atomes[i].Position[1] << "," 
                      << std::setw(14) << std::setprecision(6) << atomes[i].Position[2] <<")"
           << " vel=("<< std::setw(14) << std::setprecision(6) << atomes[i].Vitesse[0] << ","
                      << std::setw(14) << std::setprecision(6) << atomes[i].Vitesse[1] << ","
                      << std::setw(14) << std::setprecision(6) << atomes[i].Vitesse[2] << ")"
          << std::endl ;
        std::abort();
      }
      uint8_t at_type = at_type_i;

      //atom position
      Vec3d r = { atomes[i].Position[0] , atomes[i].Position[1] , atomes[i].Position[2] };
      r = r * position_scaling + position_offset;
      r = position_inv_xform * r; // account for domain aspect ratio adjustments

      //check if atom is inside the grid
      Vec3d ro = r;
      IJK loc = domain_periodic_location( domain, r ) - grid.offset();
      if( ! grid.contains(loc) )
      {
        lerr<<"Domain = "<<domain<<std::endl;
        lerr<<"Domain size = "<<domain.bounds_size()<<std::endl;
        lerr<<"particle #"<<at_id<<", ro="<<ro<<", r="<<r<<"<< in cell "<<loc<<" not in grid : offset="<<grid.offset()<<std::endl<<std::flush;
        std::abort();
      }

      if( ! is_inside(grid.cell_bounds(loc),r) )
      {
        lerr<<"particle #"<<at_id<<", ro="<<ro<<", r="<<r<<"<< in cell "<<loc<<" not inside cell bounds ="<<grid.cell_bounds(loc)<<std::endl<<std::flush;
        std::abort();
      }

      //atom velocity
      Vec3d v = {atomes[i].Vitesse[0] * vel_conv, atomes[i].Vitesse[1] * vel_conv, atomes[i].Vitesse[2] * vel_conv};
      SumVelocity += v;

      uint64_t mol_id = 0;
      uint64_t molnum = 0;
      unsigned int moltype = 0;
      //unsigned int molplace = 0;
      std::array<uint64_t,4> cmol = { 0, 0, 0, 0 };

      if( entete->bloc_molecules )
      {
        //molecular ID : note, ids are INT in stamp, not UINT.
        molnum = static_cast<uint64_t>( molecules[i].numeroMolecule);
        moltype = static_cast<unsigned int>( molecules[i].typeMolecule);
        mol_id = molnum + ( uint64_t(moltype) << 48 );

        //atom connectivity : note, ids are INT in stamp, not UINT.
        cmol = std::array<uint64_t,4>{static_cast<uint64_t>(molecules[i].Connectivite[0]),
                                      static_cast<uint64_t>(molecules[i].Connectivite[1]),
                                      static_cast<uint64_t>(molecules[i].Connectivite[2]),
                                      static_cast<uint64_t>(molecules[i].Connectivite[3])};
      }

      Quaternion orient = { 0., 0., 0., 0. };
      Vec3d angmom = { 0., 0., 0. };
      if( entete->bloc_molrig )
      {
        orient = Quaternion { rigidmols[i].quaternion[0] , rigidmols[i].quaternion[1] , rigidmols[i].quaternion[2] , rigidmols[i].quaternion[3] };
        angmom = Vec3d { rigidmols[i].momentangulaire[0] , rigidmols[i].momentangulaire[1] , rigidmols[i].momentangulaire[2] };
      }
      
      MoleculeTupleIO full_tp = MoleculeTupleIO(r.x, r.y, r.z, v.x, v.y, v.z, at_id, at_type, orient, angmom, mol_id, cmol);
      const MoleculeTuple t = full_tp;
      const size_t cell_idx = grid.cell_index(loc);
      const size_t p_idx = grid.cell(cell_idx).size();
      grid.cell(cell_idx).push_back( t );

      //particle_optional_data.push_back( tp );
      if( ( entete->bloc_molecules && !has_mol_con_fields ) || ( entete->bloc_molrig && !has_rigid_mol_fields ) )
      {
        particle_optional_place.push_back( encode_cell_particle( cell_idx , p_idx ) );
      }
      if( entete->bloc_molecules && !has_mol_con_fields )
      {
        mol_con_opt.push_back(full_tp);
      }
      if( entete->bloc_molrig && !has_rigid_mol_fields )
      {
        rigid_mol_opt.push_back(full_tp);
      }
    }

    SumVelocity /= count;
    lout << "Average velocity  = "<<SumVelocity<<std::endl;
  
    // --------- transmission des donnees -------------
    // ------------------------------------------------

    //==========================================================================================
    //==========================================================================================
    //==========================================================================================

    file.close();
    grid.rebuild_particle_offsets();
    
    // fill in optional fields if needed
    if constexpr ( ! has_mol_con_fields )
    {
      if( entete->bloc_molecules )
      {
        lout << "Allocating optional fields idmol and cmol" << std::endl;
        if( particle_optional_place.size() != mol_con_opt.size() ) fatal_error() << "Internal error : Inconsistent optional data" << std::endl;        
        auto field_idmol = grid.field_accessor( field::idmol );
        auto field_cmol = grid.field_accessor( field::cmol );
        auto cells = grid.cells_accessor();
        const size_t nb_particles = grid.number_of_particles();
        if( nb_particles != particle_optional_place.size() ) fatal_error() << "Internal error : optional data size does not match number of particles" << std::endl;        
        for(size_t i=0;i<nb_particles;i++)
        {
          field_idmol.m_flat_array_ptr[i] = std::numeric_limits<uint64_t>::max();
          field_cmol.m_flat_array_ptr[i][0] = std::numeric_limits<uint64_t>::max();
          field_cmol.m_flat_array_ptr[i][1] = std::numeric_limits<uint64_t>::max();
          field_cmol.m_flat_array_ptr[i][2] = std::numeric_limits<uint64_t>::max();
          field_cmol.m_flat_array_ptr[i][3] = std::numeric_limits<uint64_t>::max();
        }
        for(size_t i=0;i<nb_particles;i++)
        {
          size_t cell_idx=0, p_idx=0;
          decode_cell_particle( particle_optional_place[i], cell_idx , p_idx );
          cells[cell_idx][field_idmol][p_idx] = mol_con_opt[i][field::idmol];
          cells[cell_idx][field_cmol][p_idx] = mol_con_opt[i][field::cmol];
        }
      }
    }
    
    if constexpr ( ! has_rigid_mol_fields )
    {
      if( entete->bloc_molrig )
      {
        lout << "Allocating optional fields orient and angmom" << std::endl;
        if( particle_optional_place.size() != rigid_mol_opt.size() ) fatal_error() << "Internal error : Inconsistent optional data" << std::endl;        
        auto field_orient = grid.field_accessor( field::orient );
        auto field_angmom = grid.field_accessor( field::angmom );
        auto cells = grid.cells_accessor();
        const size_t nb_particles = grid.number_of_particles();
        if( nb_particles != particle_optional_place.size() ) fatal_error() << "Internal error : optional data size does not match number of particles" << std::endl;        
        for(size_t i=0;i<nb_particles;i++)
        {
          field_orient.m_flat_array_ptr[i] = Quaternion{0.,0.,0.,0.};
          field_angmom.m_flat_array_ptr[i] = Vec3d{0.,0.,0.};
        }
        for(size_t i=0;i<nb_particles;i++)
        {
          size_t cell_idx=0, p_idx=0;
          decode_cell_particle( particle_optional_place[i], cell_idx , p_idx );
          cells[cell_idx][field_orient][p_idx] = rigid_mol_opt[i][field::orient];
          cells[cell_idx][field_angmom][p_idx] = rigid_mol_opt[i][field::angmom];
        }
      }
    }

    ldbg << "*****\n** StampV4 dump file <"<< file_name <<"> read\n*****\n" << std::endl ;
    lout << "==========================================" << std::endl << std::endl;

  }


  template<class LDBG, class GridT>
  static inline void read_stamp_v4(
                       LDBG& ldbg,
                       MPI_Comm comm,
                       const std::string& file_name,
                       double enlarge_bounds,
                       ReadBoundsSelectionMode bounds_mode,
                       GridT& grid,
                       Domain& domain,
                       int64_t& iteration_number,
                       double& dt,
                       double& phystime,
                       ParticleSpecies& species,                     
                       bool enable_xform_scale,
                       bool pbc_adjust_xform,
                       int force_version ,
                       bool use_molrigid_42 )
  {
    using stampv4::IOVersion;

    // Checks that grid is empty and zero initialized
    grid.reset();

    //set some global ldbg preferences
    ldbg << std::scientific ;
    ldbg << std::showpoint ;

    // read Stamp V4 file
    ldbg << "\n\n*****\n** Reading StampV4 dump file <"<< file_name <<">\n*****\n" << std::endl ;

    // open MPIIO file
    exanb::MpiIO file;
    file.open(file_name,"r");

    //==========================================================================================
    //==========================================================================================
    //==========================================================================================

    ldbg << "\t--> Start reading the MPIIO file version - offset = " << file.current_offset() << "\n" << std::endl ;
    // -------------------------------
    // --------- version -------------
    IOVersion version = { 1 };
    file.read( &version );
    file.increment_offset( &version );
    assert( exanb::all_value_equal(comm, version ) );
    ldbg << "\t\tversion = " << version.version << "\n" << std::endl;
    if( force_version != -1 )
    {
          ldbg << "\t\tforced version = " << force_version << "\n" << std::endl;
    }
    // --------- version -------------
    // -------------------------------
    ldbg << "\t--> End reading the MPIIO file version - offset = " << file.current_offset() << "\n" << std::endl ;

    if( ( version.version == 1 && force_version==-1 ) || force_version == 41 )
    {
      using Version = std::integral_constant<unsigned int,41>;
      if( use_molrigid_42 )
      {
        read_stamp_v4(file,ldbg,comm,file_name,enlarge_bounds,bounds_mode,grid,domain,iteration_number,dt,phystime,species,enable_xform_scale, pbc_adjust_xform, StampV4Options<Version,std::true_type>{} );
      }
      else
      {
        read_stamp_v4(file,ldbg,comm,file_name,enlarge_bounds,bounds_mode,grid,domain,iteration_number,dt,phystime,species,enable_xform_scale, pbc_adjust_xform, StampV4Options<Version,std::false_type>{} );
      }
    }
    else if( ( version.version == 2 && force_version==-1 ) || force_version == 42 )
    {
      using Version = std::integral_constant<unsigned int,42>;
      if( use_molrigid_42 )
      {
        read_stamp_v4(file,ldbg,comm,file_name,enlarge_bounds,bounds_mode,grid,domain,iteration_number,dt,phystime,species,enable_xform_scale, pbc_adjust_xform, StampV4Options<Version,std::true_type>{} );
      }
      else
      {
        read_stamp_v4(file,ldbg,comm,file_name,enlarge_bounds,bounds_mode,grid,domain,iteration_number,dt,phystime,species,enable_xform_scale, pbc_adjust_xform, StampV4Options<Version,std::false_type>{} );
      }
    }
    else
    {
      lerr << "Unsupported file version number "<<version.version<<std::endl;
      std::abort();
    }
  }


}

