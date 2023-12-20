#include <exanb/core/basic_types_stream.h>
#include <exanb/core/basic_types_yaml.h>
#include <exaStamp/compute/thermodynamic_state.h>
#include <exanb/core/domain.h>
#include <exanb/core/grid.h>
#include <exanb/core/log.h>
#include <exanb/core/make_grid_variant_operator.h>
#include <exanb/core/operator.h>
#include <exanb/core/operator_slot.h>
#include <exanb/core/operator_factory.h>
#include <exaStamp/particle_species/particle_specie.h>
#include <exanb/core/unityConverterHelper.h>

#include <exanb/io/mpi_file_io.h>
#include <exaStamp/io/stampv4_io.h>

#include <exaStamp/molecule/molecule_species.h>

#include <mpi.h>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <experimental/filesystem>

//#define XSTAMP_STAMPV4_VERBOSE_DBG 1

namespace exaStamp
{
  using namespace exanb;

  template<class GridT>
  struct WriteStampV4Operator : public OperatorNode
  {
    static constexpr size_t DEFAULT_BLOCK_SIZE = 4ull * 1024ull * 1024ull; // write maximum 4M elements at a time

    ADD_SLOT(MPI_Comm           , mpi                 , INPUT , MPI_COMM_WORLD, DocString{"MPI communicator"} );
    ADD_SLOT(std::string        , filename            , INPUT                 , DocString{"File name"} );
    ADD_SLOT(GridT              , grid                , INPUT , REQUIRED      , DocString{"Particle grid"} );
    ADD_SLOT(long               , timestep            , INPUT , REQUIRED      , DocString{"Iteration number"} );
    ADD_SLOT(double             , dt                  , INPUT , REQUIRED      , DocString{"Time step"} );
    ADD_SLOT(double             , physical_time       , INPUT , REQUIRED      , DocString{"Physical time"} );
    ADD_SLOT(ParticleSpecies    , species             , INPUT , REQUIRED      , DocString{"Particle species data block"} );
    ADD_SLOT(ThermodynamicState , thermodynamic_state , INPUT , REQUIRED      , DocString{"Current thermodynamic data"} );
    ADD_SLOT(long               , block_size          , INPUT , DEFAULT_BLOCK_SIZE , DocString{"Block size for collective writing"} );
    ADD_SLOT(long               , max_file_size       , INPUT , exanb::MpiIO::DEFAULT_MAX_FILE_SIZE , DocString{"Maximum file size"} );
    ADD_SLOT(Domain             , domain              , INPUT , REQUIRED      , DocString{"Deformation box matrix"} );
    ADD_SLOT(bool               , xsv2_compatibility  , INPUT , false         , DocString{"set to true to read dump with exaStampV2, false to read it with StampV4"} );
    ADD_SLOT(long               , version             , INPUT , 42            , DocString{"variant to use (41 or 42)"} );
    ADD_SLOT(bool               , gen_user_input      , INPUT , false         , DocString{"print user input to include in order to read file"} );

    // -----------------------------------------------
    // ----------- Operator documentation ------------
    inline std::string documentation() const override final
    {
      return R"EOF(
        Creates a protection file in StampV4 format: 
        this format is fully compatible with the Stamp4 MD code (read/write)
        )EOF";
    }
  
    inline void execute () override final
    {
      if( *version == 41 ) write_file( std::integral_constant<unsigned int,41> {} );
      else if( *version == 42 ) write_file( std::integral_constant<unsigned int,42> {} );
      else
      {
        lerr << "unspported file version "<< (*version) / 10.0 << std::endl;
        std::abort();
      }
      
      if( *gen_user_input )
      {
        lout << std::endl << std::endl << "# include the following in user input file" << std::endl << std::endl << "# Domain properties (some will be overwritten by data reader)" << std::endl;
        exanb::print_user_input( *domain , lout );
        lout << std::endl << "# Species description" << std::endl << "species:"<< std::endl ;
        exaStamp::print_user_input( *species , lout , 1 );
      }
    }

    template<class Version = std::integral_constant<unsigned int,41> >
    inline void write_file(Version = {})
    {
      using IndexT = std::conditional_t< Version::value==41 , int32_t , int64_t >;
      using stampv4::IOVersion;
      //using IOGraines = stampv4::IOGraines;
      using IOMolecules = stampv4::IOMolecules;
      using IOMonomeres = stampv4::IOMonomeres;
      using IOEntete = stampv4::IOEntete<IndexT>;
      using IOAtomes = stampv4::IOAtomes<IndexT>;
      using IOMolRig = std::conditional_t< Version::value==41 , stampv4::IOMolRigV4_1 , stampv4::IOMolRigV4_2 >;

      using std::ostringstream;

      namespace fs = std::experimental::filesystem;

      using has_field_id_t     = typename GridT:: template HasField <field::_id>;
      static constexpr bool has_field_id = has_field_id_t::value;

      using has_field_charge_t = typename GridT::CellParticles::template HasField < field::_charge > ;
      static constexpr bool has_field_charge = has_field_charge_t::value;

      using has_field_orient_t = typename GridT::CellParticles::template HasField < field::_orient > ;
      using has_field_angmom_t = typename GridT::CellParticles::template HasField < field::_angmom > ;
      static constexpr bool has_rigidmol = has_field_orient_t::value && has_field_angmom_t::value;

      using has_field_idmol_t = typename GridT::CellParticles::template HasField < field::_idmol > ;
      using has_field_cmol_t = typename GridT::CellParticles::template HasField < field::_cmol > ;
      static constexpr bool has_molecule = has_field_idmol_t::value && has_field_cmol_t::value;

      using ParticleTupleIO = onika::soatl::FieldTuple<field::_rx,field::_ry,field::_rz, 
                                                field::_vx,field::_vy,field::_vz,
                                                field::_fx,field::_fy,field::_fz,
                                                field::_virial,
                                                field::_id,
                                                field::_orient,
                                                field::_angmom,
                                                field::_idmol,
                                                field::_cmol,
                                                field::_orient,
                                                field::_angmom,
                                                field::_type>;
      
      //static const double coord_conv = UnityConverterHelper::convert(1.0, "m"); // conversion factor between stampv4 unit system and ExaStampV2 internal units
      static const double vel_conv = UnityConverterHelper::convert(1.0, "m/s");
      static const double coord_conv = UnityConverterHelper::convert(1.0, "m");
      static const double time_conv = UnityConverterHelper::convert(1.0, "s");
      static const double e_conv = UnityConverterHelper::convert(1.0, "J");

      ldbg << "vel_conv="<<vel_conv<<std::endl;
      ldbg << "coord_conv="<<coord_conv<<std::endl;
      ldbg << "time_conv="<<time_conv<<std::endl;
      ldbg << "e_conv="<<e_conv<<std::endl;

      //set some global ldbg preferences
      ldbg << std::scientific ;
      ldbg << std::showpoint ;

      IJK dims = grid->dimension();
      ssize_t gl = grid->ghost_layers();
      IJK gstart { gl, gl, gl };
      IJK gend = dims - IJK{ gl, gl, gl };
      auto cells = grid->cells();
      Mat3d xform = domain->xform();

      //if( ! is_diagonal(xform) )
      //{
      //  lerr<<"StampV4 writer does not support non diagonal transform matrix"<<std::endl;
      //  std::abort();
      //}

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
      size_t particle_offset_start   = 0;
      size_t particle_offset_end     = 0;
      for(int i=0;i<np;i++)
      {
        all_number_of_particles += proc_particles[i];
        if( i <  rank ) { particle_offset_start += proc_particles[i]; }
        if( i <= rank ) { particle_offset_end   += proc_particles[i]; }
      }
 
      ldbg << "\n*****\n** Writing StampV4 dump file <" << (*filename) << ">\n*****" << "\n" << std::endl ;

      // open MPIIO file
      exanb::MpiIO file;
      // file.set_compression( *compression );
      file.open( (*filename) ,"w",*max_file_size) ;

      //==========================================================================================
      //==========================================================================================
      //==========================================================================================
 
      ldbg << "\t-> Start writing version section - offset = " << file.current_offset() << "\n" << std::endl ;
      // -------------------------------
      // --------- version -------------
      IOVersion version { Version::value==41 ? 1 : 2 };
      ldbg << "\t\tversion = " << version.version << "\n" << std::endl;
      if( rank == 0 ) file.write( &version );
      file.increment_offset( &version );
      // --------- version -------------
      // -------------------------------
      ldbg << "\t-> End writing version section - offset = " << file.current_offset() << "\n" << std::endl ;

      //==========================================================================================
      //==========================================================================================
      //==========================================================================================

      ldbg << "\t-> Start writing entete section - offset = " << file.current_offset() << "\n" << std::endl ;
      // -------------------------------
      // --------- entete -------------
      IOEntete entete;
      std::memset( &entete , 0 , sizeof(IOEntete) );

      /* donnees structurales */
      entete.Natomes = all_number_of_particles;
      Vec3d a = xform * Vec3d{domain->extent().x - domain->origin().x,0.,0.} ;
      Vec3d b = xform * Vec3d{0.,domain->extent().y - domain->origin().y,0.} ;
      Vec3d c = xform * Vec3d{0.,0.,domain->extent().z - domain->origin().z} ;
      double A = norm(a) ;
      double B = norm(b) ;
      double C = norm(c) ;

      //on convertit les distances en m
      entete.long_a          = A / UnityConverterHelper::convert(1,"m") ;
      entete.long_b          = B / UnityConverterHelper::convert(1,"m") ;
      entete.long_c          = C / UnityConverterHelper::convert(1,"m") ;

      //on convertit les angles en degres
      entete.angle_a         = acos(dot(b,c)/(B*C))/acos(-1.)*180. ;
      entete.angle_b         = acos(dot(c,a)/(C*A))/acos(-1.)*180. ;
      entete.angle_g         = acos(dot(a,b)/(A*B))/acos(-1.)*180. ;

      const bool xsv2_mode = *xsv2_compatibility;
      if( xsv2_mode )
      {
        entete.set_xsv2_extension(); // version majeure exaStamp + 1000
        entete.MatriceCR[0][0] = xform.m11;
        entete.MatriceCR[1][0] = xform.m21;
        entete.MatriceCR[2][0] = xform.m31;
        entete.MatriceCR[0][1] = xform.m12;
        entete.MatriceCR[1][1] = xform.m22;
        entete.MatriceCR[2][1] = xform.m32;
        entete.MatriceCR[0][2] = xform.m13;
        entete.MatriceCR[1][2] = xform.m23;
        entete.MatriceCR[2][2] = xform.m33;
	      entete.bloc_Rjohndoe01 = domain->bounds().bmin.x;
	      entete.bloc_Rjohndoe02 = domain->bounds().bmin.y;
	      entete.bloc_Rjohndoe03 = domain->bounds().bmin.z;
	      entete.bloc_Rjohndoe04 = domain->bounds().bmax.x;
	      entete.bloc_Rjohndoe05 = domain->bounds().bmax.y;
	      entete.bloc_Rjohndoe06 = domain->bounds().bmax.z;
	      entete.bloc_Rjohndoe07 = domain->cell_size();
      }
      else
      {
        entete.MatriceCR[0][0] = a.x/A ;
        entete.MatriceCR[1][0] = a.y/A ;
        entete.MatriceCR[2][0] = a.z/A ;
        entete.MatriceCR[0][1] = b.x/B ;
        entete.MatriceCR[1][1] = b.y/B ;
        entete.MatriceCR[2][1] = b.z/B ;
        entete.MatriceCR[0][2] = c.x/C ;
        entete.MatriceCR[1][2] = c.y/C ;
        entete.MatriceCR[2][2] = c.z/C ;
      }

      entete.XCGeo=entete.YCGeo=entete.ZCGeo = 0.;

      /* donnees temporelles */
      entete.NumeroIterationAbsolu = *timestep ;
      entete.tempsPhysique         = (*physical_time) / UnityConverterHelper::convert(1,"s") ;
      entete.dt_adaptatif          = *dt / UnityConverterHelper::convert(1,"s") ;

      /* donnees energetiques */
      /* CONSTANTES A HOMOGENEISER AVEC LE CODE */
      entete.EnergieTotale        = thermodynamic_state->total_energy()        * 1.03641882007443324881e-02 * 1.602176565e-19 ;
      entete.EnergiePotentielle   = thermodynamic_state->potential_energy()    * 1.03641882007443324881e-02 * 1.602176565e-19 ;
      entete.EnergieCinetique     = thermodynamic_state->kinetic_energy_scal() * 1.03641882007443324881e-02 * 1.602176565e-19 ;
      entete.EnergieRotationnelle = 0.                                         * 1.03641882007443324881e-02 * 1.602176565e-19 ;

      /* invariants */
      entete.invariant = 0.;

      /* variables dynamiques de simulation */
      entete.LNVhug_T            = -1. ;
      entete.LNVhug_Tref         = -1. ;
      entete.NVT_gamma       [0] = entete.NVT_gamma [1] = entete.NVT_gamma  [2] = 0. ;
      entete.NVT_gammap      [0] = entete.NVT_gammap[1] = entete.NVT_gammap [2] = 0. ;
      entete.NPT_qsi         [0] = entete.NPT_qsi   [1] = entete.NPT_qsi    [2] = 0. ;
      entete.NPT_qsip        [0] = entete.NPT_qsip  [1] = entete.NPT_qsip   [2] = 0. ;
      entete.NPH_omega       [0] = entete.NPH_omega [1] = entete.NPH_omega  [2] = 0. ;
      entete.NPH_omegap      [0] = entete.NPH_omegap[1] = entete.NPH_omegap [2] = 0. ;
      entete.NPH_pi          [0] = entete.NPH_pi    [1] = entete.NPH_pi     [2] = 0. ;
      entete.bloc_molecules      = has_molecule ;
      entete.bloc_dpd            = 0 ;
      entete.bloc_graines        = 0 ;
      entete.bloc_molrig         = has_rigidmol ;
      entete.bloc_monomeres      = 0 ;
      entete.bloc_polymerisation = 0 ;
      entete.bloc_posfiltre      = 0 ;

      if( rank == 0 ) file.write( &entete );
      file.increment_offset( &entete );

      ldbg << "\t\tnumber of atoms     = " << std::setw(12) << entete.Natomes               << "\n";
      ldbg << "\t\t*****\n";
      ldbg << "\t\tmaille length x     = " << std::setw(14) << entete.long_a                << "\n";
      ldbg << "\t\tmaille length y     = " << std::setw(14) << entete.long_b                << "\n";
      ldbg << "\t\tmaille length z     = " << std::setw(14) << entete.long_c                << "\n";
      ldbg << "\t\tmaille angle x      = " << std::setw(14) << entete.angle_a               << "\n";
      ldbg << "\t\tmaille angle y      = " << std::setw(14) << entete.angle_b               << "\n";
      ldbg << "\t\tmaille angle z      = " << std::setw(14) << entete.angle_g               << "\n";
      ldbg << "\t\t*****\n";
      ldbg << "\t\t                    [ " << std::setw(14) << entete.MatriceCR[0][0]       << " " 
                                           << std::setw(14) << entete.MatriceCR[0][1]       << " " 
                                           << std::setw(14) << entete.MatriceCR[0][2]       << "\n";
      ldbg << "\t\tmatrice             = " << std::setw(14) << entete.MatriceCR[1][0]       << " " 
                                           << std::setw(14) << entete.MatriceCR[1][1]       << " " 
                                           << std::setw(14) << entete.MatriceCR[1][2]       << "\n";
      ldbg << "\t\t                      " << std::setw(14) << entete.MatriceCR[2][0]       << " " 
                                           << std::setw(14) << entete.MatriceCR[2][1]       << " " 
                                           << std::setw(14) << entete.MatriceCR[2][2]       << "]\n";
      ldbg << "\t\t*****\n";
      ldbg << "\t\titeration number    = " << std::setw(12) << entete.NumeroIterationAbsolu << "\n";
      ldbg << "\t\tphysical time       = " << std::setw(14) << entete.tempsPhysique         << "\n";
      ldbg << "\t\tcurrent time step   = " << std::setw(14) << entete.dt_adaptatif          << "\n";
      ldbg << "\t\t*****\n";
      ldbg << "\t\tenergie totale      = " << std::setw(14) << entete.EnergieTotale         << "\n";
      ldbg << "\t\tenergie potentielle = " << std::setw(14) << entete.EnergiePotentielle    << "\n";
      ldbg << "\t\tenergie cinetique   = " << std::setw(14) << entete.EnergieCinetique      << "\n";
      ldbg << "\t\t*****\n";
      ldbg << "\t\tbloc molecule       = " << std::setw( 6) << entete.bloc_molecules        << "\n";
      ldbg << "\t\tbloc DPD            = " << std::setw( 6) << entete.bloc_dpd              << "\n";
      ldbg << "\t\tbloc graines        = " << std::setw( 6) << entete.bloc_graines          << "\n";
      ldbg << "\t\tbloc molrig         = " << std::setw( 6) << entete.bloc_molrig           << "\n";
      ldbg << "\t\tbloc monomeres      = " << std::setw( 6) << entete.bloc_monomeres        << "\n";
      ldbg << "\t\tbloc polymerisation = " << std::setw( 6) << entete.bloc_polymerisation   << "\n";
      ldbg << "\t\tbloc posfiltre      = " << std::setw( 6) << entete.bloc_posfiltre        << "\n";
      ldbg << std::endl;

      // --------- entete -------------
      // -------------------------------
      ldbg << "\t-> End writing entete section - offset = " << file.current_offset() << "\n" << std::endl ;

      //==========================================================================================
      //==========================================================================================
      //==========================================================================================

      ldbg << "\t-> Start writing atoms section - offset = " << file.current_offset() << "\n" << std::endl ;
      // ------------------------------
      // --------- atomes -------------
      IOAtomes * atomes = new IOAtomes[*block_size]; std::memset(atomes,0,sizeof(IOAtomes)*(*block_size));
      unsigned long long pcounter = 0;
      unsigned long long tcounter = 0;
      size_t id_counter = particle_offset_start + 1;
      file.increment_offset( atomes, particle_offset_start ); //displace the offset to the beginning of the specific location assigned to this process
      
      GRID_BLOCK_FOR_BEGIN(dims,gstart,gend,i,_)
      {
        size_t n = cells[i].size();
        const auto* __restrict__ cell_charges = cells[i].field_pointer_or_null(field::charge);
        for(size_t j=0;j<n;j++)
        {
          ParticleTupleIO pt = cells[i][j];

          int type_id = pt[field::type];

          strncpy( atomes[pcounter].typeAtome , (*species)[type_id].m_name, 16 ) ;

          if( has_field_charge ) { atomes[pcounter].charge = cell_charges[j] ; }
          else                   { atomes[pcounter].charge = (*species)[type_id].m_charge ; }

          if( has_field_id )     { atomes[pcounter].numeroAtome = pt[field::id] ; }
          else                   { atomes[pcounter].numeroAtome = id_counter++ ; }
          
          Vec3d r = { pt[field::rx] , pt[field::ry] , pt[field::rz] };
          if( xsv2_mode )
          {
            atomes[pcounter].Position[0] = r.x; //atomic positions should be in "STAMP-compatible" reduced coordinates
            atomes[pcounter].Position[1] = r.y; //
            atomes[pcounter].Position[2] = r.z; //
          }
          else
          {
            //r = xform * r;
            atomes[pcounter].Position[0] = ( r.x - domain->origin().x ) / ( domain->extent().x - domain->origin().x ) - 0.5; //atomic positions should be in "STAMP-compatible" reduced coordinates
            atomes[pcounter].Position[1] = ( r.y - domain->origin().y ) / ( domain->extent().y - domain->origin().y ) - 0.5;
            atomes[pcounter].Position[2] = ( r.z - domain->origin().z ) / ( domain->extent().z - domain->origin().z ) - 0.5;
            if( fabs(atomes[pcounter].Position[0]) > 0.5 || fabs(atomes[pcounter].Position[1]) > 0.5 || fabs(atomes[pcounter].Position[2]) > 0.5 )
            {
              lerr<<"atom #"<<atomes[pcounter].numeroAtome<<" : r = "<<atomes[pcounter].Position[0]<<" , "<<atomes[pcounter].Position[1]<<" , "<<atomes[pcounter].Position[2]<< std::endl;
              lerr<<"Particule coordinates not in [-0.5;0.5] as expected, StampV4 dump may be unreadable"<<std::endl;
            }
          }

	        Vec3d v = Vec3d{pt[field::vx],pt[field::vy],pt[field::vz]} / vel_conv;
          atomes[pcounter].Vitesse[0] = v.x ;
          atomes[pcounter].Vitesse[1] = v.y ;
          atomes[pcounter].Vitesse[2] = v.z ;

#         ifdef XSTAMP_STAMPV4_VERBOSE_DBG
          ldbg 
                    << "\t\tatom[" << std::setw( 6) << pcounter << "]"
                    << " - id "    << std::setw( 6) << atomes[pcounter].numeroAtome
                    << " - type "  << std::setw( 6) << atomes[pcounter].typeAtome   << " [" << std::setw(6) << type_id << "]"
                    << " - pos "   << std::setw(14) << std::setprecision(6) << atomes[pcounter].Position[0] << " " 
                                   << std::setw(14) << std::setprecision(6) << atomes[pcounter].Position[1] << " " 
                                   << std::setw(14) << std::setprecision(6) << atomes[pcounter].Position[2]
                    << " - vel "   << std::setw(14) << std::setprecision(6) << atomes[pcounter].Vitesse [0] << " " 
                                   << std::setw(14) << std::setprecision(6) << atomes[pcounter].Vitesse [1] << " " 
                                   << std::setw(14) << std::setprecision(6) << atomes[pcounter].Vitesse [2]
                    << std::endl ;
#         endif

          //mise a jour du nombre d'atomes stockes dans le bloc _atomes_
          ++ pcounter;

          //impression du bloc _atomes_ lorsque la taille de bloc max est atteinte
          if( pcounter >= size_t(*block_size) )
          {
            file.write( atomes, pcounter );
            file.increment_offset( atomes, pcounter ); //displace the offset to the beginning of the next block
            tcounter += pcounter;
            pcounter = 0;
          }
        }
      }
      GRID_BLOCK_FOR_END

      //impression du bloc _atomes_ residuel
      if( pcounter > 0 )
      {
        file.write( atomes, pcounter );
        file.increment_offset( atomes, pcounter ); //displace the offset to the end of the specific location assigned to this process;
        tcounter += pcounter;
        pcounter = 0;
      }
      ldbg << std::endl;

      //recalage de l'offset d'ecriture en fin de bloc _atomes_
      file.increment_offset( atomes , all_number_of_particles - particle_offset_end );

      MPI_Allreduce(MPI_IN_PLACE,&tcounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,comm);

      //test sur le nombre d'atomes imprimes dans le bloc _atomes_
      if( tcounter != static_cast<size_t>(entete.Natomes) )
      {
        lerr << "Wrong number of particles in output _atomes_: expected " 
             << entete.Natomes 
             << " but " 
             << tcounter 
             << " has been written" 
             << std::endl;

      }
      // --------- atomes -------------
      // ------------------------------
      ldbg << "\t-> End writing atomes section - offset = " << file.current_offset() << "\n" << std::endl ;
      delete[]  atomes;

      //==========================================================================================
      //==========================================================================================
      //==========================================================================================

      if constexpr ( has_molecule )
      {
        ldbg << "\t-> Start writing molecules section - offset = " << file.current_offset() << "\n" << std::endl ;
        // ---------------------------------
        // --------- molecules -------------
        IOMolecules * molecules = new IOMolecules[*block_size];
        pcounter = 0;
        tcounter = 0;
        file.increment_offset( molecules, particle_offset_start ); //displace the offset to the beginning of the specific location assigned to this process
        GRID_BLOCK_FOR_BEGIN(dims,gstart,gend,i,_)
        {
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            ParticleTupleIO pt = cells[i][j];
	          const auto& sp = species->at(pt[field::type]);
	          auto molid = pt[field::idmol];
            strncpy( molecules[pcounter].typeAtomeFF, sp.m_name, 16 );
            molecules[pcounter].typeMolecule    = molecule_type_from_id( molid );
            molecules[pcounter].numeroMolecule  = molecule_instance_from_id( molid );
            molecules[pcounter].Connectivite[0] = pt[field::cmol][0] ;
            molecules[pcounter].Connectivite[1] = pt[field::cmol][1] ;
            molecules[pcounter].Connectivite[2] = pt[field::cmol][2] ;
            molecules[pcounter].Connectivite[3] = pt[field::cmol][3] ;
            molecules[pcounter].Connectivite[4] = -1 ; //la connectivite est de dimension 4 dans exaStamp, 5 dans Stamp4 : on fournit -1 par defaut

#           ifdef XSTAMP_STAMPV4_VERBOSE_DBG
            //Impression des donnees moleculaires
            int at_id=0;
            if( has_field_id ) {at_id = pt[field::id] ;}
            else               {at_id = id_counter++  ;}
            ldbg 
                      << "\t\tatom["         << std::setw(6) << pcounter << "]"
                      << " - id "            << std::setw(6) << at_id
                      << " - molecule id "   << std::setw(6) << molecules[pcounter].numeroMolecule
                      << " - molecule type " << std::setw(6) << molecules[pcounter].typeMolecule
                      << " - connectivity [" << std::setw(6) << molecules[pcounter].Connectivite[0] << " ; "
                                             << std::setw(6) << molecules[pcounter].Connectivite[1] << " ; "
                                             << std::setw(6) << molecules[pcounter].Connectivite[2] << " ; "
                                             << std::setw(6) << molecules[pcounter].Connectivite[3] << " ; "
                                             << std::setw(6) << molecules[pcounter].Connectivite[4] << "]"
                      << std::endl;
  #         endif

            //mise a jour du nombre d'atomes stockes dans le bloc _molecules_
            ++ pcounter;

            //impression du bloc _molecules_ lorsque la taille de bloc max est atteinte
            if( pcounter >= size_t(*block_size) )
            {
              file.write( molecules, pcounter );
              file.increment_offset( molecules, pcounter ); //displace the offset to the beginning of the next block
              tcounter += pcounter;
              pcounter = 0;
            }
          }
        }
        GRID_BLOCK_FOR_END

        //impression du bloc _molecules_ residuel
        if( pcounter > 0 )
        {
          file.write( molecules, pcounter );
          file.increment_offset( molecules, pcounter ); //displace the offset to the end of the specific location assigned to this process;
          tcounter += pcounter;
          pcounter = 0;
        }
        ldbg << std::endl;

        //recalage de l'offset d'ecriture en fin de bloc _molecules_
        file.increment_offset( molecules, all_number_of_particles - particle_offset_end );

        MPI_Allreduce(MPI_IN_PLACE,&tcounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,comm);

        //test sur le nombre d'atomes imprimes dans le bloc _molecules_
        if( tcounter != static_cast<size_t>(entete.Natomes) )
        {
          lerr << "Wrong number of particles in output _molecules_: expected " 
               << entete.Natomes 
               << " but " 
               << tcounter 
               << " has been written" 
               << std::endl;
          std::abort();
        }
        // --------- molecules -------------
        // ---------------------------------
        ldbg << "\t-> End writing molecules section - offset = " << file.current_offset() << "\n" << std::endl ;

        delete[]  molecules;
      }

      //==========================================================================================
      //==========================================================================================
      //==========================================================================================

      if constexpr (has_rigidmol)
      {
        
        ldbg << "\t-> Start writing rigid molecules section - offset = " << file.current_offset() << "\n" << std::endl ;
        // ---------------------------------
        // --------- monomeres -------------
        IOMolRig * rmol = new IOMolRig[*block_size];
        pcounter = 0;
        tcounter = 0;
        file.increment_offset( rmol, particle_offset_start ); //displace the offset to the beginning of the specific location assigned to this process
        GRID_BLOCK_FOR_BEGIN(dims,gstart,gend,i,_)
        {
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            ParticleTupleIO pt = cells[i][j];
            auto quat = pt[field::orient];
            auto angmom = pt[field::angmom];
            rmol[pcounter].quaternion[0] = quat.w;
            rmol[pcounter].quaternion[1] = quat.x;
            rmol[pcounter].quaternion[2] = quat.y;
            rmol[pcounter].quaternion[3] = quat.z;
            rmol[pcounter].momentangulaire[0]  = angmom.x;
            rmol[pcounter].momentangulaire[1]  = angmom.y;
            rmol[pcounter].momentangulaire[2]  = angmom.z;
            
            //mise a jour du nombre d'atomes stockes dans le bloc _monomeres_
            ++ pcounter;

            //impression du bloc _monomeres_ lorsque la taille de bloc max est atteinte
            if( pcounter >= size_t(*block_size) )
            {
              file.write( rmol, pcounter );
              file.increment_offset( rmol, pcounter ); //displace the offset to the beginning of the next block
              tcounter += pcounter;
              pcounter = 0;
            }
          }
        }
        GRID_BLOCK_FOR_END
        
        //ecriture du bloc _rigidmol_ residuel
        if( pcounter > 0 )
        {
          file.write( rmol, pcounter );
          file.increment_offset( rmol, pcounter ); //displace the offset to the beginning of the next block
          tcounter += pcounter;
          pcounter = 0;
        }

        //recalage de l'offset d'ecriture en fin de bloc _rigidmol_
        file.increment_offset( rmol, all_number_of_particles - particle_offset_end );
        MPI_Allreduce(MPI_IN_PLACE,&tcounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,comm);

        //test sur le nombre d'atomes imprimes dans le bloc _rigidmol_
        if( tcounter != static_cast<size_t>(entete.Natomes) )
        {
          lerr << "Wrong number of particles in output _rigidmol_: expected " 
               << entete.Natomes 
               << " but " 
               << tcounter 
               << " has been written" 
               << std::endl;
          std::abort();
        }
        ldbg << "\t-> End writing rigidmol section - offset = " << file.current_offset() << "\n" << std::endl ;
        
        delete [] rmol;
      }

      //####################################################
      //##  this MONOMER section is useless to exaStamp  ##
      //##  BUT it is necessary to provide compatibility  ##
      //##  with Stamp4.                                  ##
      //####################################################
      if( entete.bloc_monomeres )
      {
        ldbg << "\t-> Start writing monomers section - offset = " << file.current_offset() << "\n" << std::endl ;
        // ---------------------------------
        // --------- monomeres -------------
        IOMonomeres * monomeres = new IOMonomeres[*block_size];
        pcounter = 0;
        tcounter = 0;
        file.increment_offset( monomeres, particle_offset_start ); //displace the offset to the beginning of the specific location assigned to this process
        GRID_BLOCK_FOR_BEGIN(dims,gstart,gend,i,_)
        {
  #       ifndef NDEBUG
	  auto* __restrict__ cell_ids = cells[i].field_pointer_or_null(field::id); if(cell_ids==nullptr){}
  #       endif
          size_t n = cells[i].size();
          for(size_t j=0;j<n;j++)
          {
            // ParticleTupleIO pt = cells[i][j];

            monomeres[pcounter].reactivite      = -1 ;
            monomeres[pcounter].numeroMonomere  = -1 ;

#         ifdef XSTAMP_STAMPV4_VERBOSE_DBG
            //Impression des donnees moleculaires
            int at_id=0;
            if( has_field_id ) {at_id = cell_ids[j] ;}
            else               {at_id = id_counter++  ;}
            ldbg 
                      << "\t\tatom["       << std::setw( 6) << pcounter << "]"
                      << " - id "          << std::setw( 6) << at_id
                      << " - monomere id " << std::setw( 6) << monomeres[pcounter].numeroMonomere
                      << " - reactivity "  << std::setw( 6) << monomeres[pcounter].reactivite
                      << std::endl;
  #         endif

            //mise a jour du nombre d'atomes stockes dans le bloc _monomeres_
            ++ pcounter;

            //impression du bloc _monomeres_ lorsque la taille de bloc max est atteinte
            if( pcounter >= size_t(*block_size) )
            {
              file.write( monomeres, pcounter );
              file.increment_offset( monomeres, pcounter ); //displace the offset to the beginning of the next block
              tcounter += pcounter;
              pcounter = 0;
            }
          }
        }
        GRID_BLOCK_FOR_END

        //impression du bloc _monomeres_ residuel
        if( pcounter > 0 )
        {
          file.write( monomeres, pcounter );
          file.increment_offset( monomeres, pcounter ); //displace the offset to the end of the specific location assigned to this process;
          tcounter += pcounter;
          pcounter = 0;
        }
        ldbg << std::endl;

        //recalage de l'offset d'ecriture en fin de bloc _monomeres_
        file.increment_offset( monomeres, all_number_of_particles - particle_offset_end );

        MPI_Allreduce(MPI_IN_PLACE,&tcounter,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,comm);

        //test sur le nombre d'atomes imprimes dans le bloc _monomeres_
        if( tcounter != static_cast<size_t>(entete.Natomes) )
        {
          lerr << "Wrong number of particles in output _monomeres_: expected " 
               << entete.Natomes 
               << " but " 
               << tcounter 
               << " has been written" 
               << std::endl;
          std::abort();
        }
        // --------- monomeres -------------
        // ---------------------------------
        ldbg << "\t-> End writing monomers section - offset = " << file.current_offset() << "\n" << std::endl ;
       
        delete[]  monomeres;
      }
      //==========================================================================================
      //==========================================================================================
      //==========================================================================================

      ldbg << "*****\n** StampV4 dump file <" << (*filename) << "> written\n*****\n" << std::endl ;

      //close MPIIO file
      file.close();
    }

    //==========================================================================================
    //==========================================================================================
    //==========================================================================================

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
  CONSTRUCTOR_FUNCTION
  {
    OperatorNodeFactory::instance()->register_factory( "write_stamp_v4", make_grid_variant_operator< WriteStampV4Operator > );
  }

}
