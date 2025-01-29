#include "eam_alloy.h"
#include <onika/file_utils.h>
#include <onika/log.h>
#include <onika/physics/constants.h>
#include <fstream>


namespace exaStamp
{
  namespace EamAlloyTools
  {      
    void interpolate(int n, double delta, const std::vector<double>& f, std::vector< std::vector<double> > & spline)
    {
      for (int m = 1; m <= n; m++) spline[m][6] = f[m];

      spline[1][5] = spline[2][6] - spline[1][6];
      spline[2][5] = 0.5 * (spline[3][6]-spline[1][6]);
      spline[n-1][5] = 0.5 * (spline[n][6]-spline[n-2][6]);
      spline[n][5] = spline[n][6] - spline[n-1][6];

      for (int m = 3; m <= n-2; m++)
        spline[m][5] = ((spline[m-2][6]-spline[m+2][6]) +
                        8.0*(spline[m+1][6]-spline[m-1][6])) / 12.0;

      for (int m = 1; m <= n-1; m++) {
        spline[m][4] = 3.0*(spline[m+1][6]-spline[m][6]) -
          2.0*spline[m][5] - spline[m+1][5];
        spline[m][3] = spline[m][5] + spline[m+1][5] -
          2.0*(spline[m+1][6]-spline[m][6]);
      }

      spline[n][4] = 0.0;
      spline[n][3] = 0.0;

      for (int m = 1; m <= n; m++) {
        spline[m][2] = spline[m][5]/delta;
        spline[m][1] = 2.0*spline[m][4]/delta;
        spline[m][0] = 3.0*spline[m][3]/delta;
      }
      
    }

  }
}

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  bool convert<exaStamp::EamAlloyParameters>::decode(const Node& _node, exaStamp::EamAlloyParameters& v)
  {
    using exaStamp::EamAlloyTools::interpolate;
    using exanb::lout;

    Node node = _node;
    std::string file_to_load;
    if( node.IsScalar() )
    {
      file_to_load = onika::onika::data_file_path( node.as<std::string>() );
    }
    else
    {
      if( !node.IsMap() ) { return false; }
      if( node["file"] )
      {
        file_to_load = onika::onika::data_file_path( node["file"].as<std::string>() );
      }
    }

    std::string line;
    std::ifstream file( file_to_load );

    // Ignore first 3 lines (always considered as comments in input file)
    for (int i = 1; i < 4; ++i)
      {
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      }    

    // First line to be read contains the number of materials N followed by the N types (ex: "2 Al Cu" or "1 Ta")
    size_t nelements = 0;
    file >> nelements;
    std::vector<std::string> matnames;
    matnames.resize(nelements);
    for (size_t i=0;i<nelements;i++) {
      file >> matnames[i];
    }
    
    // std::getline(file,line);
    // char split_char = ' ';
    // std::istringstream split(line);
    // std::vector<std::string> tokens;
    // for (std::string each; std::getline(split, each, split_char); tokens.push_back(each));
    // int nelements = tokens.size()-1;

    // std::vector<std::string> matnames;
    // matnames.resize(nelements);
    // for (int i=0;i<nelements;i++) {
    //   matnames[i] = tokens[i+1];
    // }
    lout << "====== EAM Dynamo reader =====" << std::endl;
    lout << "Materials =";
    for (size_t i=0;i<nelements;i++) {
      lout << (i==0?" ":",") << matnames[i];
    }
    lout <<" ("<< nelements<<")" << std::endl;
    
    // Second line to be read contains information on the tabulation discretization for phi, rho and f(rho) functions
    std::getline(file,line);
    size_t nrho,nr;
    double drho, dr;
    double rcut;

    file >> nrho >> drho >> nr >> dr >> rcut;
    //    std::istringstream(line) >> nrho >> drho >> nr >> dr >> rcut;
    lout << "nrho  = " << nrho << std::endl;
    lout << "nr    = " << nr << std::endl;
#   ifndef NDEBUG
    lout << "drho  = " << drho << std::endl;    
    lout << "dr    = " << dr << std::endl;    
    lout << "rcut  = " << rcut << std::endl;    
#   endif
    v.nelements = nelements;
    v.nr = nr;
    v.nrho = nrho;
    v.rhomax = (nrho - 1) * drho;
    v.rc = rcut;
    size_t nfrho = nelements;//+1;
    size_t nrhor = nelements;
    // If N materials there are N(N+1)/2 pair potentials to read
    size_t nz2r = nelements*(nelements+1)/2;
    
#   ifndef NDEBUG
    lout << "nfrho = " << nfrho << std::endl;
    lout << "nrhor = " << nrhor << std::endl;
    lout << "nz2r  = " << nz2r  << std::endl;  
#   endif
    
    // When multimaterial --> frho and rhor need to be vector<vector<double>> bc 1 function per type
    std::vector< std::vector<double> > frho;
    std::vector< std::vector<double> > rhor;
    frho.resize(nelements);
    rhor.resize(nelements);
    for (size_t i = 0; i < nelements; i++) {
      frho[i].assign( nrho+1 , 0.0 );
      rhor[i].assign( nr+1 , 0.0 );
    }
    
    // Read f(rho) and rho(r) by block for each material
    for (size_t i = 0; i < nelements; i++) {
      int zmat;
      double massmat,azeromat;
      std::string structmat;

      // std::getline(file,line);
      //      std::istringstream(line) >> zmat >> massmat >> azeromat >> structmat;
      file >> std::ws;
      file >> zmat >> massmat >> azeromat >> structmat;

#     ifndef NDEBUG
      lout << "Material # "<<i<<" :" <<line<<std::endl;
      lout << "\tzmat  = " << zmat << std::endl;
      lout << "\tmass  = " << massmat << std::endl;    
      lout << "\tazero = " << azeromat << std::endl;
      lout << "\tstruc = " << structmat << std::endl;    
#     endif

      for(size_t cnt_frho=0 ; cnt_frho < nrho ; cnt_frho++) file >> frho[i][cnt_frho+1];
      for(size_t cnt_rhor=0 ; cnt_rhor < nr ; cnt_rhor++) file >> rhor[i][cnt_rhor+1];
    }
    
    // Read z2r(r) by block for each material
    // This is in fact directly r * phi(r) in eV * ang
    // If N materials : blocks are put this way :
    // 1-1; 1-2; 1-3; 1-4; 1-..; 2-2; 2-3; 2-4; 2-..; 3-3; 3-4; 3-..; 4-4; etc.
    std::vector< std::vector<double> > z2r;
    z2r.resize(nz2r);    
    for (size_t i = 0; i<nz2r; i++) {
      assert( i < z2r.size() );
      z2r[i].assign( nr+1 , 0.0 );
      for(size_t j=0;j<nr;j++)
      {
        assert( (j+1) < z2r[i].size() );
        file >> z2r[i][j+1];
      }
    }
    
//    lout << std::endl << std::endl;
    
    //    std::abort();
    // LAMMPS : in pair_eam.cpp, section equivalent to PairEAM::array2spline()
    // Now that the file is read, we need to spline the functions and store them in the EAM parameter coefs
    // First, declare the spline arrays 
    v.rdr = 1.0/dr;
    v.rdrho = 1.0/drho;
    static constexpr auto nspl = exaStamp::EamAlloyParameters::N_SPLINE_POINTS;

    // Spline array for f(rho)
    v.frho_spline.resize(nfrho);
    for (size_t i = 0; i<nfrho; i++) {
      v.frho_spline[i].resize(nrho+1);
      for (size_t j = 0; j<(nrho+1); j++) {
        v.frho_spline[i][j].assign( nspl , 0.0 );
      }
      interpolate(nrho,drho,frho[i],v.frho_spline[i]);
    }
    
    // Spline array for rho(r)    
    v.rhor_spline.resize(nrhor);
    for (size_t i = 0; i<nrhor; i++) {
      v.rhor_spline[i].resize(nr+1);
      for (size_t j = 0; j<(nr+1); j++) {
        v.rhor_spline[i][j].assign( nspl , 0.0 );
      }
      interpolate(nr,dr,rhor[i],v.rhor_spline[i]);
    }

    // Spline array for phi(r)
    v.z2r_spline.resize(nz2r);
    for (size_t i = 0; i<nz2r; i++) {
      v.z2r_spline[i].resize(nr+1);
      for (size_t j = 0; j<(nr+1); j++) {
        v.z2r_spline[i][j].assign( nspl , 0.0 );
      }
      interpolate(nr,dr,z2r[i],v.z2r_spline[i]);
    }


    // === flatten spline data to Cuda accessible arrays ===    

    // z2r          size = nz2r X (nr+1)
    // frho_spline  size = nfrho X (nrho+1) X nspl
    // rhor_spline  size = nrhor X (nr+1) X nspl
    // z2r_spline   size = nz2r X (nr+1) X nspl
    
    v.frho_spline_data.resize( nfrho * (nrho+1) );
    for (size_t i = 0; i<nfrho; i++) {
      for (size_t j = 0; j<(nrho+1); j++) {
        auto & coeffs = v.frho_spline_data[ i*(nrho+1) + j ].coeffs;
        for(size_t k=0;k<nspl;k++) coeffs[k] = v.frho_spline[i][j][k];
      }
    }
    
    v.rhor_spline_data.resize( nrhor * (nr+1) );
    for (size_t i = 0; i<nrhor; i++) {
      for (size_t j = 0; j<(nr+1); j++) {
        auto & coeffs = v.rhor_spline_data[ i*(nr+1) + j ].coeffs;
        for(size_t k=0;k<nspl;k++) coeffs[k] = v.rhor_spline[i][j][k];
      }
    }

    v.z2r_spline_data.resize( nz2r * (nr+1) );
    for (size_t i = 0; i<nz2r; i++) {
      for (size_t j = 0; j<(nr+1); j++) {
        auto & coeffs = v.z2r_spline_data[ i*(nr+1) + j ].coeffs;
        for(size_t k=0;k<nspl;k++) coeffs[k] = v.z2r_spline[i][j][k];
      }
    }

    lout << "==============================" << std::endl << std::endl;

    return true;
  }
    
}

