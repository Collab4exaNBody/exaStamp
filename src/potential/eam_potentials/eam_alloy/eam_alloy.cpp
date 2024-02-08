#include "eam_alloy.h"
#include <exanb/core/file_utils.h>
#include <exanb/core/log.h>
#include <exanb/core/physics_constants.h>
#include <fstream>

// Yaml conversion operators, allows to read potential parameters from config file
namespace YAML
{
  bool convert<exaStamp::EamAlloyParameters>::decode(const Node& _node, exaStamp::EamAlloyParameters& v)
  {
    using exaStamp::EamAlloyTools::interpolate;

    Node node = _node;
    std::string file_to_load;
    if( node.IsScalar() )
    {
      file_to_load = node.as<std::string>();
    }
    else
    {
      if( !node.IsMap() ) { return false; }
      if( node["file"] )
      {
        file_to_load = node["file"].as<std::string>();
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
    std::getline(file,line);
    char split_char = ' ';
    std::istringstream split(line);
    std::vector<std::string> tokens;
    for (std::string each; std::getline(split, each, split_char); tokens.push_back(each));
    int nelements = tokens.size()-1;

    std::vector<std::string> matnames;
    matnames.resize(nelements);
    for (int i=0;i<nelements;i++) {
      matnames[i] = tokens[i+1];
    }
    std::cout << "Nmaterials = " << nelements << std::endl;
    for (int i=0;i<nelements;i++) {
      std::cout << "\tMaterial #" << i << " = " << matnames[i] << std::endl;
    }
    
    // Second line to be read contains information on the tabulation discretization for phi, rho and f(rho) functions
    std::getline(file,line);
    size_t nrho,nr;
    double drho, dr;
    double rcut;
    std::istringstream(line) >> nrho >> drho >> nr >> dr >> rcut;
    std::cout << "\tnrho  = " << nrho << std::endl;
    std::cout << "\tdrho  = " << drho << std::endl;    
    std::cout << "\tnr    = " << nr << std::endl;
    std::cout << "\tdr    = " << dr << std::endl;    
    std::cout << "\trcut  = " << rcut << std::endl;    
    v.nelements = nelements;
    v.nr = nr;
    v.nrho = nrho;
    v.rhomax = (nrho - 1) * drho;
    v.rc = rcut;
    int nfrho = nelements;//+1;
    int nrhor = nelements;
    // If N materials there are N(N+1)/2 pair potentials to read
    int nz2r=nelements*(nelements+1)/2;
    
    std::cout << "nfrho = " << nfrho << std::endl;
    std::cout << "nrhor = " << nrhor << std::endl;
    std::cout << "nz2r  = " << nz2r  << std::endl;  
    
    // When multimaterial --> frho and rhor need to be vector<vector<double>> bc 1 function per type
    std::vector< std::vector<double> > frho;
    std::vector< std::vector<double> > rhor;
    frho.resize(nelements);
    rhor.resize(nelements);
    for (int i = 0; i < nelements; i++) {
      frho[i].assign( nrho+1 , 0.0 );
      rhor[i].assign( nr+1 , 0.0 );
    }
    
    // Read f(rho) and rho(r) by block for each material
    for (int i = 0; i < nelements; i++) {
      file >> std::ws;
      std::getline(file,line);
      std::cout << "For Material # "<<i<<" :" <<line<<std::endl;
      int zmat;
      double massmat,azeromat;
      std::string structmat;
      std::istringstream(line) >> zmat >> massmat >> azeromat >> structmat;
      std::cout << "\tzmat  = " << zmat << std::endl;
      std::cout << "\tmass  = " << massmat << std::endl;    
      std::cout << "\tazero = " << azeromat << std::endl;
      std::cout << "\tstruc = " << structmat << std::endl;    

      for(size_t cnt_frho=0 ; cnt_frho < nrho ; cnt_frho++) file >> frho[i][cnt_frho+1];
/*
      size_t cnt_frho=0;
      while(cnt_frho != nrho) {
        std::getline(file,line);
        std::istringstream split(line);
        std::vector<std::string> tokens;
        for (std::string each; std::getline(split, each, split_char); tokens.push_back(each));
        int nvals=tokens.size();
        for (int j = 0; j<nvals; j++) {
          frho[i][cnt_frho+j+1] = std::stod(tokens[j]);
        }
        cnt_frho+=nvals;
      }
*/      
      for(size_t cnt_rhor=0 ; cnt_rhor < nr ; cnt_rhor++) file >> rhor[i][cnt_rhor+1];
/*
      size_t cnt_rhor=0;
      while(cnt_rhor != nr) {
        std::getline(file,line);
        std::istringstream split(line);
        std::vector<std::string> tokens;
        for (std::string each; std::getline(split, each, split_char); tokens.push_back(each));
        int nvals=tokens.size();
        for (int j = 0; j<nvals; j++) {
          rhor[i][cnt_rhor+j+1] = std::stod(tokens[j]);
        }
        cnt_rhor+=nvals;
      }
*/    
    }
    
    // Read z2r(r) by block for each material
    // This is in fact directly r * phi(r) in eV * ang
    // If N materials : blocks are put this way :
    // 1-1; 1-2; 1-3; 1-4; 1-..; 2-2; 2-3; 2-4; 2-..; 3-3; 3-4; 3-..; 4-4; etc.
/*
    std::vector<double> phivals;
    do
      {
        std::getline(file,line);
        std::istringstream split(line);
        for (std::string each; std::getline(split, each, split_char); phivals.push_back(std::stod(each)));
      } while( !file.eof() );
*/

    std::vector< std::vector<double> > z2r;
    z2r.resize(nz2r);    
    for (int i = 0; i<nz2r; i++) {
      int start = i * nr;
      int end = start + nr;
      assert(i<z2r.size());
      z2r[i].assign( nr+1 , 0.0 );
      for(int j=0;j<nr;j++)
      {
        assert( (j+1) < z2r[i].size() );
        file >> z2r[i][j+1];
      }
      //z2r[i].insert(z2r[i].begin()+1,phivals.begin() + start, phivals.begin() + end + 1);
    }
    
    std::cout << "File EAM read" << std::endl;

    // LAMMPS : in pair_eam.cpp, section equivalent to PairEAM::array2spline()
    // Now that the file is read, we need to spline the functions and store them in the EAM parameter coefs
    // First, declare the spline arrays 
    v.rdr = 1.0/dr;
    v.rdrho = 1.0/drho;
    int nspl = 7;

    // Spline array for f(rho)
    v.frho_spline.resize(nfrho);
    for (int i = 0; i<nfrho; i++) {
      v.frho_spline[i].resize(nrho+1);
      for (size_t j = 0; j<(nrho+1); j++) {
        v.frho_spline[i][j].assign( nspl , 0.0 );
      }
      interpolate(nrho,drho,frho[i],v.frho_spline[i]);
    }
    
    // Spline array for rho(r)    
    v.rhor_spline.resize(nrhor);
    for (int i = 0; i<nrhor; i++) {
      v.rhor_spline[i].resize(nr+1);
      for (size_t j = 0; j<(nr+1); j++) {
        v.rhor_spline[i][j].assign( nspl , 0.0 );
      }
      interpolate(nr,dr,rhor[i],v.rhor_spline[i]);
    }

    // Spline array for phi(r)
    v.z2r_spline.resize(nz2r);
    for (int i = 0; i<nz2r; i++) {
      v.z2r_spline[i].resize(nr+1);
      for (size_t j = 0; j<(nr+1); j++) {
        v.z2r_spline[i][j].assign( nspl , 0.0 );
      }
      interpolate(nr,dr,z2r[i],v.z2r_spline[i]);
    }

    return true;
  }
    
}

