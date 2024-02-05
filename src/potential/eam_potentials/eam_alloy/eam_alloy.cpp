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

    // First line to be read contains the number of materials N followed by the N types (ex: "2 Al Cu")
    std::getline(file,line);
    char split_char = ' ';
    std::istringstream split(line);
    std::vector<std::string> tokens;
    for (std::string each; std::getline(split, each, split_char); tokens.push_back(each));
    int nmats = tokens.size()-1;
    std::vector<std::string> matnames;
    matnames.resize(nmats);
    for (int i=0;i<nmats;i++) {
      matnames[i] = tokens[i+1];
    }
    std::cout << "Nmats = " << nmats << std::endl;
    for (int i=0;i<nmats;i++) {
      std::cout << "\t Material #" << i << " = " << matnames[i] << std::endl;
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
    v.nmats = nmats;
    v.nr = nr;
    v.nrho = nrho;
    v.rhomax = (nrho - 1) * drho;
    v.rc = rcut;
    // When multimaterial --> frho and rhor need to be vector<vector<double>> bc 1 function per type
    std::vector<std::vector<double>> frho;
    std::vector<std::vector<double>> rhor;
    frho.resize(nmats);
    rhor.resize(nmats);
    for (int i = 0; i < nmats; i++) {
      frho[i].resize(nrho+1);
      rhor[i].resize(nr+1);
    }
    
    // Read f(rho) and r * rho(r) by block for each material
    for (int i = 0; i < nmats; i++) {
      std::getline(file,line);
      std::cout << "For material # "<<i<<" :" <<line<<std::endl;
      int zmat;
      double massmat,azeromat;
      std::string structmat;
      std::istringstream(line) >> zmat >> massmat >> azeromat >> structmat;
      std::cout << "\tzmat  = " << zmat << std::endl;
      std::cout << "\tmass  = " << massmat << std::endl;    
      std::cout << "\tazero = " << azeromat << std::endl;
      std::cout << "\tstruc = " << structmat << std::endl;    

      size_t cnt_frho=0;
      while(cnt_frho != nrho) {
        std::getline(file,line);
        std::istringstream split(line);
        std::vector<std::string> tokens;
        for (std::string each; std::getline(split, each, split_char); tokens.push_back(each));
        int nvals=tokens.size();
        for (int j = 0; j<nvals; j++) {
          frho[i][cnt_frho+j] = std::stod(tokens[j]);
        }
        cnt_frho+=nvals;
        //        std::cout << "cnt_frho = " << cnt_frho << std::endl;;
      }

      size_t cnt_rhor=0;
      while(cnt_rhor != nr) {
        std::getline(file,line);
        std::istringstream split(line);
        std::vector<std::string> tokens;
        for (std::string each; std::getline(split, each, split_char); tokens.push_back(each));
        int nvals=tokens.size();
        for (int j = 0; j<nvals; j++) {
          rhor[i][cnt_rhor+j] = std::stod(tokens[j]);
        }
        cnt_rhor+=nvals;
        //        std::cout << "cnt_rhor = " << cnt_rhor << std::endl;;
      }

      // std::cout << "frho0 start = " << frho[0][0] << std::endl;
      // std::cout << "frho0 end   = " << frho[0][nrho-1] << std::endl;

      // std::cout << "rhor0 start = " << rhor[0][0] << std::endl;
      // std::cout << "rhor0 end   = " << rhor[0][nr-1] << std::endl;      

      // std::cout << "frho1 start = " << frho[1][0] << std::endl;
      // std::cout << "frho1 end   = " << frho[1][nrho-1] << std::endl;

      // std::cout << "rhor1 start = " << rhor[1][0] << std::endl;
      // std::cout << "rhor1 end   = " << rhor[1][nr-1] << std::endl;
      
    }
    
    // Read phi(r) by block for each material
    // If N materials first read the monotype interactions by block then read the cross terms 
    // If N materials there are N(N+1)/2 pair potentials to read
    int nphis=nmats*(nmats+1)/2;
    std::vector<double> phivals;
    do
      {
        std::getline(file,line);
        std::istringstream split(line);
        for (std::string each; std::getline(split, each, split_char); phivals.push_back(std::stod(each)));
      } while( !file.eof() );

    std::vector<std::vector<double>> z2r;
    z2r.resize(nphis);    
    for (int i = 0; i<nphis; i++) {
      int start = i * nr;
      int end = start + nr;
      z2r[i].resize(nr);
      z2r[i].insert(z2r[i].begin(),phivals.begin() + start, phivals.begin() + end + 1);
    }
    std::cout << "File EAM read" << std::endl;
    
    // Now that the file is read, we need to spline the functions and store them in the EAM parameter coefs
    // First, declare the spline arrays 
    v.rdr = 1.0/dr;
    v.rdrho = 1.0/drho;
    int nspl = 7;

    // Spline array for f(rho)
    //    std::vector<std::vector<std::vector<double>>> frho_spline;
    v.frho_spline.resize(nmats);
    for (int i = 0; i<nmats; i++) {
      v.frho_spline[i].resize(nrho);
      for (int j = 0; j<nrho; j++) {
        v.frho_spline[i][j].resize(nspl);
      }
      interpolate(nrho,drho,frho[i],v.frho_spline[i]);      
    }
    
    // Spline array for rho(r)    
    //    std::vector<std::vector<std::vector<double>>> rhor_spline;
    v.rhor_spline.resize(nmats);
    for (int i = 0; i<nmats; i++) {
      v.rhor_spline[i].resize(nr);
      for (int j = 0; j<nr; j++) {
        v.rhor_spline[i][j].resize(nspl);
      }
      interpolate(nr,dr,rhor[i],v.rhor_spline[i]);
    }

    // Spline array for phi(r)
    // std::vector<std::vector<std::vector<double>>> z2r_spline;
    v.z2r_spline.resize(nphis);
    for (int i = 0; i<nphis; i++) {
      v.z2r_spline[i].resize(nr);
      for (int j = 0; j<nr; j++) {
        v.z2r_spline[i][j].resize(nspl);
      }
      interpolate(nr,dr,z2r[i],v.z2r_spline[i]);
    }

    return true;
  }
    
}

