#include "domain/domainInfo.hpp"
#include "domain/domainInterface.hpp"
#include "parallel/node.hpp"

#include <fstream>
#include <string>
#include <set>
#include <assert.h>
#include <vector>
#include <algorithm>


static inline size_t check_double_array(const double* array, size_t N)
{
  for(size_t i=0; i<N; i++)
  {
    if( std::isnan(array[i]) || std::isinf(array[i]) )
    {
        return i;
    }
  }
  return N;
}

static void assert_array_not_corrupt(const double* array, size_t N, const std::string& message)
{
  size_t i = check_double_array(array,N);
  if( i < N )
  {
    std::cerr<<message<<" at position "<<i<<" / "<<N<< std::endl;
    abort();
  }
}



/// @brief local structure to store and sort particles force according to their Ids
struct ParticleIdAndForce
{
	uint64_t m_id;
	vec3<double> m_force;
	inline bool operator < (const ParticleIdAndForce& rhs) const { return m_id < rhs.m_id ; }
};


void NodeSingleDomain::verifySimulationYAML(uint64_t timeStep)
{
  YAML::Node node = YAML::LoadFile( inputOutput()->checkFile() );
  std::set<ParticleReferenceValue> reference_values_set;
  {
    std::vector<ParticleReferenceValue> reference_values;
    reference_values = node.as< std::vector<ParticleReferenceValue> >();
    reference_values_set.insert( reference_values.begin(), reference_values.end() );
  }

  const double maxErrorL2Norm = inputOutput()->checkEpsilon();

  int nproc = commManager.getNumberOfNodes();
  int rank = commManager.getRank();
  
  int64_t nbParticles = domain->getNumberOfParticles();
  
  double ae = 0.0;
  double re = 0.0;

 
  domain->ctest(reference_values_set, ae, re);


  // sum errors across all processors
  if( nproc > 1 )
  {
    double tmp[2] = { re, ae };
    MPI_Allreduce( MPI_IN_PLACE , tmp , 2 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD);
    re = tmp[0];
    ae = tmp[1];
  }

  std::cout << "Acceleration error L2 norm = "<<ae<<std::endl;
  std::cout << "Position error L2 norm = "<<re<<std::endl;
  std::cout << "error threshold = "<< maxErrorL2Norm << std::endl;
  
  if( ae > maxErrorL2Norm || re > maxErrorL2Norm )
  {
    std::abort();
  }  
}

/// @brief compare computed values with reference values stored in "check_values" file
void NodeSingleDomain::verifySimulationValues(uint64_t timeStep)
{
  if( this->inputOutput()->checkValues() && ( timeStep==this->inputOutput()->checkIteration() || this->inputOutput()->checkIteration()==-1 ) )
  {
    std::cout<<"CHECK VALUES: iteration = "<<timeStep<<std::endl<<std::flush;
        verifySimulationYAML(timeStep);
  }
}

