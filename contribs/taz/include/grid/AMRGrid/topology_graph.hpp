class topology_graph {

public :

  std::vecto<std::vector<uint>> Wave;

public : 

  void find_source(bool isRectilinear);


};

 

void find_source(RectilinearGridInfo & info)
{

  Wave.resize(8);

  vec3<int> getNumberOfCellsPerDim = info.getNumberOfCellsPerDim();

  for(int k=0 ; k<numberOfCellsPerDim.z ; k+=2)
    for(int j=0 ; j<numberOfCellsPerDim.y ; j+=2)
      for(int i=0 ; i<numberOfCellsPerDim.x ; i+=2) 
        Wave[0].push_back(info.convert(vec3<int>(i,j,k));

      
}
