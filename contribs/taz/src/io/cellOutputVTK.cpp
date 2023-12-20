/// @file
/// @brief Implementations for the cellWriterVTK


#include "domain/domainInfo.hpp"

#include "io/cellOutputVTK.hpp"

#include "parallel/commManager.hpp"


/// @brief Gather cell output from all nodes on one node
/// @param [in] comm Communication manager
/// @param [in] sendBuffer Cells on each node
/// @param [in] counts Number of cell to get on each node
/// @param [in] disps Position where to put the data from each node
void CellOutput::gather(CommManager* comm, CellOutput* sendBuffer, Array<int>& counts, Array<int>& disps) {

  // Gather positions
  comm->gatherV(sendBuffer->id, this->id, counts, disps);

  // Gather positions
  comm->gatherV(sendBuffer->r, this->r, counts, disps);
  
  // Gather pressures
  comm->gatherV(sendBuffer->p, this->p, counts, disps);

  // Gather temperatures
  comm->gatherV(sendBuffer->t, this->t, counts, disps);

  // Gather densities 
  comm->gatherV(sendBuffer->d, this->d, counts, disps);

}





/// @brief Create the name of the output file
/// @param [in] step Step where the output is written
/// @param [out] filename Output file name
void CellWriterVTK::setFilename(uint step, std::string& filename) {

  std::string tmp = std::to_string(step);

  while (tmp.size()<LOG_MAX_STEPS) tmp.insert(0, "0");

  filename = rootName + "_" + tmp + extension;

}


/// @brief Write the cells
/// @param [in] step Step where the output is written
void CellWriterVTK::write(uint step) {

  // Get the file name and open it
  std::string tmp;
  this->setFilename(step, tmp);

  std::ofstream out(tmp.c_str(), std::ios::out);

  // Write the header
  writeHeader(out);

  // Sort indices list
  uint size = cells->id.size();
  Array< uint > id(size);
  for (uint i = 0; i<size; i++)
    id[cells->id[i]] = i;
  
  // Write the positions
  writePositions(out,id);

  // Write the number of printed objects
  out<< "POINT_DATA " << cells->r.size() << std::endl;

  // Write the pressures
  writePressures(out,id);

  // Write the temepratures
  writeTemperatures(out,id);

  // Write the densities
  writeDensities(out,id);

  // Close the file
  out.close();

}


/// @brief Write the header in the specified file
/// @param [in] flux Output file
void CellWriterVTK::writeHeader(std::ofstream& flux) {

  flux<< "# vtk DataFile Version 2.0" << std::endl 
      << "A useless comment" << std::endl 
      << "ASCII"<< std::endl << std::endl;

}


/// @brief Write the cells positions in the specified file
/// @param [in] flux Output file
/// @param [in] id Cell indices
void CellWriterVTK::writePositions(std::ofstream& flux, const Array< uint >& id) {

  Array< vec3<double> >& r = cells->r;
  uint size = r.size();

  flux<< "DATASET STRUCTURED_GRID" << std::endl << "DIMENSIONS " << Global::domainInfo.getNumberOfCellsPerDim() << std::endl;
  flux<< "POINTS " <<size << " float"  << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << (r[ id[p] ]*10.0) << std::endl;
  
}


/// @brief Write the cells pressures in the specified file
/// @param [in] flux Output file
/// @param [in] id Cell indices
void CellWriterVTK::writePressures(std::ofstream& flux, const Array< uint >& id) {

  Array< mat3<double> >& pr = cells->p;
  uint size = pr.size();

  flux<< "SCALARS Pxx float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << pr[ id[p] ].m11 << std::endl;

  flux<< "SCALARS Pyy float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;
  
  for (uint p=0; p<size; ++p) 
    flux << pr[ id[p] ].m22 << std::endl;
  
  flux<< "SCALARS Pzz float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;
  
  for (uint p=0; p<size; ++p) 
    flux << pr[ id[p] ].m33 << std::endl;

  flux<< "SCALARS Pxy float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << 0.5*(pr[ id[p] ].m12+pr[ id[p] ].m21) << std::endl;

  flux<< "SCALARS Pyz float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << 0.5*(pr[ id[p] ].m23+pr[ id[p] ].m32) << std::endl;

  flux<< "SCALARS Pzx float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << 0.5*(pr[ id[p] ].m31+pr[ id[p] ].m13) << std::endl;

}


/// @brief Write the cells temperatures in the specified file
/// @param [in] flux Output file
/// @param [in] id Cell indices
void CellWriterVTK::writeTemperatures(std::ofstream& flux, const Array< uint >& id) {

  Array<double>& tp = cells->t;
  uint size = tp.size();

  flux<< "SCALARS Temperature float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << tp[ id[p] ] << std::endl;

}

/// @brief Write the cells densities in the specified file
/// @param [in] flux Output file
/// @param [in] id Cell indices
void CellWriterVTK::writeDensities(std::ofstream& flux, const Array< uint >& id) {

  Array<double>& ds = cells->d;
  uint size = ds.size();

  flux<< "SCALARS Density float 1" << std::endl<< "LOOKUP_TABLE default" << std::endl;

  for (uint p=0; p<size; ++p) 
    flux << ds[ id[p] ] << std::endl;

}

